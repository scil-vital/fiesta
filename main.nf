#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count",
                "device":"$params.device",
                "registration_speed":"$params.registration_speed"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "FINTA Flow"
log.info "==============================================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Model]"
log.info ""
log.info "Model file: $params.model"
log.info ""
log.info "[Atlas]"
log.info ""
log.info "Atlas Config: $params.atlas_config"
log.info "Atlas Anat: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info "Atlas Thresholds: $params.atlas_thresholds"
log.info ""
log.info "[Processing]"
log.info ""
log.info "Device: $params.device"
log.info "Registration speed: $params.registration_speed"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

log.info "Input: $params.input"
root = file(params.input)

Channel
    .fromFilePairs("$root/**/*.trk", size: -1) { it.parent.name }
    .into{ tractogram; tractogram_for_check } // [sid, tractogram.trk]

Channel
    .fromPath("$root/**/*t1.nii.gz")
    .map{[it.parent.name, it]}
    .into{ reference; reference_for_display_binta; reference_for_display_gesta; reference_for_display_binta_gesta; reference_for_check } // [sid, t1.nii.gz]

Channel
    .fromPath("$root/**/*wm.nii.gz")
    .map{[it.parent.name, it]}
    .into{ wm; wm_for_check } // [sid, wm.nii.gz]

Channel
    .fromPath("$root/**/*peaks.nii.gz")
    .map{[it.parent.name, it]}
    .into{ peaks; peaks_for_check } // [sid, peaks.nii.gz]

if (!(params.atlas_anat) || !(params.atlas_config) || !(params.atlas_directory) || !(params.atlas_thresholds)) {
    error "You must specify all 4 atlas related input. --atlas_anat, " +
    "--atlas_config, --atlas_directory and atlas_thresholds all are mandatory."
}

if (!params.device.equals("cpu") && !params.device.equals("cuda")){
    error "Device must either be cpu or cuda. "
}

if (params.registration_speed == 1){
    registration_script = Channel.value("antsRegistrationSyNQuick.sh")
}
else if (params.registration_speed == 0){
    registration_script = Channel.value("antsRegistrationSyN.sh")
}
else {
    error "Registration speed must be 0 or 1"
}

atlas_anat = Channel.fromPath("$params.atlas_anat")

Channel.fromPath("$params.atlas_config")
    .into{atlas_config_for_binta; atlas_config_for_concatenation; atlas_config_for_gesta; atlas_config_for_finta_gesta}

Channel.fromPath("$params.atlas_directory")
    .into{atlas_directory_for_binta; atlas_directory_for_gesta}
Channel.fromPath("$params.model")
    .into{model_for_binta; model_for_gesta}
Channel.fromPath("$params.atlas_thresholds")
    .into{atlas_thresholds_for_binta; atlas_thresholds_for_gesta}
Channel.value("$params.device")
    .into{device_for_binta; device_for_gesta}

tractogram_for_check
    .join(reference_for_check)
    .join(wm_for_check)
    .join(peaks_for_check)
    .set{compatibility_check}


process Check_Files_Compatibility {
    errorStrategy 'ignore'

    input:
    set sid, 
        file(tractogram), 
        file(reference),
        file(wm),
        file(peaks) from compatibility_check // [sid, tractogram.trk, t1.nii.gz, wm.nii.gz, peaks.nii.gz]

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, atlas.nii.gz, t1.nii.gz]
    set sid into sid_registration, sid_apply_registration

    script:
    """
    compatibility=\$(scil_verify_space_attributes_compatibility.py ${tractogram} ${reference} ${wm} ${peaks})
    if [[ \$compatibility != "All input files have compatible headers." ]]
    then
        exit 1
    fi    
    """
}

sid_registration
    .join(reference)
    .combine(atlas_anat)
    .combine(registration_script)
    .set{reference_atlas_anat}


process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(reference), file(atlas_anat), val(registration_script) from reference_atlas_anat // [sid, t1.nii.gz, atlas.nii.gz]

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, atlas.nii.gz, t1.nii.gz]
    set sid, "${sid}__output0GenericAffine.mat", 
        "${sid}__output1InverseWarp.nii.gz", 
        "${atlas_anat}", 
        "${sid}__native_anat.nii.gz" into transformation_for_tractogram, transformation_for_gesta
    set sid, "${sid}__output1Warp.nii.gz" into warp_for_gesta
    file "${sid}__outputWarped.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    ${registration_script} -d 3 -f ${atlas_anat} -m ${reference} -o ${sid}__output -t s -n ${params.register_processes}
    mv ${reference} ${sid}__native_anat.nii.gz

    """
}

// [sid, tractogram.trk, affine.mat, inverseWarp.nii.gz, atlas.nii.gz, t1.nii.gz]
sid_apply_registration
    .join(tractogram)
    .join(transformation_for_tractogram)
    .set{tractogram_registration} 

process Register_Streamlines {
    memory '10 GB'

    input:
    set sid, file(tractogram), file(affine), file(inverse_warp), file(atlas_anat), file(native_anat) from tractogram_registration

    output:
    set sid, "${sid}_output.trk", "${atlas_anat}" into tractogram_registered // [sid, output.trk, atlas.nii.gz]
    set sid, "${sid}_native.trk", "${native_anat}" into tractogram_native // [sid, native.trk, t1.nii.gz]

    script:
    """
    files="${tractogram}"
    if [[ \$( wc -w <<< \$files ) -gt 1 ]]
    then 
        echo \$files
        scil_streamlines_math.py concatenate \$files ${sid}_native.trk -f -vv
        scil_apply_transform_to_tractogram.py ${sid}_native.trk ${atlas_anat} \
        ${affine} ${sid}_output.trk \
        --inverse --in_deformation ${inverse_warp} -f --keep_invalid
    else
        echo \$files
        scil_apply_transform_to_tractogram.py ${tractogram} ${atlas_anat} \
        ${affine} ${sid}_output.trk \
        --inverse --in_deformation ${inverse_warp} -f --keep_invalid
        mv ${tractogram} ${sid}_native.trk
    fi

    """
}

// [sid, output.trk, atlas.nii.gz, native.trk, t1.nii.gz, model.pt, thresholds.json, atlas_dir/, config.json]
tractogram_registered
    .join(tractogram_native)
    .combine(model_for_binta)
    .combine(atlas_thresholds_for_binta)
    .combine(atlas_directory_for_binta)
    .combine(atlas_config_for_binta)
    .combine(device_for_binta)
    .set{filtering_channels}

process Filter_Streamlines {
    //cache false
    memory '30 GB'

    input:
    set sid, 
        file(tractogram), 
        file(atlas_anat), 
        file(native_tractogram), 
        file(native_anat), 
        file(model), 
        file(thresholds), 
        file(atlas_directory), 
        file(atlas_config), 
        val(device) from filtering_channels

    output:
    set sid, "*.trk" into bundles // [sid, AC_0.trk, AC_1.trk, ..., AF_L_0.trk, AF_L_1.trk, ...]

    script:
    """
    filter_streamline.py ${tractogram} ${atlas_directory} \
        ${model} ${atlas_anat} \
        ${thresholds} ${atlas_config} . \
        --original_tractogram ${native_tractogram} --original_reference ${native_anat}  -d ${device} -b 500000 -f -vv
    """
}

// [sid, AC_0.trk, AC_1.trk, ..., AF_L_0.trk, AF_L_1.trk, ..., config.json]
bundles.combine(atlas_config_for_concatenation).set{file_for_concatenation}

process Clean_Bundles {
    memory '5 GB'

    input:
    set sid, file(bundles), file(atlas_config) from file_for_concatenation

    output:
    // [sid, AC.trk, AF_L.trk, ...]
    set sid, "*.trk" into bundles_concatenated, bundle_for_gesta, bundles_finta

    script:
    """
    mkdir -p tmp 
    mv ${bundles} tmp
    cat "${atlas_config}" | jq -r '. | keys[]' |
    while IFS= read -r value; do
        echo Concatenating tmp/*\${value}* | echo "Done"
        scil_streamlines_math.py concatenate tmp/*\${value}* \${value}.trk -vv | echo "Done"
        mv \${value}.trk ${sid}__\${value}.trk | echo "Done"
    done
    """
}

bundle_for_gesta
    .join(transformation_for_gesta)
    .join(warp_for_gesta)
    .combine(model_for_gesta)
    .combine(atlas_directory_for_gesta)
    .combine(atlas_thresholds_for_gesta)
    .combine(atlas_config_for_gesta)
    .combine(device_for_gesta)
    .join(wm)
    .join(peaks)
    .set{files_for_gesta}

process GESTA {
    memory '10 GB'

    input:
    set sid, 
        file(bundles), 
        file(affine), file(inverse_warp), file(atlas_anat), file(native_anat),
        file(warp),
        file(model),
        file(atlas),
        file(thresholds),
        file(config),
        val(device),
        file(wm),
        file(peaks) from files_for_gesta

    output:
    // [sid, AC.trk, AF_L.trk, ...]
    set sid, "${sid}_*.trk" into bundles_augmented, bundles_gesta
    file "${sid}_wm_mni.nii.gz" 

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    mkdir -p mni 

    for b in ${bundles_list}
    do
        scil_apply_transform_to_tractogram.py \${b} ${atlas_anat} \
        ${affine} mni_\${b} \
        --inverse --in_deformation ${inverse_warp} -f --keep_invalid

        mv mni_\${b} mni/
    done

    antsApplyTransforms -d 3 -e 0 -i ${wm} -r ${atlas_anat} -o ${sid}_wm_mni.nii.gz -n NearestNeighbor -t ${warp} -t ${affine} -v 1
 
    generate_streamline.py \
	--in_bundles_MNI mni/*.trk \
	--in_bundles_native ${bundles} \
	--model ${model} \
	--reference_MNI ${atlas_anat} \
	--reference_native ${native_anat} \
	--thresholds_file ${thresholds} \
	--anatomy_file ${config} \
	--output . \
	--atlas_path ${atlas} \
	-d ${device} \
	-f -vv -n 20000 --ratio 0.5 0.5 -m 10000 \
	--wm_parc ${sid}_wm_mni.nii.gz \
	--peaks ${peaks} \
	--in_transfo ${affine} \
	--in_deformation ${warp} \
    --use_rs \
    --minL 20 \
    --maxL 220

    for f in *fodf_mask_20_220.trk;
    do
        scil_apply_transform_to_tractogram.py \$f ${native_anat} ${affine} ${sid}_\$f --in_deformation ${warp} --reverse_operation -f -vv
    done

    """
}

bundles_finta
    .join(bundles_gesta)
    .combine(atlas_config_for_finta_gesta).set{bundles_binta}


process BINTA {
    memory '5 GB'

    input:
    set sid, file(bundles_finta), file(bundles_gesta), file(config) from bundles_binta

    output:
    set sid, "${sid}__*concat.trk" into bundles_binta_gesta

    script:
    """
    cat ${config} | jq -r '. | keys[]' >keys.txt

    while IFS= read -r value; do
        count=`ls -1 *\${value}* 2>/dev/null | wc -l`
        if [[ \$count != 0 ]]; then
            echo Concatenating *\${value}*
            scil_streamlines_math.py concatenate *\${value}* ${sid}__\${value}_concat.trk -vv --no_metadata &>warnings.txt
        fi
    done <keys.txt
    """
}

// [sid, t1.nii.gz, AC.trk, AF_L.trk, ...]
reference_for_display_binta.join(bundles_concatenated).set{bundles_for_display_binta}

process Visualize_Bundles_BINTA {
    //errorStrategy 'retry'
    //maxRetries 3
    //memory '20 GB'

    input:
    set sid, file(anat), file(bundles) from bundles_for_display_binta

    output:
    file "*png"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_visualize_bundles_mosaic.py ${anat} ${bundles} ${sid}__bundles.png
    """
}

reference_for_display_gesta.join(bundles_augmented).set{bundles_for_display_gesta}

process Visualize_Bundles_GESTA {
    //errorStrategy 'retry'
    //maxRetries 3
    //memory '20 GB'

    input:
    set sid, file(anat), file(bundles) from bundles_for_display_gesta

    output:
    file "*png"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_visualize_bundles_mosaic.py ${anat} ${bundles} ${sid}__bundles.png
    """
}

reference_for_display_binta_gesta.join(bundles_binta_gesta).set{bundles_for_display_binta_gesta}

process Visualize_Bundles_BINTA_GESTA {
    //errorStrategy 'retry'
    //maxRetries 3
    //memory '20 GB'

    input:
    set sid, file(anat), file(bundles) from bundles_for_display_binta_gesta

    output:
    file "*png"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_visualize_bundles_mosaic.py ${anat} ${bundles} ${sid}__bundles.png
    """
}