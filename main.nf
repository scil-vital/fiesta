#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count",
                "device":"$params.device",
                "registration_speed":"$params.registration_speed",
                "ratio_atlas_bundle":"$params.ratio_atlas_bundle_config",
                "number_rejection_sampling":"$params.number_rejection_sampling_config",
                "parzen_window_seeds": "$params.parzen_window_seeds",
                "max_total_sampling": "$params.max_total_sampling_config",
                "batch_sampling": "$params.batch_sampling",
                "degree": "$params.degree_config",
                "fa_threshold": "$params.fa_threshold",
                "run_gesta": "$params.run_gesta",
                "white_matter": "$params.white_matter_config"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "BINTA Flow"
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
log.info "Model file: $params.model"
log.info ""
log.info "[Atlas]"
log.info "Atlas Config: $params.atlas_config"
log.info "Atlas Anat: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info "Atlas Thresholds: $params.atlas_thresholds"
log.info ""
log.info "[Processing]"
log.info "Device: $params.device"
log.info "Registration speed: $params.registration_speed"
log.info ""
log.info "[Gesta]"
log.info "Ratio Atlas Bundle Config: $params.ratio_atlas_bundle_config"
log.info "Number Rejection Sampling Config: $params.number_rejection_sampling_config"
log.info "Number of Parzen Window seeds: $params.parzen_window_seeds"
log.info "Max total sampling config: $params.max_total_sampling_config"
log.info "Batch sampling: $params.batch_sampling"
log.info "Acceptance angle config: $params.degree_config"
log.info "WM config $params.white_matter_config"
log.info "FA Threshold : $params.fa_threshold"
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
    .into{ reference; reference_for_display_cinta; reference_for_display_gesta; reference_for_display_binta; reference_for_check } // [sid, t1.nii.gz]

Channel
    .fromPath("$root/**/*wm.nii.gz")
    .map{[it.parent.name, it]}
    .into{ wm; wm_for_check; wm_for_filtering } // [sid, wm.nii.gz]

Channel
    .fromPath("$root/**/*fa.nii.gz")
    .map{[it.parent.name, it]}
    .into{ fa; fa_for_check; fa_for_filtering } // [sid, fa.nii.gz]

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
    .into{atlas_config_for_cinta; atlas_config_for_concatenation; atlas_config_for_gesta; atlas_config_for_finta_gesta}

Channel.fromPath("$params.atlas_directory")
    .into{atlas_directory_for_cinta; atlas_directory_for_gesta}

Channel.fromPath("$params.model")
    .into{model_for_cinta; model_for_gesta}

Channel.fromPath("$params.atlas_thresholds")
    .into{atlas_thresholds_for_cinta; atlas_thresholds_for_gesta}

Channel.fromPath("$params.number_rejection_sampling_config")
    .set{number_rejection_sampling}

Channel.fromPath("$params.max_total_sampling_config")
    .set{max_total_sampling}

Channel.fromPath("$params.ratio_atlas_bundle_config")
    .set{ratio_atlas_bundle}

Channel.fromPath("$params.degree_config")
    .set{degree}

Channel.fromPath("$params.white_matter_config")
    .set{white_matter}

Channel.value("$params.device")
    .into{device_for_cinta; device_for_gesta}

tractogram_for_check
    .join(reference_for_check)
    .join(wm_for_check)
    .join(peaks_for_check)
    .join(fa_for_check)
    .set{compatibility_check}


process Check_Files_Compatibility {
    errorStrategy 'ignore'

    input:
    set sid, 
        file(tractogram), 
        file(reference),
        file(wm),
        file(peaks),
        file(fa) from compatibility_check // [sid, tractogram.trk, t1.nii.gz, wm.nii.gz, peaks.nii.gz, fa.nii.gz]

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, atlas.nii.gz, t1.nii.gz]
    val sid into sid_registration, sid_apply_registration

    script:
    """
    compatibility=\$(scil_verify_space_attributes_compatibility.py ${tractogram} ${reference} ${wm} ${peaks} ${fa})
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
    memory '24 GB'

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
    .combine(model_for_cinta)
    .combine(atlas_thresholds_for_cinta)
    .combine(atlas_directory_for_cinta)
    .combine(atlas_config_for_cinta)
    .combine(device_for_cinta)
    .set{files_for_cinta}

process FINTA_Multibundle {
    cache false
    memory '15 GB'

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
        val(device) from files_for_cinta

    output:
    set sid, "*.trk" into bundles // [sid, AC_0.trk, AC_1.trk, ..., AF_L_0.trk, AF_L_1.trk, ...]

    script:
    """
    ae_bundle_streamlines ${tractogram} ${atlas_directory} \
        ${model} ${atlas_anat} \
        ${thresholds} ${atlas_config} . \
        --original_tractogram ${native_tractogram} --original_reference ${native_anat} -d ${device} -b 500000 -f -vv

    mkdir tmp
    mv *trk tmp
    mv -t . tmp/${native_tractogram} tmp/${tractogram}

    if [ -z "\$(ls -A tmp)" ]; then
        touch empty.trk
    else
        rm tmp/*implausible* | echo "Done"
        mv tmp/* .
    fi
    """
}

// [sid, AC_0.trk, AC_1.trk, ..., AF_L_0.trk, AF_L_1.trk, ..., config.json]
bundles.combine(atlas_config_for_concatenation).set{file_for_concatenation}

process Concatenating_FINTA_Multibundle {
    memory '2 GB'

    input:
    set sid, file(bundles), file(atlas_config) from file_for_concatenation

    output:
    // [sid, AC.trk, AF_L.trk, ...]
    set sid, "*.trk" into bundles_concatenated, bundle_for_gesta, bundles_cinta

    script:
    """
    mkdir cat_bundles
    mkdir -p tmp 
    mv ${bundles} tmp
    cwd=\$(pwd)
    cat "${atlas_config}" | jq -r '. | keys[]' |
    while IFS= read -r value; do
        if [[ -n "\$(ls -A tmp/*\${value}* 2>/dev/null)" ]]; then
            echo Concatenating tmp/*\${value}*
            cd tmp
            for v in *\${value}*; do
                scil_count_streamlines.py \${v}
                scil_remove_invalid_streamlines.py \${v} no_invalid_\${v}
                scil_count_streamlines.py no_invalid_\${v}
            done
            cd \${cwd}
            scil_streamlines_math.py concatenate tmp/no_invalid*\${value}* \${value}.trk -vv
            mv \${value}.trk cat_bundles/${sid}__\${value}.trk
        fi
    done

    if [ -z "\$(ls -A cat_bundles)" ]; then
        echo "Empty !"
        touch empty_0.trk
    else
        echo "Not Empty !"
        mv cat_bundles/* .
    fi

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
    .combine(number_rejection_sampling)
    .combine(max_total_sampling)
    .combine(degree)
    .combine(ratio_atlas_bundle)
    .combine(white_matter)
    .join(wm)
    .join(peaks)
    .join(fa)
    .set{files_for_gesta}

process GESTA_Fast {
    memory '15 GB'

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
        file(number_rejection_sampling),
        file(max_total_sampling),
        file(degree),
        file(ratio_atlas_bundle),
        file(white_matter),
        file(wm),
        file(peaks),
        file(fa) from files_for_gesta

    when:
    params.run_gesta

    output:
    // [sid, AC.trk, AF_L.trk, ...]
    set sid, "${sid}_*.trk" into bundles_augmented, bundles_gesta
    file "${sid}_wm_mni.nii.gz" 
    file "${sid}_fa_mni.nii.gz" 

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    echo $number_rejection_sampling
    mkdir -p mni 

    for b in ${bundles_list}
    do
        if [[ \$b != *"empty"* ]]; then
            echo \${b}
            scil_apply_transform_to_tractogram.py \${b} ${atlas_anat} \
            ${affine} mni_\${b} \
            --inverse --in_deformation ${inverse_warp} --keep_invalid -f 
            scil_remove_invalid_streamlines.py mni_\${b} mni_\${b} -f
            mv mni_\${b} mni/
        fi
    done

    antsApplyTransforms -d 3 -e 0 -i ${wm} -r ${atlas_anat} -o ${sid}_wm_mni.nii.gz -n NearestNeighbor -t ${warp} -t ${affine} -v 1
    antsApplyTransforms -d 3 -e 0 -i ${fa} -r ${atlas_anat} -o ${sid}_fa_mni.nii.gz -n Linear -t ${warp} -t ${affine} -v 1
    
    ae_generate_streamlines \
	--in_bundles_common_space mni/*.trk \
	--model ${model} \
	--reference_common_space ${atlas_anat} \
	--reference_native ${native_anat} \
	--anatomy_file ${config} \
	--output . \
	--atlas_path ${atlas} \
	-d ${device} \
	-f -vv \
    -n ${number_rejection_sampling} \
    --max_total_sampling ${max_total_sampling} \
    --ratio ${ratio_atlas_bundle} \
    -m $params.parzen_window_seeds \
	--wm_parc_common_space ${sid}_wm_mni.nii.gz \
    --fa_common_space ${sid}_fa_mni.nii.gz \
	--peaks ${peaks} \
	--in_transfo ${affine} \
	--in_deformation ${warp} \
    --use_rs \
    --batch_sampling $params.batch_sampling \
    --minL 20 \
    --maxL 220 \
    --degree ${degree} \
    -a \
    --threshold_fa $params.fa_threshold \
    --white_matter_config ${white_matter}

    for f in *fodf_mask_20_220*.trk;
    do
        scil_count_streamlines.py \$f > count.json
        count=\$(jq ".[].streamline_count" count.json)

        if [[ \$count -gt 0 ]]; then
            scil_apply_transform_to_tractogram.py \$f ${native_anat} ${affine} to_concatenate_\$f --in_deformation ${warp} --reverse_operation --keep_invalid -f
            scil_remove_invalid_streamlines.py to_concatenate_\$f to_concatenate_\$f -f
        fi
    done

    if find . -name "to_concatenate*fodf_mask_20_220*trk" | grep -q .; then
        mkdir -p tmp
        mv to_concatenate*fodf_mask_20_220*trk tmp
        cat "${config}" | jq -r '. | keys[]' |
        while IFS= read -r value; do
            echo Concatenating tmp/*\${value}* | echo "Done"
            scil_streamlines_math.py concatenate tmp/*\${value}* \${value}.trk -vv | echo "Done"
            mv \${value}.trk ${sid}__\${value}_gesta.trk | echo "Done"
        done
    else
        touch ${sid}__empty_gesta.trk
    fi

    """
}

bundles_cinta
    .join(bundles_gesta)
    .combine(atlas_config_for_finta_gesta).set{bundles_for_binta}


process FIESTA {
    memory '5 GB'

    input:
    set sid, file(bundles_finta), file(bundles_gesta), file(config) from bundles_for_binta

    output:
    set sid, "*binta*.trk" into bundles_binta

    script:
    finta_bundle_list = bundles_finta.join(", ").replace(',', '')
    """
    cat ${config} | jq -r '. | keys[]' >keys.txt

    if find . -name "*__empty_gesta.trk" | grep -q .; then
        echo "No generated file"
        for b in ${finta_bundle_list};
        do
            mv \${b} binta_\${b}
        done
    else
        while IFS= read -r value; do
            count=`ls -1 *\${value}* 2>/dev/null | wc -l`
            if [[ \$count != 0 ]]; then
                echo Concatenating *\${value}*
                scil_streamlines_math.py concatenate *\${value}* ${sid}__\${value}_binta.trk -vv --no_metadata --ignore_invalid &>warnings.txt
            fi
        done <keys.txt
    fi
    """
}

bundles_binta.join(wm_for_filtering).set{ files_for_filtering }

process FIESTA_Clean_Bundles {
    memory '5 GB'

    input:
    set sid, file(bundles_binta), file(wm) from files_for_filtering

    output:
    set sid, "*binta_cleaned*.trk" into bundles_binta_cleaned
    file "*eroded_${wm}"

    script:
    String bundles_list = bundles_binta.join(", ").replace(',', '')
    """
    if find . -name "*empty*.trk" | grep -q .; then
        echo "No file"
        for b in ${bundles_list};
        do
            mv \${b} binta_cleaned_\${b}
        done
        touch empty_eroded_${wm}
    else
        scil_image_math.py convert ${wm} ${wm} --data_type int16 -f
        for bundle in ${bundles_list};
        do
            filename=\$(basename -- "\$bundle")
            filename="\${filename%.*}"

            scil_count_streamlines.py \${bundle} >streamline_count.json
            streamline_count=\$(jq --arg keyvar "\${filename}" '.[\$keyvar].streamline_count' streamline_count.json)

            if ((streamline_count>1)); then
                echo \${streamline_count}
                scil_remove_invalid_streamlines.py \${bundle} \${filename}_no_invalid.trk
                scil_count_streamlines.py \${filename}_no_invalid.trk

                scil_detect_streamlines_loops.py \${filename}_no_invalid.trk \${filename}_no_loops.trk -a 360 -f
                scil_image_math.py erosion ${wm} 3 eroded_${wm} -f
                scil_filter_tractogram.py \${filename}_no_loops.trk \${filename}_cleaned.trk --drawn_roi eroded_${wm} either_end exclude -f
            fi
        done
    fi
    """
}

// [sid, t1.nii.gz, AC.trk, AF_L.trk, ...]
reference_for_display_cinta.join(bundles_concatenated).set{bundles_for_display_cinta}

process Visualize_Bundles_FINTA_Multibundle {
    // errorStrategy 'ignore'
    //maxRetries 3
    //memory '20 GB'

    input:
    set sid, file(anat), file(bundles) from bundles_for_display_cinta

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

process Visualize_Bundles_GESTA_Fast {
    // errorStrategy 'ignore'
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

reference_for_display_binta.join(bundles_binta_cleaned).set{bundles_for_display_binta}

process Visualize_Bundles_FIESTA_Clean_Bundles {
    // errorStrategy 'ignore'
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