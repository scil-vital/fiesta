#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "FINTA Multibundles Flow"
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
log.info "Atlas Config: $params.atlas_config"
log.info "Atlas Anat: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info "Atlas Thresholds: $params.atlas_thresholds"
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
    .set{ tractogram } // [sid, tractogram.trk]

Channel
    .fromPath("$root/**/*.nii.gz")
    .map{[it.parent.name, it]}
    .set{ reference } // [sid, t1.nii.gz]

if (!(params.atlas_anat) || !(params.atlas_config) || !(params.atlas_directory)) {
    error "You must specify all 3 atlas related input. --atlas_anat, " +
    "--atlas_config and --atlas_directory all are mandatory."
}

atlas_anat = Channel.fromPath("$params.atlas_anat")

Channel.fromPath("$params.atlas_config").into{atlas_config; atlas_config_for_concatenation}
atlas_directory = Channel.fromPath("$params.atlas_directory")
model = Channel.fromPath("$params.model")
atlas_thresholds = Channel.fromPath("$params.atlas_thresholds")

reference
    .combine(atlas_anat)
    .set{reference_atlas_anat} 

process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(reference), file(atlas_anat) from reference_atlas_anat

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, fixed_t1.nii.gz]
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1InverseWarp.nii.gz", "${atlas_anat}" into transformation_for_tractogram
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1Warp.nii.gz", "${sid}__outputInverseWarped.nii.gz" into inverse_transformation_for_tractogram
    file "${sid}__outputWarped.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyN.sh -d 3 -f ${atlas_anat} -m ${reference} -o ${sid}__output -t s -n ${params.register_processes}
    """
}

// [sid, tractogram.trk, affine.mat, inverseWarp.nii.gz, outputWarped.nii.gz]
tractogram.join(transformation_for_tractogram).set{tractogram_registration} 

process Register_Streamlines {
    memory '20 GB'

    input:
    set sid, file(tractogram), file(affine), file(inverse_warp), file(atlas_anat) from tractogram_registration

    output:
    set sid, "${sid}_output.trk", "${atlas_anat}" into tractogram_registered // [sid, output.trk]

    script:
    """
    files="${tractogram}"
    if [[ \$( wc -w <<< \$files ) -gt 1 ]]
    then 
        echo \$files
        scil_streamlines_math.py concatenate \$files out.trk -f -vv
        scil_apply_transform_to_tractogram.py out.trk ${atlas_anat} \
        ${affine} ${sid}_output.trk \
        --inverse --in_deformation ${inverse_warp} -f -vv
    else
        echo \$files
        scil_apply_transform_to_tractogram.py ${tractogram} ${atlas_anat} \
        ${affine} ${sid}_output.trk \
        --inverse --in_deformation ${inverse_warp} -f -vv
    fi

    """
}

tractogram_registered
    .combine(model)
    .combine(atlas_thresholds)
    .combine(atlas_directory)
    .combine(atlas_config)
    .set{filtering_channels}

process Filter_Streamlines {
    memory '20 GB'

    input:
    set sid, file(tractogram), file(atlas_anat), file(model), file(thresholds), file(atlas_directory), file(atlas_config) from filtering_channels

    output:
    set sid, "*.trk" into bundles // [sid, output.trk]

    script:
    """
    filter_streamline.py ${tractogram} ${atlas_directory} \
        ${model} ${atlas_anat} \
        ${thresholds} ${atlas_config} \
        . -d cuda -b 500000 -f -vv
    """
}

bundles.combine(atlas_config_for_concatenation).set{file_for_concatenation}

process Concatenating_Bundles {
    memory '5 GB'

    input:
    set sid, file(bundles), file(atlas_config) from file_for_concatenation

    output:
    set sid, "*.trk" into bundles_concatenated

    script:
    """
    mkdir -p tmp 
    mv ${bundles} tmp
    cat "${atlas_config}" | jq -r '. | keys[]' |
    while IFS= read -r value; do
        echo Concatenating tmp/*\${value}*
        scil_streamlines_math.py concatenate tmp/*\${value}* \${value}.trk -vv | echo "Done"
    done
    """
}

bundles_concatenated.join(inverse_transformation_for_tractogram).set{files_for_inverse_transforms}

process Registering_in_Native {
    memory '5 GB'
    input:
    set sid, file(bundles), file(affine), file(warp), file(reference) from files_for_inverse_transforms

    output:
    file "*.trk"

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list
    do
        filename=\$( basename \$bundle )
        scil_apply_transform_to_tractogram.py \$bundle ${reference} ${affine} \$filename --in_deformation ${warp} --reverse_operation -f -vv
        mv \$filename ${sid}__\$filename
    done
    
    """
}