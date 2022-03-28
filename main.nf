#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["bundles_config":"$params.bundles_config",
                "to_register":"$params.to_register",
                "register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "FINTA Multibundles Quantitative Evaluation Flow"
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
log.info "[Evaluation options]"
log.info "Register moving to fixed: $params.to_register"
log.info ""
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

log.info "Input: $params.input"
root = file(params.input)

if (!(params.bundles_config)) {
    error "You must specify --bundles_config."
}

Channel
    .fromPath("$params.bundles_config")
    .into{bundles_config_for_volumetric_metrics;
          bundles_config_for_density_metrics;
          bundles_config_for_streamline_metrics;
          bundles_config_for_weighted_density_metrics}

Channel
    .fromFilePairs("$root/**/Fixed/*.trk",
                   size: -1) { it.parent.parent.name }
    .set{fixed_bundles} // [sid, [AF_L.trk, AF_R.trk, ...]]

Channel
    .fromFilePairs("$root/**/Moving/*.trk",
                   size: -1) { it.parent.parent.name }
    .into{moving_bundles_for_registration;moving_bundles_native} // [sid, [AF_L.trk, AF_R.trk, ...]]

Channel
    .fromPath("$root/**/Fixed/*.nii.gz")
    .map{[it.parent.parent.name, it]}
    .into{fixed_t1_native; fixed_t1_for_registration; fixed_t1_for_volume} // [sid, t1.nii.gz]
Channel
    .fromPath("$root/**/Moving/*.nii.gz")
    .map{[it.parent.parent.name, it]}
    .into{moving_t1_native; moving_t1_for_registration; moving_t1_for_volume} // [sid, t1.nii.gz]

fixed_t1_for_registration.join(moving_t1_for_registration).set{fixed_moving_t1} // [sid, fixed_t1.nii.gz, moving_t1.nii.gz]

process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, "fixed_t1.nii.gz", "moving_t1.nii.gz" from fixed_moving_t1

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, fixed_t1.nii.gz]
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1InverseWarp.nii.gz", "${sid}__fixed_t1.nii.gz" into transformation_for_tractogram
    // [sid, registeredImage.nii.gz]
    set sid, "${sid}__outputWarped.nii.gz" into moving_t1_registered

    when: 
    params.to_register

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f fixed_t1.nii.gz -m moving_t1.nii.gz -n ${params.register_processes} -o ${sid}__output -t s
    cp fixed_t1.nii.gz ${sid}__fixed_t1.nii.gz
    """
}

// [sid, [AF_L.trk, AF_R.trk, ...], affine.mat, inverseWarp.nii.gz, fixed_t1.nii.gz]
moving_bundles_for_registration.join(transformation_for_tractogram).set{moving_necessities} 

process Register_Streamlines {
    input:
    set sid, file(bundles), file(affine), file(inverse_warp), file(reference) from moving_necessities

    output:
    // [sid, [AF_L.trk, AF_R.trk, ...] in fixed space]
    set sid, "*.trk" into moving_bundles_registered

    when: 
    params.to_register

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
    do
        scil_apply_transform_to_tractogram.py \$bundle ${reference} \
        ${affine} "${sid}__\$bundle" \
        --inverse --in_deformation ${inverse_warp} -f -vv
    done
    """
}

if (params.to_register) {
    // [sid, [AF_L.trk, AF_R.trk, ...], t1.nii.gz]
    fixed_bundles
        .join(fixed_t1_native)
        .into{fixed_bundles_for_streamline_metrics;fixed_bundles_for_volumetric_metrics}

    // [sid, [AF_L.trk, AF_R.trk, ...], t1.nii.gz]
    moving_bundles_registered
        .join(moving_t1_registered)
        .into{moving_bundles_for_streamline_metrics;moving_bundles_for_volumetric_metrics}
} else {
    // [sid, [AF_L.trk, AF_R.trk, ...], t1.nii.gz]
    fixed_bundles
        .join(fixed_t1_for_volume)
        .into{fixed_bundles_for_streamline_metrics;fixed_bundles_for_volumetric_metrics}

    // [sid, [AF_L.trk, AF_R.trk, ...], t1.nii.gz]
    moving_bundles_native
        .join(moving_t1_for_volume)
        .into{moving_bundles_for_streamline_metrics;moving_bundles_for_volumetric_metrics}
}

process Compute_Bundle_Volumes_Fixed {
    input:
    set sid, file(fixed_bundles), "fixed_t1.nii.gz" from fixed_bundles_for_volumetric_metrics

    output:
    set sid, "*binary.nii.gz" into fixed_volumes // [sid, [AF_L_fixed_binary.nii.gz, AF_R_fixed_binary.nii.gz, ...]]
    set sid, "*density.nii.gz" into fixed_volumes_density // [sid, [AF_L_fixed_density.nii.gz, AF_R_fixed_density.nii.gz, ...]]


    script:
    String fixed_bundles_list = fixed_bundles.join(", ").replace(',', '')
    """
    for bundle in $fixed_bundles_list;
    do
        bname=\$(basename \$bundle .trk)
        scil_compute_streamlines_density_map.py \$bundle ${sid}__\${bname}_fixed_density.nii.gz --reference fixed_t1.nii.gz
        scil_image_math.py upper_threshold_eq ${sid}__\${bname}_fixed_density.nii.gz 0 ${sid}__\${bname}_fixed_binary_invert.nii.gz
        scil_image_math.py invert ${sid}__\${bname}_fixed_binary_invert.nii.gz ${sid}__\${bname}_fixed_binary.nii.gz
    done
    """
}

process Compute_Bundle_Volumes_Moving {
    input:
    set sid, file(moving_bundles), "moving_t1.nii.gz" from moving_bundles_for_volumetric_metrics

    output:
    set sid, "*binary.nii.gz" into moving_volumes // [sid, [AF_L_moving_binary.nii.gz, AF_R_moving_binary.nii.gz, ...]]
    set sid, "*density.nii.gz" into moving_volumes_density // [sid, [AF_L_moving_density.nii.gz, AF_R_moving_density.nii.gz, ...]]

    script:
    String moving_bundles_list = moving_bundles.join(", ").replace(',', '')
    """
    mrinfo moving_t1.nii.gz
    for bundle in $moving_bundles_list;
    do
        bname=\$(basename \$bundle .trk)

        if [[ "\$bundle" == *"${sid}"* ]]; then
            name_density=\${bname}_moving_density.nii.gz
            name_binary_invert=\${bname}_moving_binary_invert.nii.gz
            name_binary=\${bname}_moving_binary.nii.gz
        else
            name_density=${sid}__\${bname}_moving_density.nii.gz
            name_binary_invert=${sid}__\${bname}_moving_binary_invert.nii.gz
            name_binary=${sid}__\${bname}_moving_binary.nii.gz
        fi

        scil_compute_streamlines_density_map.py \$bundle \$name_density --reference moving_t1.nii.gz
        scil_image_math.py upper_threshold_eq \$name_density 0 \$name_binary_invert
        scil_image_math.py invert \$name_binary_invert \$name_binary
    done
    """
}

// [sid, [AF_L_fixed_binary.nii.gz, AF_R_fixed_binary.nii.gz, ...], [AF_L_moving_binary.nii.gz, AF_R_moving_binary.nii.gz, ...], bundles_config.json]
fixed_volumes
    .join(moving_volumes)
    .combine(bundles_config_for_volumetric_metrics)
    .set{file_for_volumes_metrics} 

process Volumetric_Metrics {
    input:
    set sid, file(fixed_bundles), file(moving_bundles), file(config) from file_for_volumes_metrics

    output:
    file "${sid}__volumetric_metrics.png"
    file "${sid}__volumetric_metrics.json"
    file "${sid}__fixed_bundles.json"
    file "${sid}__moving_bundles.json"
    file "*cm*"

    script:
    """
    metrics.py --bundles_path_A ${fixed_bundles} \
               --bundles_path_B ${moving_bundles} \
               --bundles_dict ${config} -t volumetric -vv
    cp *metrics.png ${sid}__volumetric_metrics.png
    cp *metrics.json ${sid}__volumetric_metrics.json
    cp *bundles_a.json ${sid}__fixed_bundles.json
    cp *bundles_b.json ${sid}__moving_bundles.json
    """
}

fixed_volumes_density
    .join(moving_volumes_density)
    .combine(bundles_config_for_density_metrics)
    .set{file_for_density_metrics} 

process Density_Metrics {
    input:
    set sid, file(fixed_bundles), file(moving_bundles), file(config) from file_for_density_metrics

    output:
    file "${sid}__density_metrics.png"
    file "${sid}__density_metrics.json"
    file "${sid}__fixed_bundles.json"
    file "${sid}__moving_bundles.json"
    file "*cm*"

    script:
    """
    metrics.py --bundles_path_A ${fixed_bundles} \
               --bundles_path_B ${moving_bundles} \
               --bundles_dict ${config} -t density -vv
    cp *metrics.png ${sid}__density_metrics.png
    cp *metrics.json ${sid}__density_metrics.json
    cp *bundles_a.json ${sid}__fixed_bundles.json
    cp *bundles_b.json ${sid}__moving_bundles.json
    """
}

process Prepare_Bundle_Streamline_Fixed {
    input:
    set sid, file(fixed_bundles), "fixed_t1.nii.gz" from fixed_bundles_for_streamline_metrics

    output:
    // [sid, [AF_L_fixed.trk, AF_R_fixed.trk, ...]]
    set sid, "fixed_t1.nii.gz", "*trk" into fixed_streamlines, fixed_bundles_for_display, fixed_bundles_weighted_density

    script:
    String fixed_bundles_list = fixed_bundles.join(", ").replace(',', '')
    """
    for bundle in $fixed_bundles_list;
    do
        bname=\$(basename \$bundle .trk)
        cp \$bundle ${sid}__\${bname}_fixed.trk
    done
    """
}

process Prepare_Bundle_Streamline_Moving {
    input:
    set sid, file(moving_bundles), "moving_t1.nii.gz" from moving_bundles_for_streamline_metrics

    output:
    // [sid, [AF_L_moving.trk, AF_R_moving.trk, ...]]
    set sid, "moving_t1.nii.gz", "*trk" into moving_streamlines, moving_bundles_for_display, moving_bundles_weighted_density

    script:
    String moving_bundles_list = moving_bundles.join(", ").replace(',', '')
    """
    for bundle in $moving_bundles_list;
    do
        bname=\$(basename \$bundle .trk)

        if [[ "\$bundle" == *"${sid}"* ]]; then
            name=\${bname}_moving.trk
        else
            name=${sid}__\${bname}_moving.trk
        fi

        cp \$bundle \$name
    done
    """
}

fixed_bundles_weighted_density
    .join(moving_bundles_weighted_density)
    .combine(bundles_config_for_weighted_density_metrics)
    .set{file_for_weighted_density_metrics} 

process Weighted_Density_Metrics {
    input:
    set sid, "fixed_t1.nii.gz", file(fixed_bundles), "moving_t1.nii.gz", file(moving_bundles), file(config) from file_for_weighted_density_metrics

    output:
    file "${sid}__weighted_density_metrics.png"
    file "${sid}__weighted_density_metrics.json"
    file "${sid}__fixed_bundles.json"
    file "${sid}__moving_bundles.json"
    file "*cm*"

    script:
    """
    metrics.py --bundles_path_A ${fixed_bundles} \
               --bundles_path_B ${moving_bundles} \
               --bundles_dict ${config} -t weighted_density -vv
    cp *metrics.png ${sid}__weighted_density_metrics.png
    cp *metrics.json ${sid}__weighted_density_metrics.json
    cp *bundles_a.json ${sid}__fixed_bundles.json
    cp *bundles_b.json ${sid}__moving_bundles.json
    """
}

fixed_bundles_for_display.join(moving_bundles_for_display).set{bundles_for_display}

process Visualize_Bundles {
    input:
    set sid, "fixed_t1.nii.gz", file(fixed_bundles), "moving_t1.nii.gz", file(moving_bundles) from bundles_for_display

    output:
    file "*png"


    script:
    String fixed_bundles_list = fixed_bundles.join(", ").replace(',', '')
    """
    scil_visualize_bundles_mosaic.py fixed_t1.nii.gz ${fixed_bundles} ${sid}__fixed_bundles.png
    scil_visualize_bundles_mosaic.py moving_t1.nii.gz ${moving_bundles} ${sid}__moving_bundles.png
    """
}

fixed_streamlines
    .join(moving_streamlines)
    .combine(bundles_config_for_streamline_metrics)
    .set{file_for_streamline_metrics} 

process Streamline_Metrics {
    input:
    set sid, "fixed_t1.nii.gz", file(fixed_bundles), "moving_t1.nii.gz", file(moving_bundles), file(config) from file_for_streamline_metrics

    output:
    file "${sid}__streamline_metrics.png"
    file "${sid}__streamline_metrics.json"
    file "${sid}__fixed_bundles.json"
    file "${sid}__moving_bundles.json"
    file "*cm*"

    script:
    """
    metrics.py --bundles_path_A ${fixed_bundles} \
               --bundles_path_B ${moving_bundles} \
               --bundles_dict ${config} -t streamline -vv
    cp *metrics.png ${sid}__streamline_metrics.png
    cp *metrics.json ${sid}__streamline_metrics.json
    cp *bundles_a.json ${sid}__fixed_bundles.json
    cp *bundles_b.json ${sid}__moving_bundles.json
    """
}