# FIESTA Flow

[![DOI tractolearn](https://zenodo.org/badge/DOI/10.5281/zenodo.7562790.svg)](https://doi.org/10.5281/zenodo.7562790)
[![DOI RBX](https://zenodo.org/badge/DOI/10.5281/zenodo.7562635.svg)](https://doi.org/10.5281/zenodo.7562635)

FIESTA Flow is a nextflow pipeline for white matter bundles recognition in tractography. 

To use this pipeline you must have access to 6 documents - namely: 
* best_model_contrastive_tractoinferno_hcp.pt
* rbx_atlas_v10.json
* mni_masked.nii.gz
* rbx_atlas (folder containing WM atlas bundles, AC.trk, AF_L.tr, ...)
* thresholds_contrastive_tractoinferno_hcp.json
* number_rejection_sampling.json
* max_total_sampling.json
* ratio.json
* degree.json
* white_matter_mask.json
* docker image felixdumais1/fiesta:`<version>`


## Requirements

* Nextflow
* Singularity/Docker

## Launch pipeline

To launch the pipeline, replace `<version>` by the current fiesta version, and run:

```
nextflow run main.nf --input inputs --model /path/to/best_model_contrastive_tractoinferno_hcp.pt --atlas_config /path/to/rbx_atlas_v10.json --atlas_anat /path/to/mni_masked.nii.gz --atlas_directory /path/to/data/atlas/pop_average --atlas_thresholds /path/to/thresholds_contrastive_tractoinferno_hcp.json --device cuda --number_rejection_sampling_config /path/to/number_rejection_sampling.json --max_total_sampling_config /path/to/max_total_sampling.json --batch_sampling 5000 --ratio_atlas_bundle_config /path/to/ratio.json --register_processes 24 --processes 24 --degree_config /path/to/degree.json --output_dir results --white_matter_config /path/to/white_matter_mask.json --fa_threshold 0.1 -resume -with-singularity felixdumais1/fiesta:<version>
```

## USAGE

See *USAGE* or run `nextflow run main.nf --help`

## How to cite

Please refer to the `tractolearn` ["How to cite"](https://github.com/scil-vital/tractolearn#how-to-cite)
section to know how to cite this work.

## Patent

J. H. Legarreta, M. Descoteaux, and P.-M. Jodoin. “PROCESSING OF TRACTOGRAPHY
RESULTS USING AN AUTOENCODER”. Filed 03 2021. Imeka Solutions Inc. United States
Patent #17/337,413. Pending.

## License

This software is distributed under a particular license. Please see the
[*LICENSE*](LICENSE) file for details.
