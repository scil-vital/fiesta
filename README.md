# FINTA Flow
===================

FINTA Flow is a nextflow pipeline using the FINTA approach for white matter bundles recognition in tractography. 

To use this pipeline you must have access to 6 documents - namely: 
* best_model_ae.pt
* finta_multibundles.sif
* mni_masked.nii.gz
* rbx_atlas (folder containing WM atlas bundles, AC.trk, AF_L.tr, ...)
* rbx.json
* thresholds_ae.json

If you use this pipeline, please cite:

```
Legarreta, J. H. et al. Filtering in tractography using autoencoders (FINTA). Medical Image Analysis 72, 102126 (2021)
```


Requirements
------------

* Nextflow

Singularity/Docker
-----------

To launch the pipeline you can run:

```
nextflow run main.nf --input inputs/ --model /path/to/best_model_ae.pt --atlas_config /path/to/rbx.json --atlas_anat /path/to/mni_masked.nii.gz --atlas_directory /path/to/rbx_atlas/ --atlas_thresholds /path/to/thresholds_ae.json -with-singularity /path/to/finta_multibundles.sif -resume
```


Usage
-----

See *USAGE* or run `nextflow run main.nf --help`

