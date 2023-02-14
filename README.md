# FIESTA Flow

[![documentation](https://readthedocs.org/projects/tractolearn/badge/?version=latest)](https://tractolearn.readthedocs.io/en/latest/?badge=latest)
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
* docker image felixdumais1/tractolearn-docker:`<version>`


## Requirements

* Nextflow
* Singularity/Docker

## Launch pipeline

To launch the pipeline, replace `<version>` by the current fiesta version, and run:

```
nextflow run main.nf --input inputs --model 
/path/to/best_model_contrastive_tractoinferno_hcp.pt --atlas_config 
/path/to/rbx_atlas_v10.json --atlas_anat /path/to/mni_masked.nii.gz 
--atlas_directory /path/to/data/atlas/pop_average --atlas_thresholds 
/path/to/thresholds_contrastive_tractoinferno_hcp.json --device cuda 
--number_rejection_sampling_config /path/to/number_rejection_sampling.json 
--max_total_sampling_config /path/to/max_total_sampling.json 
--batch_sampling 5000 --ratio_atlas_bundle_config /path/to/ratio.json 
--register_processes 24 --processes 24 --degree_config /path/to/degree.json 
--output_dir results --white_matter_config /path/to/white_matter_mask.json 
--fa_threshold 0.1 -resume -with-singularity 
felixdumais1/tractolearn-docker:<version>
```

## USAGE

See *USAGE* or run `nextflow run main.nf --help`

## How to cite

If you use this toolkit in a scientific publication or if you want to cite
our previous works, we would appreciate if you considered the following aspects:
- If you use `tractolearn`, please add a link to the appropriate code, data or
  related resource hosting service (e.g., repository, PyPI) from where you
  obtained `tractolearn`. You may want to include the specific version or commit
  hash information for the sake of reproducibility.
- Please, cite the appropriate scientific works:
  - If you use `tractolearn` to filter implausible streamlines or you want to
    cite our work in tractography filtering, cite [FINTA] and [FIESTA].
  - If you want to cite our work in tractography bundling, cite [CINTA] and
    [FIESTA].
    - If you use `tractolearn` to bundle streamlines using a k-nearest neighbor
      label approach, cite [CINTA].
    - If you use `tractolearn` to bundle streamlines using a thresholding
      approach, cite [FINTA] and [FIESTA].
  - If you use `tractolearn` for generative purposes or you want to cite our
    work in generative models for tractography, cite [GESTA] and [FIESTA].
  - If you use parts of `tractolearn` for other purposes, please generally cite
    [FINTA] and [FIESTA].

The corresponding `BibTeX` files are contained in the above links.

If you use the [data](https://zenodo.org/record/7562790) made available by the
authors, please cite the appropriate Zenodo record.

Please reach out to us if you have related questions.

## Patent

J. H. Legarreta, M. Descoteaux, and P.-M. Jodoin. “PROCESSING OF TRACTOGRAPHY
RESULTS USING AN AUTOENCODER”. Filed 03 2021. Imeka Solutions Inc. United States
Patent #17/337,413. Pending.

## License

This software is distributed under a particular license. Please see the
[*LICENSE*](LICENSE) file for details.


[FINTA]: ./doc/bibtex/Legarreta21_-_MIA_-_FINTA.bib "Filtering in tractography using autoencoders (FINTA)"
[CINTA]: ./doc/bibtex/Legarreta22_-_MICCAI-CDMRI_-_CINTA.bib "Clustering in Tractography Using Autoencoders (CINTA)"
[GESTA]: ./doc/bibtex/Legarreta22_-_arXiv_-_GESTA.bib "Generative sampling in tractography using autoencoders (GESTA)"
[FIESTA]: ./doc/bibtex/Dumais22_-_arXiv_-_FIESTA.bib "FIESTA: Autoencoders for accurate fiber segmentation in tractography"

