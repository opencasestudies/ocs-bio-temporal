# Case Study 5 Reproduction Scripts

These scripts document the analysis workflow used to create the small analysis
tables and figures in the case study. They are included for provenance and
optional full reproduction. The main case-study analysis uses the precomputed
files in `data/processed/` and does not run these scripts.

The full reproduction workflow requires large files distributed through the
case-study Zenodo record, including preprocessed H5AD objects, CoGAPS input
objects, model outputs, and sweep summaries. Keep those large files outside the
GitHub repository.

## Main Groups

- Sweep setup and aggregation:
  - `gse154386_cogaps_sweep_config.sh`
  - `gse154386_make_cogaps_jobs_tsv.py`
  - `gse154386_cogaps_sweep_array_singleprocess_rhino_light.sbatch`
  - `gse154386_cogaps_prep_cache.py`
  - `gse154386_run_one_singleprocess.py`
  - `gse154386_aggregate_results.py`
  - `cs5_generate_revised_k_selection_report.py`
- Selected K = 10 model runners:
  - `cs5_run_selected_model_r.R`
  - `cs5_run_selected_model_python.py`
  - `cs5_cogaps_python_utils.py`
- Interpretation-layer tables:
  - `cs5_build_k10_interpretation_layer.py`
- Directionality analysis:
  - `gse154386_pattern_directionality.py`
- Shared source workflow used by several reproduction scripts:
  - `gse154386_sparse_distributed_cogaps.py`

The script names preserve the development history. In particular,
`gse154386_sparse_distributed_cogaps.py` contains the local source workflow
used by the non-distributed Case Study 5 analysis even though its filename
contains the word `distributed`.
