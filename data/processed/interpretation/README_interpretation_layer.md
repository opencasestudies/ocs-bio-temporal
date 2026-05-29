# Case Study 5 K=10 Interpretation Layer

This directory contains additive interpretation tables for the revised K=10 selected model.
The tables summarize CoGAPS pattern usage at the subject level before aggregating across time or cell type.

## Assumptions

- The R K=10 selected-model output is used as the primary table source.
- R/Python agreement is tracked separately from the comparison directory.
- Subjects, not cells, are treated as the replicate unit for the interpretation summaries.
- CoGAPS weights remain nonnegative; expression direction is supplied only when the directionality output is available.

## Inputs

- Model output: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/cogaps_r_revised_model_K10_seed2_iter2000`
- Run stem: `cogaps_K10_seed2_iter2000`
- R/Python comparison: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/cogaps_revised_model_comparison_K10_seed2_iter2000`
- Revised K-selection: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10`
- Directionality: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/pattern_directionality_revised_model_K10_seed2_iter2000`

## Dataset Snapshot

- Discovery cells: `4800`
- Subjects: `3`
- Timepoints: `8`
- Broad cell types: `6`
- IFN candidate pattern: `Pattern3`

## Key Output Tables

- `subject_level/subject_pattern_means.csv`
- `subject_level/subject_time_pattern_means.csv`
- `subject_level/subject_celltype_pattern_means.csv`
- `subject_level/subject_time_celltype_pattern_means.csv`
- `pattern_summaries/pattern_by_time_subject_summary.csv`
- `pattern_summaries/pattern_by_cell_type_subject_summary.csv`
- `pattern_summaries/pattern_by_time_cell_type_subject_summary.csv`
- `rq_tables/RQ1_temporal_programs.csv`
- `rq_tables/RQ2_identity_vs_activity.csv`
- `rq_tables/RQ3_ifn_associated_program.csv`
- `pattern_annotation_table.csv`
