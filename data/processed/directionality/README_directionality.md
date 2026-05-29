# GSE154386 Selected CoGAPS Pattern Directionality

This post-sweep analysis asks whether the genes defining the selected CoGAPS patterns are up- or down-regulated relative to D0.

## Inputs

- Expression object: `/Users/othomas/Desktop/CS5_sweep_results/cache/gse154386_preprocessed_hvg.h5ad`
- Selected top genes: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/cogaps_r_revised_model_K10_seed2_iter2000/cogaps_K10_seed2_iter2000.top_genes.csv`
- Pattern summary: `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/cogaps_r_revised_model_K10_seed2_iter2000/cogaps_K10_seed2_iter2000.pattern_summary.csv`
- Run stem: `cogaps_K10_seed2_iter2000`
- Top genes per pattern tested: `25`
- Baseline timepoint: `D0`
- Direction threshold: `abs(log2FC) >= 0.25`

## Expression Matrix

- Cells analyzed before cohort filtering: `180859`
- Cells analyzed after `experimental` cohort filtering: `156182`
- HVGs available: `5000`

## Interpretation Notes

- CoGAPS gene weights are nonnegative, so they identify pattern-defining genes but do not encode up/down regulation by themselves.
- Directionality here is estimated from pseudobulk expression contrasts using raw counts from `layers['counts']`.
- The primary contrast is each post-baseline experimental timepoint versus the same subject's D0 baseline.
- Cell-type-stratified outputs should be used to avoid confusing cell-type composition with within-cell-type regulation.

## Strongest All-Cell Pattern-Timepoint Summaries

| pattern | timepoint | median gene log2FC | frac genes up | frac genes down | consensus |
|---|---:|---:|---:|---:|---|
| Pattern3 | D14/15 | 1.562 | 0.96 | 0.04 | mostly_up |
| Pattern3 | D10 | 1.457 | 0.88 | 0.00 | mostly_up |
| Pattern10 | D14/15 | -0.567 | 0.16 | 0.76 | mostly_down |
| Pattern7 | D10 | -0.474 | 0.12 | 0.72 | mostly_down |
| Pattern9 | D14/15 | 0.462 | 0.56 | 0.36 | mixed |
| Pattern10 | D10 | -0.439 | 0.08 | 0.68 | mostly_down |
| Pattern7 | D14/15 | -0.428 | 0.24 | 0.64 | mostly_near_zero |
| Pattern1 | D10 | -0.397 | 0.04 | 0.68 | mostly_down |
| Pattern6 | D10 | -0.397 | 0.12 | 0.64 | mostly_down |
| Pattern1 | D14/15 | -0.339 | 0.16 | 0.56 | mostly_near_zero |

## Output Files

- `gene_directionality_by_subject.csv`: subject-level log2FC contrasts for each selected gene.
- `gene_directionality_summary.csv`: gene-level summaries across subjects.
- `pattern_directionality_summary.csv`: pattern-level summaries of top-gene direction.
- `figures/all_patterns_top_gene_log2fc_all_cells.svg`: combined heatmap of every selected top-gene row across all patterns.
- `figures/`: additional per-pattern and cell-type-stratified SVG heatmaps.
