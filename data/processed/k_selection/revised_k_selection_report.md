# Revised K Selection for Case Study 5

Generated: `2026-05-27T04:10:56+00:00`

## Decision

The revised selected model for the main Case Study 5 analysis is:

- `K = 10`
- `seed = 2`
- `nIterations = 2000`

The previous sweep selected `K = 5` because the original rule preferred the smallest K on the seed-to-seed stability plateau. That rule is useful for a compact baseline, but it does not fully match the proposal's primary biological goal: resolving temporal dengue-response programs, especially a D10 interferon-associated activity program.

The revised rule is:

> Select the smallest stability-plateau K that resolves at least one IFN-candidate, activity-like temporal program.

This selects `K = 10`.

## Why Reassess K?

The proposal asks whether CoGAPS can recover temporal immune programs, separate cell identity from infection activity, and identify interferon-stimulated antiviral responses. A low-rank model can be highly stable while compressing several biological processes into broad lineage-associated factors. The K-selection criterion therefore needs to include biological resolution, not only numerical reproducibility.

This follows the identity/activity framing from Kotliar et al. (2019), where the useful rank depends on the resolution desired by the analyst, and the CoGAPS protocol emphasis on annotating patterns by biological process rather than selecting dimensionality mechanically. Fertig et al. (2014) also motivates using CoGAPS specifically to interpret time-course patterns as overlapping biological processes. Waickman et al. (2021) provides the biological prior that the strongest experimental DENV-1 host response occurs around day 10 and includes interferon/inflammatory genes.

## Stability-Plateau Candidates

| K | selected_n_iter | stability_core | overall_score | max_eta_timepoint | candidate_ifn_pattern_count_mean | activity_like_pattern_count_mean | ifn_secondary_score | redundancy_penalty |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5 | 8000 | 0.996 | 0.888 | 0.013 | 0.000 | 0.000 | 0.000 | 0.390 |
| 6 | 4000 | 0.993 | 0.916 | 0.026 | 0.000 | 0.000 | 0.100 | 0.304 |
| 8 | 6000 | 0.989 | 0.908 | 0.066 | 1.000 | 0.000 | 0.800 | 0.340 |
| 9 | 4000 | 0.986 | 0.910 | 0.067 | 2.000 | 0.000 | 0.800 | 0.326 |
| 10 | 2000 | 0.979 | 0.914 | 0.461 | 1.000 | 1.000 | 1.000 | 0.312 |
| 12 | 4000 | 0.979 | 0.895 | 0.094 | 2.000 | 0.000 | 0.600 | 0.337 |

## Revised Ranking

The goal-aligned score is a transparent summary used for review, not a black-box selector:

- 35% seed-to-seed stability
- 15% low redundancy
- 20% maximum temporal effect size
- 15% mean IFN-candidate pattern count
- 15% mean activity-like pattern count

| K | selected_n_iter | stability_core | low_redundancy_score | max_eta_timepoint | candidate_ifn_pattern_count_mean | activity_like_pattern_count_mean | goal_aligned_score | on_stability_plateau |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 40 | 8000 | 0.827 | 1.000 | 0.567 | 1.000 | 3.600 | 0.864 | False |
| 36 | 6000 | 0.816 | 0.951 | 0.508 | 1.000 | 3.200 | 0.816 | False |
| 20 | 8000 | 0.917 | 0.738 | 0.524 | 1.000 | 2.000 | 0.775 | False |
| 28 | 8000 | 0.846 | 0.873 | 0.511 | 1.000 | 2.200 | 0.774 | False |
| 24 | 4000 | 0.861 | 0.809 | 0.513 | 1.000 | 2.200 | 0.771 | False |
| 16 | 2000 | 0.963 | 0.629 | 0.526 | 1.000 | 1.600 | 0.759 | False |
| 32 | 8000 | 0.882 | 0.818 | 0.509 | 1.000 | 1.600 | 0.753 | False |
| 10 | 2000 | 0.979 | 0.621 | 0.461 | 1.000 | 1.000 | 0.715 | True |
| 11 | 2000 | 0.878 | 0.521 | 0.466 | 1.000 | 1.400 | 0.683 | False |
| 9 | 4000 | 0.986 | 0.511 | 0.067 | 2.000 | 0.000 | 0.595 | True |

## Representative K=10 Pattern Evidence

Representative run: `K=10`, `seed=2`, `nIterations=2000`.

| pattern | eta_timepoint | eta_broad_cell_type | pattern_class | spearman_ifn_score_rho | peak_timepoint | ifn_top_gene_overlap_top15 | candidate_ifn_pattern |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Pattern5 | 0.461 | 0.161 | activity-like | 0.757 | D10 | 6 | True |
| Pattern10 | 0.021 | 0.164 | identity-like | -0.007 | D14/15 | 0 | False |
| Pattern3 | 0.019 | 0.164 | identity-like | -0.075 | D6 | 0 | False |
| Pattern9 | 0.010 | 0.049 | identity-like | -0.091 | D14/15 | 0 | False |
| Pattern8 | 0.009 | 0.911 | identity-like | -0.171 | D0 | 0 | False |
| Pattern7 | 0.009 | 0.779 | identity-like | 0.206 | D14/15 | 0 | False |
| Pattern6 | 0.006 | 0.498 | identity-like | -0.147 | D0 | 0 | False |
| Pattern4 | 0.005 | 0.548 | identity-like | -0.086 | D14/15 | 0 | False |
| Pattern2 | 0.003 | 0.645 | identity-like | -0.136 | D14/15 | 0 | False |
| Pattern1 | 0.002 | 0.164 | identity-like | 0.075 | D28 | 0 | False |

The leading candidate pattern is `Pattern5`:

- temporal effect size: `0.461`
- broad-cell-type effect size: `0.161`
- class: `activity-like`
- IFN-score Spearman rho: `0.757`
- peak timepoint: `D10`
- top-15 IFN gene overlap: `6`
- top genes: `LY6E, IFITM1, ISG15, IFITM2, IFI6, IFITM3, MX1, PSME2, IRF7, XAF1, HLA-A, B2M, MT2A, MYL12A, IFI44L`

This aligns with the proposed case-study biology: an interferon-associated program peaking at D10, rather than a purely lineage-dominant low-rank factor.

## Interpretation

`K=5` remains useful as a compact stability baseline. It demonstrates that the strongest low-rank structure in this balanced PBMC subset is lineage/identity-associated. However, because the case study's central question is temporal immune activity, `K=10` is the better primary analysis model: it remains on the stability plateau and resolves a D10 IFN-like activity program.

The case study should present this as a methodological lesson: the best K depends on the analysis goal. Stability is necessary, but biological resolution determines whether the chosen model can answer the research questions.

## Local Supporting References

- `references/literature/fertig_2014_temporal_cogaps.pdf`
- `references/literature/johnson_2023_cogaps_protocol_nature_protocols.pdf`
- `references/literature/kotliar_2019_identity_activity_cnmf.pdf`
- `references/literature/stein_obrien_2019_sc_cogaps_transfer_learning.pdf`
- `references/literature/sharma_2020_projectr.pdf`
- `references/literature/waickman_2021_dengue_single_cell_pbmc.pdf`

## Output Files

- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/revised_k_selection_summary_by_k.csv`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/stability_plateau_candidates.csv`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/representative_K10_seed2_iter2000_patterns.csv`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/selected_pairwise_reproducibility.csv`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/revised_k_selection_manifest.json`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/figures/stability_vs_temporal_effect.png`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/figures/plateau_biology_signals.png`
- `/Users/othomas/Desktop/CASE STUDY 5 SWEEP/GSE154386/revised_k_selection_K10/figures/revised_k_goal_aligned_score.png`

