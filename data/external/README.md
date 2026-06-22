# External Data Notes

Large files are not staged for GitHub commit. The reproduction archive will be
published through Zenodo, and the final DOI will be added after the record is
public. The case-study pages use small included evidence files from
`data/processed/`; the larger files below are for learners, instructors, or
reviewers who want to audit or rebuild the full workflow.

For the learner-facing full reproduction path, mount the downloaded Zenodo
reproduction archive inside Docker at:

`/home/rstudio/project/data/external/reproduction_archive`

| Artifact | Status | Expected local source | Notes |
|---|---|---|---|
| `GSE154386_RAW.tar` | GEO download or Zenodo archive | `GSE154386/GSE154386_RAW.tar` or `reproduction_archive/source/GSE154386_RAW.tar` | Full GEO raw download. Not committed to GitHub. |
| `gse154386_preprocessed_hvg.h5ad` | forthcoming Zenodo archive | `reproduction_archive/cache/gse154386_preprocessed_hvg.h5ad` | Large optional preprocessing cache. Not required for GitHub render. |
| `gse154386_experimental_discovery_cells_x_genes.h5ad` | forthcoming Zenodo archive | `reproduction_archive/cache/gse154386_experimental_discovery_cells_x_genes.h5ad` | Primary R selected-model input for optional local rerun. |
| `gse154386_experimental_discovery_genes_x_cells.h5ad` | forthcoming Zenodo archive | `reproduction_archive/cache/gse154386_experimental_discovery_genes_x_cells.h5ad` | Primary Python selected-model input for optional local rerun. |
| `gse154386_natural_projection_target_cells_x_genes.h5ad` | forthcoming Zenodo archive | `reproduction_archive/cache/gse154386_natural_projection_target_cells_x_genes.h5ad` | Optional natural projection target; not required for the core case-study lesson. |
