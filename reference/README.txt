CreMAP reference data (optional local bundle)
==============================================

Place two Seurat objects here (or run the download script to create them):

  mouse_spleen_reference.rds
  human_spleen_reference.rds

Each must be a Seurat object with an RNA assay and a reduction named "umap".

To build these files from public atlases once (network + Bioconductor required):

  Rscript scripts/download_reference_data.R

The RDS files are gitignored by *.rds; keep them on disk for offline / reliable app startup.

Override directory with environment variable CREMAP_REFERENCE_DIR (absolute path).
