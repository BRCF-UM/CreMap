# CreMAP

Shiny application for **mouse and human** spleen single-cell RNA-seq: **side-by-side** gene expression (UMAP + mean by cell type), **MouseMine (MGI)** Cre driver lookup mapped onto mouse and human where symbols match, cell-type labeling, and **Wilcoxon** differential expression (mouse object preferred when both are loaded).

## Requirements

- R 4.3+ recommended  
- [Seurat](https://satijalab.org/seurat/) 5.x (see `DESCRIPTION` for full imports)

The repo uses **[renv](https://rstudio.github.io/renv/)** so packages install into `renv/library` and do not touch your global R library (once `.Rprofile` is loaded).

## Test environment (easiest)

From the project root:

```bash
bash run_test.sh
```

That runs `scripts/install_deps.R` (first time can take several minutes while Seurat installs; **Bioconductor** packages for reference spleen defaults may also install and download hub data on first app launch), then starts the app at **http://127.0.0.1:8787**.

Or use Make:

```bash
make install   # once (or when dependencies change)
make run       # each time you want the app
```

Environment variables (optional):

| Variable | Default | Meaning |
|----------|---------|---------|
| `SHINY_HOST` | `127.0.0.1` | Bind address |
| `SHINY_PORT` | `8787` | Port |
| `SHINY_LAUNCH_BROWSER` | `true` | Open a browser automatically |
| `CREMAP_REFERENCE_DIR` | (empty) | Absolute path to a folder containing `mouse_spleen_reference.rds` and `human_spleen_reference.rds` (default: `<project>/reference`) |
| `CREMAP_ENABLE_CELLXGENE_HUMAN` | `0` | Set to `1` to allow CELLxGENE human fallback if He atlas fails |
| `CREMAP_ENABLE_REMOTE_REFERENCES` | local: `1`; shinyapps.io: `0` | Allow live atlas downloads when local RDS files are absent |
| `CREMAP_MAX_UPLOAD_MB` | `100` | Maximum Seurat RDS or Cell Ranger ZIP upload size |
| `CREMAP_MAX_UNCOMPRESSED_MB` | `500` | Maximum expanded size accepted from a Cell Ranger ZIP |

### Docker (optional)

```bash
docker compose up --build
```

Then open **http://localhost:8787**. The first image build compiles/installs Seurat and Bioconductor reference packages. For reliable defaults without hub traffic at runtime, run **`Rscript scripts/download_reference_data.R`** inside the container (or mount a host folder with the two RDS files and set **`CREMAP_REFERENCE_DIR`**). Otherwise the first launch may download Tabula Muris Senis and He atlas hub files (large, cached under the container user’s ExperimentHub cache).

### If you see “Missing packages: shiny, Seurat” after `run_test.sh`

That means packages are installed globally but **not** under `renv/library`. Re-run **`Rscript scripts/install_deps.R`** (the script was updated to always install into the project library). You should see a long install the first time, then `run_app.R` should work.

## Run manually in R

```r
setwd("/path/to/CreMAP")  # project root; loads renv via .Rprofile
shiny::runApp(".", launch.browser = TRUE)
```

Or open `app.R` in RStudio and click **Run App** (accept loading the project `.Rprofile` so renv activates).

## Data

- **Built-in default** (selected on launch):
  - **Recommended:** run **`Rscript scripts/download_reference_data.R`** once from the project root (network + Bioconductor; can take a long time). That writes **`reference/mouse_spleen_reference.rds`** and **`reference/human_spleen_reference.rds`**. The app **reads those files first**—no hub calls at runtime, which avoids most `DelayedArray` / `as.Seurat` / download issues.
  - If those two files are **missing**, the app builds **synthetic** data immediately, then loads Tabula Muris Senis (mouse) and He organ atlas (human) **in the background** (`later`). Optional CELLxGENE human fallback: set **`CREMAP_ENABLE_CELLXGENE_HUMAN=1`** before starting R.
- **Seurat RDS**: choose **Mouse** or **Human** and upload a `.rds` with `RNA` + `umap`. You can load one or both species over separate uploads (switch species, upload again).
- **Cell Ranger 10x matrix**: choose **species** (mouse/human), then folder path or ZIP as before. Load twice to fill both slots if desired.

Install reference dependencies after a pull (from R, with Bioconductor configured as usual):

```r
install.packages("BiocManager")
BiocManager::install(c("TabulaMurisSenisData", "scRNAseq", "cellxgenedp", "zellkonverter"))
```

Or from the shell: `Rscript scripts/install_deps.R` (installs into the project `renv` library when `.Rprofile` is active).

**ExperimentHub:** Only needed when **building** the local RDS files with `scripts/download_reference_data.R`, or when the app falls back to live hub loading. After the two reference RDS files exist under **`reference/`**, the running app does not need hub access for defaults. CreMAP sets **`EXPERIMENT_HUB_ASK`** / **`ANNOTATION_HUB_ASK`** to **`FALSE`** so cache prompts go to the default location. **`CREMAP_REFERENCE_DIR`** can point to another folder if you store the RDS files elsewhere.

## MGI / MouseMine

Cre alleles are queried live from [MouseMine](https://www.mousemine.org/mousemine) (Mouse Genome Informatics). A network connection is required for that tab. Driver gene symbols come from the `Allele.drivenBy.feature` association when present in MGI.

## Public deployment

CreMAP can be hosted on [shinyapps.io](https://www.shinyapps.io/) and embedded
in Google Sites. See **[DEPLOYMENT.md](DEPLOYMENT.md)** for the public GitHub
checklist, credentials setup, reference-data choices, bundle inspection, and
the included deployment command.

For hosted deployments, live atlas downloads are disabled by default and the
synthetic demonstration data is used unless both `reference/*.rds` files are
included in the deployment bundle. This avoids large, repeated ExperimentHub
downloads on ephemeral application workers.

## Notes

- Differential expression uses `Seurat::FindMarkers` with `ident.1` = reference cell type and `ident.2` = union of selected comparison types.
- Default reference loaders and paths live in `R/real_spleen_data.R`; local bundle file names are fixed (`mouse_spleen_reference.rds`, `human_spleen_reference.rds`). Re-run `scripts/download_reference_data.R` after changing Bioconductor versions if you need to refresh.
