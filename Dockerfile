# Optional containerized test environment (first build installs Seurat; can take a while).
FROM rocker/shiny:4.4.1

WORKDIR /srv/shiny-server/CreMAP
COPY . .

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libhdf5-dev libxml2-dev libcurl4-openssl-dev libssl-dev libgsl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript --vanilla -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages(c('shiny', 'bslib', 'dplyr', 'plotly', 'httr2', 'Matrix', 'DT', 'Seurat', 'SeuratObject'), dependencies = TRUE) \
"

RUN Rscript --vanilla -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
  BiocManager::install( \
    c('SingleCellExperiment', 'SummarizedExperiment', 'ExperimentHub', 'TabulaMurisSenisData', 'scRNAseq', 'cellxgenedp', 'zellkonverter'), \
    ask=FALSE, update=FALSE \
  ) \
"

EXPOSE 8787
ENV SHINY_HOST=0.0.0.0
ENV SHINY_PORT=8787
ENV SHINY_LAUNCH_BROWSER=false

CMD ["R", "--vanilla", "-e", "shiny::runApp('/srv/shiny-server/CreMAP', host='0.0.0.0', port=8787L, launch.browser=FALSE)"]
