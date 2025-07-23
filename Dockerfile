# STEP 1: Base image 
FROM rocker/r-ver:4.3.1

# STEP 2: System dependencies
# This installs a comprehensive list of system libraries needed for R package compilation.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libzstd-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    build-essential \
    pandoc \
    pandoc-citeproc \
    libglpk-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# STEP 3: Set working directory 
WORKDIR /analysis

# STEP 4: Copy our project files 
# We copy only the necessary files to keep the container clean.
COPY R ./R
COPY DESCRIPTION .
COPY NAMESPACE .
COPY run_full_analysis.R .
COPY barcodes.tsv.gz .
COPY features.tsv.gz .
COPY matrix.mtx.gz .
COPY Homo_sapiens.GRCh38.111.gtf .

# STEP 5: Install specific R dependencies (CON DIPENDENZE BIOCONDUCTOR AGGIUNTE)
# Aggiungiamo 'SingleCellExperiment' e 'SummarizedExperiment' che sono dipendenze
# fondamentali per SingleR e celldex.
RUN R -e "options(Ncpus = parallel::detectCores()); \
          install.packages('BiocManager'); \
          BiocManager::install(c('Seurat', 'SingleR', 'celldex', 'rtracklayer', 'ggplot2', 'dplyr', 'data.table', 'remotes', 'devtools', 'SingleCellExperiment', 'SummarizedExperiment'), update=FALSE, ask=FALSE)"


# STEP 6: Install our local package 
# We already copied the source code of our package (R/, DESCRIPTION, NAMESPACE).
# Now we install it directly from the local files within the container.
# This removes the dependency on GitHub during the build.
RUN R -e "devtools::install('/analysis')"

# STEP 8: Final execution command (RENAMED from 7 to 8)
# This will run your main analysis script when the container starts.
CMD ["Rscript", "run_full_analysis.R"]