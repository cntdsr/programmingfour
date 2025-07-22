# Analysis of a 10X Genomics Single-Cell RNA-seq Dataset

This repository contains an R package, `programmingfour`, designed to perform a complete analysis pipeline on a single-cell RNA-seq dataset. The project fulfills requirements including data loading, gene annotation, filtering, dimensionality reduction (PCA, UMAP), clustering, and automated cell type annotation.

The project is structured to be fully reproducible using `renv` for dependency management.

---

## Project Structure

- **/R:** Contains all the custom functions developed for the analysis pipeline.
- **/man:** Documentation for each function.
- **/vignettes:** A detailed tutorial (`analysis_walkthrough.Rmd`) demonstrating how to use the package.
- **run_full_analysis.R:** The main script to execute the entire analysis from start to finish.
- **DESCRIPTION:** Defines the package and its dependencies.
- **renv.lock:** The `renv` lockfile, which ensures that the exact same package versions can be re-installed to guarantee reproducibility.

---

## 1. Data Acquisition

The raw data files are not included in this repository due to GitHub's file size limits. To run the analysis, you must first download the required data files and place them in the root directory of this project.

- **Single-cell Dataset (matrix, features, barcodes):** Download from (https://biotec.i-learn.unito.it/pluginfile.php/184389/mod_resource/content/0/exam_july.pdf)
- **Gene Information (GTF File):** Download from https://biotec.i-learn.unito.it/pluginfile.php/184389/mod_resource/content/0/exam_july.pdf

After downloading, ensure the following files are present in the project's root folder:
- `matrix.mtx.gz`
- `features.tsv.gz`
- `barcodes.tsv.gz`
- `Homo_sapiens.GRCh38.111.gtf`

---

## 2. Installation and Dependencies

This project uses `renv` to manage dependencies. To set up the environment and install the exact package versions used for this analysis, follow these steps in R or RStudio:

1.  **Clone or download this repository.**
2.  **Open the `programming4.Rproj` file in RStudio.**
3.  **Restore the environment using `renv`.** Run the following command in the R console. This will read the `renv.lock` file and install all required packages automatically.

    ```R
    renv::restore()
    ```
    This process might take a while.

---

## 3. How to Run the Analysis

Once the data is in place and the `renv` environment has been restored, you can run the entire analysis pipeline with a single command.

1.  Make sure you are in the RStudio project.
2.  Source the main script from the R console:

    ```R
    source("run_full_analysis.R")
    ```

The script will execute all analysis steps sequentially, printing status messages to the console and generating plots in the "Plots" pane of RStudio. The final, annotated Seurat object (`scrna_obj`) will be available in the R environment upon completion.
