# --- PASSO 0: INSTALLAZIONE E CARICAMENTO LIBRERIE ---

# Installazione di Seurat
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

# Installazione di tidyverse
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# Installazione di BiocManager e rtracklayer
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

# Caricamento delle librerie
library(Seurat)
library(tidyverse)
library(rtracklayer)

print("Librerie installate e caricate con successo.")

# --- TASK 1: GENE ANNOTATION ---

# Definisci il percorso della cartella contenente i dati
# Assicurati che il percorso sia corretto rispetto alla posizione dello script
data_dir <- "." 

# Carica i dati 10X
# La funzione Read10X cerca i file matrix.mtx.gz, barcodes.tsv.gz, e features.tsv.gz
expression_matrix <- Read10X(data.dir = data_dir)

# Crea l'oggetto Seurat
# L'oggetto Seurat è il contenitore centrale per i dati e le analisi
scrna_obj <- CreateSeuratObject(counts = expression_matrix, project = "scRNA_Project")

# Visualizza un riassunto dell'oggetto
print(scrna_obj)

# Percorso del file GTF
gtf_file <- file.path(data_dir, "Homo_sapiens.GRCh38.111.gtf")

# Importa il file GTF
gtf_data <- rtracklayer::import(gtf_file)

# Converti l'oggetto GRanges in un data.frame per una più facile ispezione
gtf_df <- as.data.frame(gtf_data)

# Filtra per mantenere solo i geni di tipo 'protein_coding'
# La colonna 'gene_biotype' contiene questa informazione
protein_coding_genes <- gtf_df %>%
  filter(type == "gene" & gene_biotype == "protein_coding") %>%
  pull(gene_name) %>% # Estrai i nomi dei geni
  unique() # Assicurati che ogni nome di gene sia unico

# Stampa il numero di geni protein-coding trovati
cat(paste("Numero di geni 'protein-coding' unici trovati nel GTF:", length(protein_coding_genes), "\n"))

# Verifica il numero di geni prima del filtro
n_genes_before <- nrow(scrna_obj)
cat(paste("Numero di geni nell'oggetto Seurat prima del filtro:", n_genes_before, "\n"))

# we retain only genes that are present in our count matrix
protein_coding_in_matrix <- intersect(rownames(scrna_obj), protein_coding_genes)
scrna_obj <- subset(scrna_obj, features = protein_coding_in_matrix)

# Verifica il numero di geni dopo il filtro
n_genes_after <- nrow(scrna_obj)
cat(paste("Numero di geni nell'oggetto Seurat dopo il filtro:", n_genes_after, "\n"))
cat(paste("Numero di geni rimossi:", n_genes_before - n_genes_after, "\n"))


# --- TASK 2: GENE EXPRESSION SUMMARY ---

# For each cell, calculate the number of genes with expression >= 3 UMIs.
# This metric will be stored in the Seurat object's metadata.

# 1. Get the raw count matrix from the Seurat object
counts_matrix <- GetAssayData(object = scrna_obj, layer = "counts")

# 2. For each column (cell), count how many genes have a value of 3 or more.
# The result is a vector where each element corresponds to a cell.
n_genes_ge3 <- colSums(counts_matrix >= 3)

# 3. Add this new metric as a new column in the Seurat object's metadata.
# This makes it easy to plot and access later.
scrna_obj$nGenes_ge3 <- n_genes_ge3

# 4. Display the distribution of these counts using a violin plot.
# We use Seurat's VlnPlot function, which is built on top of ggplot2.
VlnPlot(scrna_obj, features = "nGenes_ge3", pt.size = 0.1) +
  NoLegend() +
  labs(
    title = "Distribution of Genes with High Expression per Cell",
    y = "Number of Genes (UMI count >= 3)",
    x = "All Cells"
  )


# --- TASK 3: GENE FILTERING ---

# We need to exclude three categories of genes:
# 1. Mitochondrial genes
# 2. Ribosomal protein genes
# 3. Ribosomal pseudogenes

# First, let's get the list of all genes currently in our object for filtering.
current_genes <- rownames(scrna_obj)

# 1. Identify Mitochondrial genes. They are typically prefixed with "MT-".
mito_genes <- grep(pattern = "^MT-", x = current_genes, value = TRUE, ignore.case = TRUE)

# 2. Identify Ribosomal protein genes. They are prefixed with "RPS" or "RPL".
ribo_genes <- grep(pattern = "^RP[SL]", x = current_genes, value = TRUE, ignore.case = TRUE)

# 3. Identify Ribosomal pseudogenes. We will use our GTF annotation for this.
# We look for genes with a 'pseudogene' biotype whose names also suggest a ribosomal origin.
ribo_pseudo_genes_info <- gtf_df %>%
  filter(
    # Find entries whose gene biotype contains the word "pseudogene"
    grepl("pseudogene", gene_biotype, ignore.case = TRUE) &
      # And whose gene name starts with RPS or RPL
      grepl("^RP[SL]", gene_name, ignore.case = TRUE)
  ) %>%
  pull(gene_name) %>% # Get the names
  unique() # Ensure they are unique

# Find which of these pseudogenes are actually in our current dataset
ribo_pseudo_genes <- intersect(current_genes, ribo_pseudo_genes_info)


# --- Create the Summary Table ---

# Create a data frame to summarize the number of genes removed from each category.
summary_table <- data.frame(
  Category = c("Mitochondrial Genes", "Ribosomal Protein Genes", "Ribosomal Pseudogenes"),
  GenesRemoved = c(length(mito_genes), length(ribo_genes), length(ribo_pseudo_genes))
)

# Print the summary table to the console.
print("Summary of Genes to be Removed:")
print(summary_table)


# --- Perform the Filtering ---

# Combine all lists of genes to be removed into a single vector.
genes_to_remove <- unique(c(mito_genes, ribo_genes, ribo_pseudo_genes))

# Get the list of genes to KEEP by finding the difference between all genes and those to remove.
genes_to_keep <- setdiff(current_genes, genes_to_remove)

# Print the numbers before and after for verification.
cat(paste("\nNumber of genes before filtering:", length(current_genes), "\n"))
cat(paste("Total unique genes to be removed:", length(genes_to_remove), "\n"))

# Subset the Seurat object to keep only the desired genes.
scrna_obj <- subset(scrna_obj, features = genes_to_keep)

cat(paste("Number of genes after filtering:", nrow(scrna_obj), "\n"))


# --- TASK 4: PRINCIPAL COMPONENT ANALYSIS (PCA) ---

# NOTE: Before running PCA, we must perform standard pre-processing steps:
# 1. Normalization: To account for differences in sequencing depth per cell.
# 2. Identification of Highly Variable Features: To focus the analysis on genes
#    that contribute most to cell-to-cell variation.

# 1. Normalization
# We use the LogNormalize method. It normalizes the counts for each cell by the total counts,
# multiplies by a scale factor (10,000 by default), and log-transforms the result.
scrna_obj <- NormalizeData(scrna_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 2. Finding Highly Variable Features (HVGs)
# These are genes that show high cell-to-cell variation in the dataset.
# Focusing on these genes in downstream analysis helps to highlight biological signal.
scrna_obj <- FindVariableFeatures(scrna_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scaling the Data
# This is a crucial step before PCA. It shifts the expression of each gene so that the
# mean expression across cells is 0 and the variance is 1. This gives equal weight to each
# gene in the PCA, so that highly-expressed genes do not dominate the analysis.
# We only scale the variable features for efficiency.
all_genes <- rownames(scrna_obj)
scrna_obj <- ScaleData(scrna_obj, features = all_genes)


# 4. Running PCA
# Now we can run the PCA on the scaled data.
# By default, Seurat uses the previously identified variable features.
scrna_obj <- RunPCA(scrna_obj, features = VariableFeatures(object = scrna_obj))


# 5. Assess the variance explained by the first 20 principal components
# Seurat provides a built-in visualization for this, often called an "Elbow Plot".
# We will create a histogram as specifically requested by the task.

# Get the standard deviations of the principal components from the Seurat object
pca_stdev <- Stdev(scrna_obj, reduction = "pca")

# Calculate the variance explained by each PC (variance = stdev^2)
variance_explained <- pca_stdev^2

# Calculate the percentage of variance explained
total_variance <- sum(variance_explained)
percent_variance <- (variance_explained / total_variance) * 100

# Create a data frame for plotting the first 20 PCs
variance_df <- data.frame(
  PC = 1:20,
  PercentVariance = percent_variance[1:20]
)

# Visualize using a histogram (bar plot in ggplot2)
# A bar plot is more appropriate here than a true histogram for showing discrete PC values.
ggplot(variance_df, aes(x = PC, y = PercentVariance)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = 1:20) + # Ensure every PC number is shown on the x-axis
  labs(
    title = "Variance Explained by the First 20 Principal Components",
    x = "Principal Component (PC)",
    y = "Percentage of Variance Explained (%)"
  ) +
  theme_bw()


# --- TASK 5: UMAP VISUALIZATION ---

# We will now run the UMAP algorithm on the results of our PCA.
# It is crucial to run UMAP on the reduced dimensions (PCs), not on the full gene expression matrix.

# 1. Run UMAP
# We specify the Seurat object and which dimensions (PCs) to use as input.
# 'dims = 1:13' tells UMAP to use the first 13 Principal Components.
# 'reduction = "pca"' specifies that these dimensions come from the PCA results.
scrna_obj <- RunUMAP(scrna_obj, reduction = "pca", dims = 1:13)

# 2. Generate a UMAP Visualization
# We use DimPlot to create a 2D scatter plot of the UMAP results.
# Each point is a cell, and its position is determined by the UMAP coordinates.
# At this stage, we have not performed clustering, so all points will be the same color.
# The purpose is to visually inspect the overall structure of the data.

DimPlot(scrna_obj, reduction = "umap") +
  labs(
    title = "UMAP Visualization of Cell Populations",
    subtitle = "Based on the first 13 Principal Components"
  ) +
  theme_bw()


# --- TASK 6: CLUSTERING ---

# Seurat's clustering workflow is graph-based. It involves two steps:
# 1. Build a graph representing cell-to-cell similarity.
# 2. Apply a community detection algorithm to partition the graph.

# 1. Find Neighbors
# This function constructs a Shared Nearest Neighbor (SNN) graph. It first identifies the
# nearest neighbors for each cell within the PCA space we've already defined.
# We use the same 13 PCs that we used for the UMAP.
scrna_obj <- FindNeighbors(scrna_obj, reduction = "pca", dims = 1:13)

# 2. Find Clusters
# This function applies a community detection algorithm (Louvain algorithm by default)
# to the graph we just built.
# The 'resolution' parameter controls the granularity of the clustering.
# Higher values lead to more clusters. 0.5 is a standard starting point.
scrna_obj <- FindClusters(scrna_obj, resolution = 0.5)

# After running FindClusters, the cluster IDs for each cell are stored in the
# Seurat object's metadata. We can see how many cells are in each cluster.
cat("\nNumber of cells in each cluster:\n")
print(table(scrna_obj$seurat_clusters))

# 3. Visualize the Clusters
# We use DimPlot again, which will now automatically color the cells by their
# newly assigned cluster ID.
# The 'label = TRUE' argument adds a numbered label to the center of each cluster.
DimPlot(scrna_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() +
  labs(
    title = "UMAP Visualization with Cell Clusters",
    subtitle = "Louvain algorithm with resolution 0.5"
  ) +
  theme_bw()

# --- TASK 7: CELL TYPE ANNOTATION (VERSIONE FINALE E COMPLETA) ---

# 1. Installazione di tutti i pacchetti necessari in modo non interattivo
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SingleR", quietly = TRUE)) BiocManager::install("SingleR", update = FALSE, ask = FALSE)
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex", update = FALSE, ask = FALSE)
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment", update = FALSE, ask = FALSE)


# 2. Caricamento di tutte le librerie
library(Seurat)
library(SingleR)
library(celldex)
library(ExperimentHub)
library(SingleCellExperiment) 


# 3. Impostazione della cache 
hub_cache_dir <- tools::R_user_dir("ExperimentHub", which = "cache")
dir.create(hub_cache_dir, recursive = TRUE, showWarnings = FALSE)
setExperimentHubOption("CACHE", hub_cache_dir)
cat("Cache di ExperimentHub assicurata in:", hub_cache_dir, "\n")


# 4. Caricamento dei dati di riferimento
cat("Caricamento dati di riferimento (HumanPrimaryCellAtlasData)...\n")
ref_data <- HumanPrimaryCellAtlasData()
cat("Dati di riferimento caricati.\n")


# 5. Preparazione dei dati per SingleR 
sce_obj <- as.SingleCellExperiment(scrna_obj)


# 6. Esecuzione di SingleR
cat("Esecuzione di SingleR per l'annotazione dei tipi cellulari...\n")
predictions <- SingleR(
  test = sce_obj,
  ref = ref_data,
  labels = ref_data$label.main
)
cat("Annotazione completata.\n")


# 7. Aggiunta delle etichette predette all'oggetto Seurat
scrna_obj$singleR_labels <- predictions$labels[match(rownames(scrna_obj@meta.data), rownames(predictions))]


# 8. Visualizzazione e salvataggio del grafico dei risultati
output_plot_path <- "umap_singleR_annotation.png"
png(output_plot_path, width = 8, height = 7, units = "in", res = 300)
print(
  DimPlot(scrna_obj, reduction = "umap", group.by = "singleR_labels", label = TRUE, repel = TRUE) +
    NoLegend() +
    labs(
      title = "UMAP con Annotazioni Automatiche dei Tipi Cellulari (SingleR)",
      subtitle = "Riferimento: Human Primary Cell Atlas"
    ) +
    theme_bw()
)
dev.off()
cat("Grafico UMAP delle annotazioni salvato in:", output_plot_path, "\n")


# 9. Creazione e stampa della tabella di concordanza
concordance_table <- table(
  Cluster = scrna_obj$seurat_clusters,
  Predicted_Cell_Type = scrna_obj$singleR_labels
)
cat("\nConcordanza tra i Cluster di Seurat e le Annotazioni di SingleR:\n")
print(concordance_table)


# ---TASK 8: TISSUE OF ORIGIN---

# Calculate and plot the cell type proportions
cell_type_counts <- table(scrna_obj$singleR_labels)
cell_type_proportions <- prop.table(cell_type_counts) * 100
cell_type_df <- as.data.frame(cell_type_proportions)
colnames(cell_type_df) <- c("CellType", "Percentage")

ggplot(cell_type_df, aes(x = reorder(CellType, -Percentage), y = Percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "Predicted Cell Composition of the Sample",
    x = "Predicted Cell Type",
    y = "Percentage of Total Cells (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))


# Visualize canonical markers (this was saved to a file in the script, here we regenerate it)
FeaturePlot(scrna_obj, features = c("PTPRC", "CD3D", "MS4A1", "CD14"), pt.size = 0.4, order = TRUE)
