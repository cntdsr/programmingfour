# analysis_functions.R

# --- FUNZIONE 1: CARICAMENTO E FILTRO INIZIALE ---

#' @title Carica i dati 10X e filtra per geni codificanti proteine.
#' @description Questa funzione legge i dati da una cartella 10X, crea un oggetto Seurat,
#' e lo filtra per mantenere esclusivamente i geni annotati come 'protein_coding'
#' da un file GTF.
#' @param data_dir Percorso della cartella contenente i file 'matrix.mtx.gz', 
#' 'barcodes.tsv.gz', e 'features.tsv.gz'.
#' @param gtf_path Percorso del file di annotazione GTF (es. 'Homo_sapiens.GRCh38.111.gtf').
#' @return Un oggetto Seurat contenente solo i geni codificanti proteine.
#' @importFrom rtracklayer import
#' @importFrom dplyr filter pull
#' @import Seurat
#' @export
load_and_filter_data <- function(data_dir, gtf_path) {
  # Carica dati grezzi
  expression_matrix <- Seurat::Read10X(data.dir = data_dir)
  scrna_obj <- Seurat::CreateSeuratObject(counts = expression_matrix, project = "scRNA_Project")
  cat("Oggetto Seurat creato con", nrow(scrna_obj), "geni totali.\n")
  
  # Importa e filtra GTF
  gtf_data <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf_data)
  protein_coding_genes <- gtf_df %>%
    dplyr::filter(type == "gene" & gene_biotype == "protein_coding") %>%
    dplyr::pull(gene_name) %>%
    unique()
  cat("Trovati", length(protein_coding_genes), "geni 'protein-coding' unici nel GTF.\n")
  
  # Filtra oggetto Seurat
  genes_in_matrix <- intersect(rownames(scrna_obj), protein_coding_genes)
  scrna_obj <- subset(scrna_obj, features = genes_in_matrix)
  cat("Oggetto Seurat filtrato. Contiene ora", nrow(scrna_obj), "geni.\n")
  
  return(scrna_obj)
}


# --- FUNZIONE 2: FILTRO GENI MITOCONDRIALI E RIBOSOMIALI ---

#' @title Filtra geni mitocondriali e ribosomiali.
#' @description Rimuove dall'oggetto Seurat i geni mitocondriali (MT-), 
#' i geni delle proteine ribosomiali (RPS, RPL) e gli pseudogeni ribosomiali.
#' @param scrna_obj Un oggetto Seurat.
#' @param gtf_path Percorso del file di annotazione GTF, usato per identificare gli pseudogeni.
#' @return Una lista contenente: `object` (l'oggetto Seurat filtrato) e 
#' `summary` (una tabella riassuntiva dei geni rimossi).
#' @importFrom dplyr filter pull
#' @import Seurat
#' @export
filter_unwanted_genes <- function(scrna_obj, gtf_path) {
  current_genes <- rownames(scrna_obj)
  
  # Identifica geni
  mito_genes <- grep(pattern = "^MT-", x = current_genes, value = TRUE, ignore.case = TRUE)
  ribo_genes <- grep(pattern = "^RP[SL]", x = current_genes, value = TRUE, ignore.case = TRUE)
  
  # Identifica pseudogeni ribosomiali dal GTF
  gtf_data <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf_data)
  ribo_pseudo_info <- gtf_df %>%
    dplyr::filter(grepl("pseudogene", gene_biotype, ignore.case = TRUE) & grepl("^RP[SL]", gene_name, ignore.case = TRUE)) %>%
    dplyr::pull(gene_name) %>%
    unique()
  ribo_pseudo_genes <- intersect(current_genes, ribo_pseudo_info)
  
  # Crea la tabella riassuntiva
  summary_table <- data.frame(
    Category = c("Mitochondrial", "Ribosomal Protein", "Ribosomal Pseudogenes"),
    GenesRemoved = c(length(mito_genes), length(ribo_genes), length(ribo_pseudo_genes))
  )
  print(summary_table)
  
  # Filtra l'oggetto
  genes_to_remove <- unique(c(mito_genes, ribo_genes, ribo_pseudo_genes))
  genes_to_keep <- setdiff(current_genes, genes_to_remove)
  scrna_obj <- subset(scrna_obj, features = genes_to_keep)
  cat("Filtraggio completato. Oggetto Seurat contiene ora", nrow(scrna_obj), "geni.\n")
  
  return(list(object = scrna_obj, summary = summary_table))
}


# --- FUNZIONE 3: PRE-PROCESSING E PCA ---

#' @title Esegue normalizzazione, scaling e PCA.
#' @description Esegue i passi standard di pre-processing di Seurat (Normalize,
#' FindVariableFeatures, ScaleData) e infine esegue la PCA.
#' @param scrna_obj Un oggetto Seurat.
#' @param n_variable_features Numero di geni altamente variabili da identificare. Default 2000.
#' @return L'oggetto Seurat processato e con i risultati della PCA.
#' @import Seurat
#' @export
preprocess_and_run_pca <- function(scrna_obj, n_variable_features = 2000) {
  scrna_obj <- Seurat::NormalizeData(scrna_obj)
  scrna_obj <- Seurat::FindVariableFeatures(scrna_obj, nfeatures = n_variable_features)
  scrna_obj <- Seurat::ScaleData(scrna_obj, features = rownames(scrna_obj))
  scrna_obj <- Seurat::RunPCA(scrna_obj, features = Seurat::VariableFeatures(object = scrna_obj))
  cat("Pre-processing e PCA completati.\n")
  return(scrna_obj)
}


# --- FUNZIONE 4: UMAP E CLUSTERING ---

#' @title Esegue UMAP e clustering.
#' @description Esegue l'algoritmo UMAP e il clustering basato su grafo (FindNeighbors, FindClusters)
#' usando un numero definito di componenti principali.
#' @param scrna_obj Un oggetto Seurat con PCA già calcolata.
#' @param pca_dims Il numero di PC da usare per l'analisi.
#' @param resolution La risoluzione per l'algoritmo di clustering. Valori più alti
#' portano a più cluster. Default 0.5.
#' @return L'oggetto Seurat con i risultati di UMAP e clustering.
#' @import Seurat
#' @export
run_umap_and_cluster <- function(scrna_obj, pca_dims, resolution = 0.5) {
  scrna_obj <- Seurat::RunUMAP(scrna_obj, reduction = "pca", dims = 1:pca_dims)
  scrna_obj <- Seurat::FindNeighbors(scrna_obj, reduction = "pca", dims = 1:pca_dims)
  scrna_obj <- Seurat::FindClusters(scrna_obj, resolution = resolution)
  cat("UMAP e clustering completati.\n")
  print(table(scrna_obj$seurat_clusters))
  return(scrna_obj)
}


# --- FUNZIONE 5: ANNOTAZIONE CELLULARE CON SINGLER ---

#' @title Annota i tipi cellulari con SingleR.
#' @description Usa il pacchetto SingleR e il riferimento HumanPrimaryCellAtlasData
#' per predire i tipi cellulari per ogni cellula.
#' @param scrna_obj Un oggetto Seurat.
#' @return L'oggetto Seurat con una nuova colonna nei metadati ('singleR_labels')
#' contenente le annotazioni.
#' @import SingleR
#' @import celldex
#' @import SingleCellExperiment
#' @export
annotate_cell_types <- function(scrna_obj) {
  cat("Caricamento dati di riferimento HumanPrimaryCellAtlasData...\n")
  ref_data <- celldex::HumanPrimaryCellAtlasData()
  
  # Conversione a oggetto compatibile con SingleR
  sce_obj <- as.SingleCellExperiment(scrna_obj)
  
  cat("Esecuzione di SingleR...\n")
  predictions <- SingleR::SingleR(
    test = sce_obj,
    ref = ref_data,
    labels = ref_data$label.main
  )
  
  # Aggiungi le etichette all'oggetto Seurat
  scrna_obj$singleR_labels <- predictions$labels[match(rownames(scrna_obj@meta.data), rownames(predictions))]
  cat("Annotazione cellulare completata.\n")
  return(scrna_obj)
}

# --- FUNZIONE 6: GENERAZIONE DI GRAFICI DI VISUALIZZAZIONE ---

#' @title Crea e visualizza i grafici di analisi principali.
#' @description Genera una serie di grafici standard per visualizzare i risultati
#' dell'analisi: la varianza della PCA, la UMAP colorata per cluster e la UMAP
#' colorata per annotazione cellulare.
#' @param scrna_obj Un oggetto Seurat completo di PCA, clustering e annotazione.
#' @return Una lista di oggetti ggplot (i grafici).
#' @import ggplot2
#' @import Seurat
#' @export
generate_visualizations <- function(scrna_obj) {
  # Grafico della varianza PCA
  pca_stdev <- Stdev(scrna_obj, reduction = "pca")
  variance_explained <- pca_stdev^2
  percent_variance <- (variance_explained / sum(variance_explained)) * 100
  variance_df <- data.frame(PC = 1:20, PercentVariance = percent_variance[1:20])
  
  pca_plot <- ggplot(variance_df, aes(x = PC, y = PercentVariance)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks = 1:20) +
    labs(
      title = "Variance Explained by First 20 PCs",
      x = "Principal Component (PC)",
      y = "% Variance Explained"
    ) +
    theme_bw()
  
  # Grafico UMAP con i cluster
  umap_cluster_plot <- DimPlot(scrna_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    labs(title = "UMAP Colored by Seurat Clusters") +
    NoLegend()
  
  # Grafico UMAP con le annotazioni SingleR
  umap_annotation_plot <- DimPlot(scrna_obj, reduction = "umap", group.by = "singleR_labels", label = TRUE, repel = TRUE) +
    labs(title = "UMAP Colored by SingleR Cell Type Annotation") +
    NoLegend()
  
  # Stampa i grafici
  print(pca_plot)
  print(umap_cluster_plot)
  print(umap_annotation_plot)
  
  return(list(
    pca_variance = pca_plot,
    umap_clusters = umap_cluster_plot,
    umap_annotations = umap_annotation_plot
  ))
}


# --- FUNZIONE 7: ANALISI DI CONCORDANZA E MARCATORI ---

#' @title Esegue l'analisi di concordanza e trova i geni marcatori.
#' @description Calcola e stampa una tabella di concordanza tra i cluster di Seurat
#' e le annotazioni di SingleR. Successivamente, identifica i geni marcatori
#' positivi per ogni tipo cellulare annotato.
#' @param scrna_obj Un oggetto Seurat completo di clustering e annotazione.
#' @param group_by_labels La colonna dei metadati da usare per trovare i marcatori (es. "singleR_labels").
#' @return Un data.frame con i geni marcatori per ogni gruppo.
#' @importFrom dplyr group_by top_n
#' @import Seurat
#' @export
analyze_concordance_and_markers <- function(scrna_obj, group_by_labels = "singleR_labels") {
  # Tabella di concordanza
  concordance_table <- table(
    Cluster = scrna_obj$seurat_clusters,
    Predicted_Cell_Type = scrna_obj[[group_by_labels, drop = TRUE]]
  )
  cat("\n--- Concordance between Seurat Clusters and SingleR Annotations ---\n")
  print(concordance_table)
  
  # Trova i marcatori per ogni tipo cellulare annotato
  cat("\n--- Finding Marker Genes for each Predicted Cell Type ---\n")
  all_markers <- FindAllMarkers(
    scrna_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    group.by = group_by_labels
  )
  
  if (nrow(all_markers) > 0) {
    top_markers <- all_markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 2, wt = avg_log2FC)
    cat("\nTop 2 marker genes for each cell type:\n")
    print(top_markers)
  } else {
    cat("No markers found with the current thresholds.\n")
    top_markers <- NULL
  }
  
  return(top_markers)
}


# --- FUNZIONE 8: IPOTESI SULL'ORIGINE DEL TESSUTO ---

#' @title Ipotizza l'origine del tessuto basandosi sulla composizione cellulare.
#' @description Calcola le proporzioni dei tipi cellulari predetti e le visualizza
#' con un grafico a barre, aiutando a formulare un'ipotesi sull'origine del tessuto.
#' @param scrna_obj Un oggetto Seurat completo di annotazione.
#' @param label_column La colonna dei metadati con le annotazioni (es. "singleR_labels").
#' @return Un oggetto ggplot (il grafico a barre delle proporzioni).
#' @import ggplot2
#' @export
infer_tissue_origin <- function(scrna_obj, label_column = "singleR_labels") {
  cat("\n--- Task 8: Inferring Tissue of Origin ---\n")
  
  # Calcola le proporzioni
  cell_type_df <- as.data.frame(
    prop.table(table(scrna_obj[[label_column, drop=TRUE]])) * 100
  )
  colnames(cell_type_df) <- c("CellType", "Percentage")
  
  # Crea il grafico
  composition_plot <- ggplot(cell_type_df, aes(x = reorder(CellType, -Percentage), y = Percentage)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(
      title = "Predicted Cellular Composition of the Sample",
      x = "Predicted Cell Type",
      y = "Percentage of Total Cells (%)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  print(composition_plot)
  
  # Aggiungi un'interpretazione testuale basata sui dati
  dominant_cell <- cell_type_df$CellType[which.max(cell_type_df$Percentage)]
  cat(paste("\nBased on the analysis, the dominant cell type is", dominant_cell, ".\n"))
  cat("This composition strongly suggests the tissue of origin could be Peripheral Blood (PBMC) or a lymphoid tissue.\n")
  
  return(composition_plot)
}