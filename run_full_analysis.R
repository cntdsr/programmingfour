# run_full_analysis.R

# FASE 0: SETUP DELL'AMBIENTE
# -------------------------------------------------------------------

# Carica tutte le funzioni dal nostro pacchetto locale (dalla cartella R/)
# Questo comando rende disponibili funzioni come load_and_filter_data(), ecc.
devtools::load_all()

# Carichiamo anche ggplot2 per essere sicuri che le funzioni grafiche funzionino
library(ggplot2)
library(dplyr) 
cat("Pacchetto e librerie caricate con successo.\n")


# -------------------------------------------------------------------
# FASE 1: IMPOSTAZIONE DEI PARAMETRI GLOBALI
# -------------------------------------------------------------------
# Avere tutti i parametri qui rende lo script più facile da modificare e riutilizzare.

DATA_DIRECTORY   <- "."  # La directory che contiene i dati (barcodes, matrix, etc.)
GTF_FILE_PATH    <- "Homo_sapiens.GRCh38.111.gtf" # Il percorso del file di annotazione GTF
PCA_COMPONENTS   <- 13   # Numero di Principal Components da usare per UMAP e clustering
CLUSTER_RESOLUTION <- 0.5  # Risoluzione per l'algoritmo di clustering di Seurat

cat("Parametri di analisi impostati.\n")


# -------------------------------------------------------------------
# FASE 2: ESECUZIONE DELLA PIPELINE DI ANALISI
# -------------------------------------------------------------------
# Chiamiamo le nostre funzioni una dopo l'altra, passando i risultati
# di un passo come input per il passo successivo.

# TASK 1: Caricamento dati e filtro per geni codificanti proteine
cat("\n--- ESECUZIONE TASK 1: Caricamento e Filtro Geni Codificanti ---\n")
scrna_obj <- load_and_filter_data(
  data_dir = DATA_DIRECTORY,
  gtf_path = GTF_FILE_PATH
)

# TASK 3: Filtro per geni indesiderati (mitocondriali, ribosomiali)
# La task 2 sulla visualizzazione dei geni espressi è implicita in Seurat e non richiede uno step separato qui.
cat("\n--- ESECUZIONE TASK 3: Rimozione Geni Indesiderati ---\n")
filter_results <- filter_unwanted_genes(
  scrna_obj = scrna_obj,
  gtf_path = GTF_FILE_PATH
)
scrna_obj <- filter_results$object # Aggiorniamo il nostro oggetto principale

# TASK 4: Pre-processing e PCA
cat("\n--- ESECUZIONE TASK 4: Normalizzazione, Scaling e PCA ---\n")
scrna_obj <- preprocess_and_run_pca(scrna_obj = scrna_obj)

# TASK 5 & 6: UMAP e Clustering
cat("\n--- ESECUZIONE TASK 5 & 6: UMAP e Clustering ---\n")
scrna_obj <- run_umap_and_cluster(
  scrna_obj = scrna_obj,
  pca_dims = PCA_COMPONENTS,
  resolution = CLUSTER_RESOLUTION
)

# TASK 7: Annotazione dei tipi cellulari con SingleR
cat("\n--- ESECUZIONE TASK 7: Annotazione Tipi Cellulari ---\n")
scrna_obj <- annotate_cell_types(scrna_obj = scrna_obj)

cat("\n>>> Pipeline di processing dei dati completata. Inizio fase di analisi e visualizzazione. <<<\n")


# -------------------------------------------------------------------
# FASE 3: ANALISI E VISUALIZZAZIONE DEI RISULTATI
# -------------------------------------------------------------------

# Visualizzazione dei grafici principali (Varianza PCA, UMAP clusters, UMAP annotazioni)
cat("\n--- VISUALIZZAZIONE TASK 4, 6, 7: Grafici Principali ---\n")
generate_visualizations(scrna_obj = scrna_obj)

# Analisi di concordanza e ricerca dei geni marcatori
cat("\n--- ANALISI TASK 7: Concordanza e Marcatori ---\n")
top_markers_table <- analyze_concordance_and_markers(
  scrna_obj = scrna_obj,
  group_by_labels = "singleR_labels"
)

# TASK 8: Inferenza sull'origine del tessuto
cat("\n--- ANALISI TASK 8: Inferenza Origine Tessuto ---\n")
composition_plot <- infer_tissue_origin(
  scrna_obj = scrna_obj,
  label_column = "singleR_labels"
)


# -------------------------------------------------------------------
# CONCLUSIONE
# -------------------------------------------------------------------
cat("\n===================================================================\n")
cat("            SCRIPT COMPLETATO CON SUCCESSO!\n")
cat("I risultati e i grafici sono stati stampati nella console e nel pannello 'Plots'.\n")
cat("L'oggetto Seurat finale 'scrna_obj' è disponibile nell'ambiente R.\n")
cat("===================================================================\n")

