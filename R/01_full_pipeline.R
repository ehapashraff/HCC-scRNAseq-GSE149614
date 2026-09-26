# ============================================================================
# HCC scRNA-seq analysis — GSE149614
# Lu et al. 2022, Nature Communications 13:4594
#
# Purpose:
#   Reproduce the core single-cell workflow from the published HCC atlas and
#   extend it with modern Seurat, SingleR/celldex, GO enrichment, and
#   DoRothEA/decoupleR TF-activity analysis.
#
# Hardware-aware default:
#   TARGET_CELLS = 20000, stratified across original samples.
#   Set USE_FULL_DATA = TRUE only on a machine with substantially more RAM.
#
# Run this script from the repository root:
#   source("R/01_full_pipeline.R")
# ============================================================================

options(stringsAsFactors = FALSE)
set.seed(149614)

# ------------------------------- CONFIG -------------------------------------
USE_FULL_DATA <- FALSE
TARGET_CELLS <- 20000L
COUNT_CHUNK_NOT_USED <- TRUE  # kept as a documented reminder: fread select is used
QC_MIN_FEATURES <- 200L
QC_MAX_FEATURES <- 8000L
QC_MIN_UMI <- 200L
QC_MAX_UMI <- Inf
QC_MAX_MT <- 10
N_HVG <- 3000L
N_PCS <- 30L
CLUSTER_RESOLUTION <- 0.6
DE_MAX_CELLS_PER_IDENT <- 1000L
RUN_SINGLER <- TRUE
RUN_TF_ACTIVITY <- TRUE
RUN_GO_ENRICHMENT <- TRUE
DOWNLOAD_DATA <- TRUE

# Optional analyses intentionally not part of the default run because they are
# considerably heavier and require additional version-specific setup:
RUN_CELLCHAT <- FALSE
RUN_INFERCNV <- FALSE

# ---------------------------- REPOSITORY PATHS -------------------------------
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
DATA_DIR <- file.path(ROOT, "data", "raw")
RESULT_DIR <- file.path(ROOT, "results")
FIG_DIR <- file.path(ROOT, "figures")
DIRS <- c(DATA_DIR, RESULT_DIR, FIG_DIR,
          file.path(FIG_DIR, "QC"), file.path(FIG_DIR, "UMAP"),
          file.path(FIG_DIR, "DE"), file.path(FIG_DIR, "enrichment"),
          file.path(FIG_DIR, "TF_activity"))
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ PACKAGES ------------------------------------
cran_pkgs <- c("data.table", "Matrix", "Seurat", "ggplot2", "patchwork",
               "dplyr", "tidyr", "pheatmap", "R.utils")
bioc_pkgs <- c("SingleR", "celldex", "clusterProfiler", "org.Hs.eg.db", "decoupleR")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(NULL))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------- DATA ---------------------------------------
COUNT_FILE <- file.path(DATA_DIR, "GSE149614_HCC.scRNAseq.S71915.count.txt.gz")
META_FILE <- file.path(DATA_DIR, "GSE149614_HCC.metadata.updated.txt.gz")
BASE_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/"

if (DOWNLOAD_DATA) {
  if (!file.exists(COUNT_FILE)) {
    message("Downloading count matrix (~158 MB compressed)...")
    options(timeout = 3600)
    download.file(
      paste0(BASE_URL, basename(COUNT_FILE)),
      COUNT_FILE, mode = "wb", quiet = FALSE
    )
  }
  if (!file.exists(META_FILE)) {
    message("Downloading metadata...")
    options(timeout = 600)
    download.file(
      paste0(BASE_URL, basename(META_FILE)),
      META_FILE, mode = "wb", quiet = FALSE
    )
  }
}

stopifnot(file.exists(COUNT_FILE), file.exists(META_FILE))

meta <- fread(META_FILE, sep = "\t", data.table = FALSE)
meta <- as.data.frame(meta, stringsAsFactors = FALSE)
required_meta <- c("Cell", "sample", "patient", "site")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta)) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "))
}

# Normalize key metadata types.
meta$Cell <- as.character(meta$Cell)
meta$sample <- as.character(meta$sample)
meta$patient <- as.character(meta$patient)
meta$site <- as.character(meta$site)
if ("celltype" %in% colnames(meta)) meta$celltype <- as.character(meta$celltype)
if ("stage" %in% colnames(meta)) meta$stage <- as.character(meta$stage)
if ("virus" %in% colnames(meta)) meta$virus <- as.character(meta$virus)

# ------------------------ BALANCED CELL SAMPLING -----------------------------
stratified_sample_cells <- function(meta, target_n) {
  target_n <- min(as.integer(target_n), nrow(meta))
  sizes <- table(meta$sample)
  sizes_num <- as.integer(sizes)
  names(sizes_num) <- names(sizes)
  exact <- target_n * sizes_num / sum(sizes_num)
  alloc <- floor(exact)
  remainder <- target_n - sum(alloc)
  if (remainder > 0) {
    ord <- order(exact - alloc, decreasing = TRUE)
    alloc[ord[seq_len(remainder)]] <- alloc[ord[seq_len(remainder)]] + 1L
  }
  pieces <- lapply(names(alloc), function(s) {
    cells <- meta$Cell[meta$sample == s]
    n_take <- min(length(cells), alloc[[s]])
    if (n_take <= 0) character(0) else sample(cells, n_take)
  })
  unlist(pieces, use.names = FALSE)
}

if (USE_FULL_DATA) {
  cells_keep <- meta$Cell
} else {
  cells_keep <- stratified_sample_cells(meta, TARGET_CELLS)
}

message("Cells selected for analysis: ", length(cells_keep))
writeLines(cells_keep, file.path(RESULT_DIR, "selected_cells.txt"))

# -------------------------- COUNT MATRIX LOADING -----------------------------
# The GEO count file is genes x cells and has 25,712 genes x 71,915 cells.
# We use data.table::fread(select=...) so only the selected cell columns are
# materialized. This is important for 16-GB RAM machines.
read_selected_counts <- function(count_gz, cells_keep) {
  con <- gzfile(count_gz, open = "rt")
  header_line <- readLines(con, n = 1L, warn = FALSE)
  close(con)
  if (!length(header_line)) stop("Could not read the count-matrix header.")
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  cell_cols <- match(cells_keep, header)
  if (anyNA(cell_cols)) {
    missing_cells <- cells_keep[is.na(cell_cols)]
    stop("Some requested cells are absent from the count matrix: ",
         paste(head(missing_cells, 10), collapse = ", "))
  }

  # First column is the gene identifier; selected cell columns follow.
  select_cols <- c(1L, cell_cols)
  message("Reading selected count columns... this can take some time.")
  dt <- data.table::fread(
    count_gz,
    sep = "\t",
    header = FALSE,
    select = select_cols,
    data.table = TRUE,
    showProgress = TRUE,
    nThread = max(1L, parallel::detectCores(logical = FALSE) - 1L)
  )
  gene_ids <- as.character(dt[[1L]])
  mat <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(mat) <- "numeric"
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  rownames(mat) <- make.unique(gene_ids)
  colnames(mat) <- cells_keep
  rm(dt, gene_ids)
  gc()
  mat
}

counts <- read_selected_counts(COUNT_FILE, cells_keep)
message("Count matrix loaded: ", nrow(counts), " genes x ", ncol(counts), " cells")

# Match metadata exactly to count-matrix cell order.
meta_use <- meta[match(colnames(counts), meta$Cell), , drop = FALSE]
rownames(meta_use) <- meta_use$Cell
stopifnot(identical(rownames(meta_use), colnames(counts)))

# ----------------------------- SEURAT OBJECT --------------------------------
obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_use,
  min.cells = 3,
  min.features = 0,
  project = "HCC_GSE149614"
)
rm(counts, meta_use, meta)
gc()

obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

# QC summary before filtering.
p_qc_before <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
ggsave(file.path(FIG_DIR, "QC", "01_QC_before_filtering.png"),
       p_qc_before, width = 12, height = 4.5, dpi = 300)

p_scatter_qc <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(FIG_DIR, "QC", "02_QC_counts_vs_features.png"),
       p_scatter_qc, width = 6, height = 5, dpi = 300)

# Reproduce the main paper-level QC logic while keeping values explicit.
obj <- subset(
  obj,
  subset = nFeature_RNA >= QC_MIN_FEATURES &
           nFeature_RNA <= QC_MAX_FEATURES &
           nCount_RNA >= QC_MIN_UMI &
           nCount_RNA <= QC_MAX_UMI &
           percent.mt <= QC_MAX_MT
)

p_qc_after <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
ggsave(file.path(FIG_DIR, "QC", "03_QC_after_filtering.png"),
       p_qc_after, width = 12, height = 4.5, dpi = 300)

qc_counts <- data.frame(
  stage = c("selected_before_QC", "after_QC"),
  cells = c(length(cells_keep), ncol(obj))
)
write.csv(qc_counts, file.path(RESULT_DIR, "QC_cell_counts.csv"), row.names = FALSE)

write.csv(as.data.frame(obj@meta.data),
          file.path(RESULT_DIR, "cell_metadata_after_QC.csv"))

# ---------------------- NORMALIZATION / DIMENSIONALITY -----------------------
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = TRUE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG,
                            verbose = TRUE)

p_hvg <- VariableFeaturePlot(obj)
ggsave(file.path(FIG_DIR, "QC", "04_highly_variable_genes.png"),
       p_hvg, width = 7, height = 5, dpi = 300)

obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = "percent.mt",
                 verbose = TRUE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PCS, verbose = FALSE)

p_elbow <- ElbowPlot(obj, ndims = N_PCS)
ggsave(file.path(FIG_DIR, "UMAP", "01_elbow_plot.png"),
       p_elbow, width = 7, height = 5, dpi = 300)

obj <- FindNeighbors(obj, dims = 1:N_PCS, verbose = FALSE)
obj <- FindClusters(obj, resolution = CLUSTER_RESOLUTION, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:N_PCS, seed.use = 149614, verbose = FALSE)

p_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                     repel = TRUE) + NoLegend()
p_site <- DimPlot(obj, reduction = "umap", group.by = "site")
p_sample <- DimPlot(obj, reduction = "umap", group.by = "sample", label = FALSE)
p_paper <- if ("celltype" %in% colnames(obj@meta.data)) {
  DimPlot(obj, reduction = "umap", group.by = "celltype", label = FALSE)
} else NULL

ggsave(file.path(FIG_DIR, "UMAP", "02_UMAP_clusters.png"),
       p_cluster, width = 9, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "UMAP", "03_UMAP_tissue_site.png"),
       p_site, width = 9, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "UMAP", "04_UMAP_sample.png"),
       p_sample, width = 9, height = 7, dpi = 300)
if (!is.null(p_paper)) {
  ggsave(file.path(FIG_DIR, "UMAP", "05_UMAP_paper_celltype_validation.png"),
         p_paper, width = 10, height = 7, dpi = 300)
}

# Canonical markers for broad-lineage validation.
canonical_markers <- c(
  "ALB", "APOA1", "KRT19", "KRT8", "KRT18",
  "PTPRC", "CD3D", "CD3E", "TRBC2", "NKG7",
  "MS4A1", "CD79A", "LYZ", "CD68", "FCGR3A", "C1QC",
  "PECAM1", "VWF", "KDR", "COL1A1", "COL3A1", "DCN",
  "MMP9", "SPP1", "TREM2", "PPARG"
)
canonical_markers <- canonical_markers[canonical_markers %in% rownames(obj)]
if (length(canonical_markers) >= 4) {
  p_markers <- DotPlot(obj, features = canonical_markers, group.by = "seurat_clusters") +
    RotatedAxis()
  ggsave(file.path(FIG_DIR, "UMAP", "06_cluster_marker_dotplot.png"),
         p_markers, width = 14, height = 9, dpi = 300)
}

# ------------------------------ MARKERS -------------------------------------
Idents(obj) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  max.cells.per.ident = DE_MAX_CELLS_PER_IDENT,
  verbose = FALSE
)
write.csv(cluster_markers,
          file.path(RESULT_DIR, "cluster_markers.csv"), row.names = FALSE)

# --------------------------- CELL ANNOTATION --------------------------------
# SingleR is run at the cluster level, not per cell, which greatly reduces
# runtime and produces a transparent cluster-level annotation step.
if (RUN_SINGLER) {
  message("Running cluster-level SingleR annotation...")
  suppressPackageStartupMessages({
    library(SingleR)
    library(celldex)
  })

  cluster_avg <- AverageExpression(
    obj,
    assays = "RNA",
    group.by = "seurat_clusters",
    return.seurat = FALSE,
    layer = "data"
  )$RNA

  ref <- celldex::HumanPrimaryCellAtlasData()
  sr <- SingleR(
    test = as.matrix(cluster_avg),
    ref = ref,
    labels = ref$label.main
  )

  annotation_table <- data.frame(
    seurat_clusters = rownames(sr),
    SingleR_label = sr$pruned.labels,
    SingleR_label_unpruned = sr$labels,
    stringsAsFactors = FALSE
  )

  # Convert SingleR labels into stable broad categories.
  broad_from_singler <- function(x) {
    x <- tolower(ifelse(is.na(x), "", x))
    dplyr::case_when(
      grepl("t cell|nk", x) ~ "T_NK",
      grepl("b cell|plasma", x) ~ "B",
      grepl("mono|macro|dendritic|myeloid", x) ~ "Myeloid",
      grepl("hepatocyte", x) ~ "Hepatocyte",
      grepl("cholangi|epithelial", x) ~ "Cholangiocyte_Epithelial",
      grepl("endothelial", x) ~ "Endothelial",
      grepl("fibroblast|stromal|stellate", x) ~ "Fibroblast_Stromal",
      TRUE ~ "Other"
    )
  }
  annotation_table$broad_celltype <- broad_from_singler(annotation_table$SingleR_label_unpruned)
  write.csv(annotation_table,
            file.path(RESULT_DIR, "cluster_annotation_SingleR.csv"), row.names = FALSE)

  map_broad <- setNames(annotation_table$broad_celltype, annotation_table$seurat_clusters)
  obj$broad_celltype <- unname(map_broad[as.character(obj$seurat_clusters)])

  # Validation against the paper's published cell-type labels, when supplied
  # in GEO metadata. This is validation only; those labels are not used to
  # perform the annotation.
  if ("celltype" %in% colnames(obj@meta.data)) {
    validation <- table(
      inferred = obj$broad_celltype,
      published = obj$celltype,
      useNA = "ifany"
    )
    write.csv(as.data.frame(validation),
              file.path(RESULT_DIR, "annotation_validation_vs_paper_labels.csv"),
              row.names = FALSE)
  }

  p_annot <- DimPlot(obj, reduction = "umap", group.by = "broad_celltype",
                     label = TRUE, repel = TRUE)
  ggsave(file.path(FIG_DIR, "UMAP", "07_UMAP_broad_celltypes.png"),
         p_annot, width = 11, height = 8, dpi = 300)

  rm(cluster_avg, ref, sr, annotation_table)
  gc()
} else {
  obj$broad_celltype <- "Not_annotated"
}

# Reorder broad labels for easier plots.
obj$broad_celltype <- factor(
  obj$broad_celltype,
  levels = c("Hepatocyte", "Cholangiocyte_Epithelial", "Endothelial",
             "Fibroblast_Stromal", "Myeloid", "B", "T_NK", "Other",
             "Not_annotated")
)

# ----------------------- CELLULAR COMPOSITION --------------------------------
composition <- as.data.frame(table(obj$site, obj$broad_celltype))
colnames(composition) <- c("site", "broad_celltype", "cells")
composition <- composition %>%
  group_by(site) %>%
  mutate(percent = 100 * cells / sum(cells)) %>%
  ungroup()
write.csv(composition, file.path(RESULT_DIR, "cell_composition_by_site.csv"),
          row.names = FALSE)

p_comp <- ggplot(composition, aes(x = site, y = percent, fill = broad_celltype)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(x = NULL, y = "Cell proportion", fill = "Broad cell type",
       title = "Cellular composition across HCC tissue sites") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "UMAP", "08_cell_composition_by_site.png"),
       p_comp, width = 9, height = 6, dpi = 300)

# ---------------------- T/NK AND MYELOID SUBANALYSES -------------------------
run_subcluster <- function(seu, cellgroup, dims = 20L, resolution = 0.5) {
  keep <- WhichCells(seu, expression = broad_celltype == cellgroup)
  if (length(keep) < 100) {
    warning("Too few cells for ", cellgroup, "; skipping subanalysis.")
    return(NULL)
  }
  x <- subset(seu, cells = keep)
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = min(2000L, nrow(x)),
                            verbose = FALSE)
  x <- ScaleData(x, features = VariableFeatures(x), vars.to.regress = "percent.mt",
                 verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), npcs = dims, verbose = FALSE)
  x <- FindNeighbors(x, dims = 1:dims, verbose = FALSE)
  x <- FindClusters(x, resolution = resolution, verbose = FALSE)
  x <- RunUMAP(x, dims = 1:dims, seed.use = 149614, verbose = FALSE)
  x
}

if (any(obj$broad_celltype == "T_NK")) {
  tnk <- run_subcluster(obj, "T_NK", dims = 20L, resolution = 0.5)
  if (!is.null(tnk)) {
    p <- DimPlot(tnk, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
    ggsave(file.path(FIG_DIR, "UMAP", "09_TNK_subclusters.png"), p,
           width = 8, height = 6, dpi = 300)
    tnk_markers <- FindAllMarkers(tnk, only.pos = TRUE, min.pct = 0.2,
                                  logfc.threshold = 0.25, verbose = FALSE)
    write.csv(tnk_markers, file.path(RESULT_DIR, "T_NK_subcluster_markers.csv"), row.names = FALSE)
    saveRDS(tnk, file.path(RESULT_DIR, "T_NK_subanalysis.rds"))
  }
}

if (any(obj$broad_celltype == "Myeloid")) {
  myeloid <- run_subcluster(obj, "Myeloid", dims = 20L, resolution = 0.5)
  if (!is.null(myeloid)) {
    p <- DimPlot(myeloid, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
    ggsave(file.path(FIG_DIR, "UMAP", "10_myeloid_subclusters.png"), p,
           width = 8, height = 6, dpi = 300)
    my_markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.2,
                                 logfc.threshold = 0.25, verbose = FALSE)
    write.csv(my_markers, file.path(RESULT_DIR, "myeloid_subcluster_markers.csv"), row.names = FALSE)
    # The HCC paper specifically highlights MMP9/SPP1/TREM2/PPARG.
    focus_genes <- intersect(c("MMP9", "SPP1", "TREM2", "PPARG", "C1QC", "FCN1", "VCAN", "LYZ"),
                             rownames(myeloid))
    if (length(focus_genes) > 0) {
      p_focus <- FeaturePlot(myeloid, features = focus_genes, ncol = 2, order = TRUE)
      ggsave(file.path(FIG_DIR, "UMAP", "11_myeloid_HCC_focus_genes.png"), p_focus,
             width = 10, height = ceiling(length(focus_genes) / 2) * 4, dpi = 300)
    }
    saveRDS(myeloid, file.path(RESULT_DIR, "myeloid_subanalysis.rds"))
  }
}

# ----------------------- SITE-SPECIFIC EXPLORATORY DE ------------------------
# Important: cell-level DE can be anti-conservative because cells from the same
# patient are not independent replicates. A publication-grade comparison
# should use sample-level pseudobulk/edgeR or DESeq2. These results are therefore
# explicitly labelled EXPLORATORY.
run_site_de <- function(seu, cellgroup, site1 = "PT", site2 = "NTL") {
  x <- subset(seu, subset = broad_celltype == cellgroup & site %in% c(site1, site2))
  if (ncol(x) < 50) return(NULL)
  if (length(unique(x$site)) < 2) return(NULL)
  DefaultAssay(x) <- "RNA"
  Idents(x) <- x$site
  out <- tryCatch(
    FindMarkers(x, ident.1 = site1, ident.2 = site2,
                min.pct = 0.10, logfc.threshold = 0.25,
                max.cells.per.ident = DE_MAX_CELLS_PER_IDENT,
                verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  out$gene <- rownames(out)
  rownames(out) <- NULL
  out$cellgroup <- cellgroup
  out
}

for (g in intersect(c("T_NK", "Myeloid", "Hepatocyte"), unique(as.character(obj$broad_celltype)))) {
  de <- run_site_de(obj, g, "PT", "NTL")
  if (!is.null(de)) {
    write.csv(de, file.path(RESULT_DIR, paste0("EXPLORATORY_DE_PT_vs_NTL_", g, ".csv")), row.names = FALSE)
  }
}

# ----------------------------- GO ENRICHMENT -------------------------------
if (RUN_GO_ENRICHMENT) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })

  Idents(obj) <- "broad_celltype"
  broad_markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.20,
    logfc.threshold = 0.25,
    max.cells.per.ident = DE_MAX_CELLS_PER_IDENT,
    verbose = FALSE
  )
  write.csv(broad_markers, file.path(RESULT_DIR, "broad_celltype_markers.csv"), row.names = FALSE)

  enrichment_tables <- list()
  for (g in unique(broad_markers$cluster)) {
    genes <- broad_markers %>%
      filter(cluster == g, p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 200) %>%
      pull(gene) %>%
      unique()
    if (length(genes) < 10) next

    mapped <- tryCatch(
      suppressMessages(clusterProfiler::bitr(
        genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
      )),
      error = function(e) NULL
    )
    if (is.null(mapped) || nrow(mapped) < 5) next

    ego <- tryCatch(
      enrichGO(
        gene = unique(mapped$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(ego)) next
    tab <- as.data.frame(ego)
    if (!nrow(tab)) next
    tab$cellgroup <- as.character(g)
    enrichment_tables[[g]] <- tab

    p <- dotplot(ego, showCategory = 12) +
      ggtitle(paste("GO Biological Process:", g))
    safe_g <- gsub("[^A-Za-z0-9_]+", "_", as.character(g))
    ggsave(file.path(FIG_DIR, "enrichment", paste0("GO_BP_", safe_g, ".png")),
           p, width = 10, height = 7, dpi = 300)
  }

  if (length(enrichment_tables)) {
    enrichment_all <- bind_rows(enrichment_tables)
    write.csv(enrichment_all,
              file.path(RESULT_DIR, "GO_BP_enrichment_all_broad_celltypes.csv"),
              row.names = FALSE)
  }
}

# -------------------------- DoRothEA / decoupleR -----------------------------
if (RUN_TF_ACTIVITY) {
  suppressPackageStartupMessages(library(decoupleR))

  # Average expression by broad cell type keeps TF inference computationally
  # lightweight and makes the output easy to interpret.
  tf_mat <- AverageExpression(
    obj,
    assays = "RNA",
    group.by = "broad_celltype",
    return.seurat = FALSE,
    layer = "data"
  )$RNA

  dorothea <- decoupleR::get_dorothea(
    organism = "human",
    levels = c("A", "B", "C")
  )

  tf_res <- decoupleR::run_ulm(
    mat = as.matrix(tf_mat),
    network = dorothea,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = 5L
  )
  write.csv(tf_res, file.path(RESULT_DIR, "DoRothEA_ULM_TF_activity.csv"), row.names = FALSE)

  tf_wide <- tf_res %>%
    select(source, condition, score) %>%
    distinct() %>%
    pivot_wider(names_from = condition, values_from = score)

  tf_numeric <- as.matrix(tf_wide[, -1, drop = FALSE])
  rownames(tf_numeric) <- tf_wide$source
  tf_rank <- order(apply(abs(tf_numeric), 1, max, na.rm = TRUE), decreasing = TRUE)
  tf_numeric <- tf_numeric[tf_rank[seq_len(min(30L, length(tf_rank)))], , drop = FALSE]

  pheatmap::pheatmap(
    tf_numeric,
    filename = file.path(FIG_DIR, "TF_activity", "01_DoRothEA_ULM_heatmap.png"),
    width = 11,
    height = 10,
    main = "DoRothEA TF activity (ULM)"
  )
}

# -------------------------- HCC-SPECIFIC SUMMARY ----------------------------
# Focus on the genes highlighted in the paper and report their expression
# across the major broad populations. This is descriptive, not causal.
focus_genes <- intersect(
  c("MMP9", "SPP1", "TREM2", "PPARG", "MMP9", "CXCL10", "MIF", "CD74",
    "PDCD1", "CTLA4", "LAG3", "CCR7", "MS4A1"),
  rownames(obj)
)
if (length(focus_genes)) {
  avg_focus <- AverageExpression(obj, assays = "RNA", features = focus_genes,
                                 group.by = "broad_celltype", layer = "data")$RNA
  write.csv(avg_focus,
            file.path(RESULT_DIR, "HCC_focus_gene_average_expression.csv"))
}

# Save only the final metadata and small summary objects in the repository.
# The full Seurat object is intentionally NOT recommended for GitHub because
# it can exceed GitHub's practical file-size limits.
write.csv(as.data.frame(obj@meta.data),
          file.path(RESULT_DIR, "final_cell_metadata.csv"))

summary_lines <- c(
  paste0("Cells after QC: ", ncol(obj)),
  paste0("Genes: ", nrow(obj)),
  paste0("Clusters: ", length(unique(obj$seurat_clusters))),
  paste0("Sites: ", paste(sort(unique(as.character(obj$site))), collapse = ", ")),
  paste0("Broad cell types: ", paste(sort(unique(as.character(obj$broad_celltype))), collapse = ", ")),
  paste0("USE_FULL_DATA: ", USE_FULL_DATA),
  paste0("TARGET_CELLS: ", TARGET_CELLS),
  paste0("QC max mitochondrial percent: ", QC_MAX_MT)
)
writeLines(summary_lines, file.path(RESULT_DIR, "analysis_summary.txt"))

capture.output(sessionInfo(), file = file.path(RESULT_DIR, "sessionInfo.txt"))

# Optional heavy analyses are deliberately left as explicit switches rather
# than silently running potentially incompatible/version-sensitive code.
if (RUN_CELLCHAT) {
  message("RUN_CELLCHAT=TRUE: add a version-pinned CellChat workflow here after the core pipeline is validated.")
}
if (RUN_INFERCNV) {
  message("RUN_INFERCNV=TRUE: add a version-pinned inferCNV workflow here using a downloaded hg38 gene-order file and the analyzed count matrix.")
}

# Save a compact Seurat object with only the data needed for interactive review.
# This file is for local use and should remain ignored by Git.
obj@meta.data$seurat_clusters <- as.character(obj$seurat_clusters)
saveRDS(obj, file.path(RESULT_DIR, "HCC_GSE149614_analysis.rds"), compress = "xz")

message("\nAnalysis completed.")
message("Main results: ", RESULT_DIR)
message("Figures: ", FIG_DIR)
