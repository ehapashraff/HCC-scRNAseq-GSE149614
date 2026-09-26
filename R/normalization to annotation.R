pbmc = NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

pbmc = FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000)

VariableFeaturePlot(pbmc)

top10 = head(
  VariableFeatures(pbmc),
  10)

top10

pbmc = ScaleData(
  pbmc)

pbmc = RunPCA(
  pbmc,
  features = VariableFeatures(pbmc))

ElbowPlot(pbmc)

pbmc = FindNeighbors(
  pbmc,
  dims = 1:10)

pbmc = FindClusters(
  pbmc,
  resolution = 0.5)

pbmc = RunUMAP(
  pbmc,
  dims = 1:10)

DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE)


markers = FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

write.csv(
  markers,
  "cluster_markers.csv",
  row.names = FALSE)

saveRDS(
  pbmc,
  "HCC_GSE149614_after_clustering.rds")

# Step 15 — Cell-type annotation

celltype_annotations = c(
  "0" = "Treg_T_Cell",
  "1" = "SPP1_Macrophage",
  "2" = "Kupffer_Resident_Macrophage",
  "3" = "NK_Cytotoxic_Lymphocyte",
  "4" = "Hepatocyte_Zonated",
  "5" = "Cycling_Hepatocyte_Epithelial",
  "6" = "Hepatocyte_Zonated",
  "7" = "Hepatocyte_Zonated",
  "8" = "Endothelial",
  "9" = "Fibroblast_CAF",
  "10" = "Hepatocyte_Zonated",
  "11" = "Kupffer_Resident_Macrophage",
  "12" = "SPP1_Macrophage",
  "13" = "Smooth_Muscle_Pericyte",
  "14" = "Cycling_Proliferating",
  "15" = "Hepatocyte_Zonated")

pbmc$celltype = unname(
  celltype_annotations[
    as.character(pbmc$seurat_clusters)])

table(pbmc$celltype)

DimPlot(pbmc, group.by = "celltype",reduction = "umap",repel = TRUE)

ggsave(
  "figures/UMAP_annotated_celltypes.png",
  width = 10,
  height = 7,
  dpi = 300)

saveRDS(
  pbmc,
  "HCC_GSE149614_annotated.rds")
