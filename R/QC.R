library(Seurat)

pbmc = CreateSeuratObject(
  counts = counts,
  project = "HCC_GSE149614",
  min.cells = 3,
  min.features = 200)

pbmc[["percent.mt"]] = PercentageFeatureSet(
  pbmc,
  pattern = "^MT-")

VlnPlot(
  pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3)

summary(pbmc$nFeature_RNA)
summary(pbmc$nCount_RNA)
summary(pbmc$percent.mt)

pbmc = subset(
  pbmc,
  subset = nFeature_RNA >= 200 &
    nFeature_RNA <= 8000 &
    nCount_RNA >= 200 &
    percent.mt <= 10)

VlnPlot(
  pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3)

pbmc
