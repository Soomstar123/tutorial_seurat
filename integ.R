library(Seurat)
library(dplyr)
library(patchwork)
HP.data <- ReadMtx(cells='helped/barcodes.tsv', features='helped/features.tsv.gz', mtx='helped/matrix.mtx.gz')
HPS.data <- Read10X(data.dir = "helpless/")
NAI.data <- Read10X(data.dir = "naive/")


#Create the Seurat object
HP <- CreateSeuratObject(counts = HP.data, project = "HP", min.cells = 3, min.features = 200)
HPS <- CreateSeuratObject(counts = HPS.data, project = "HPS", min.cells = 3, min.features = 200)
NAI <- CreateSeuratObject(counts = NAI.data, project = "NAI", min.cells = 3, min.features = 200)

#Quality control
HP[["percent.mt"]] <- PercentageFeatureSet(HP, pattern = "^mt-")
VlnPlot(HP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HP <- subset(HP, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)



HPS[["percent.mt"]] <- PercentageFeatureSet(HPS, pattern = "^mt-")
VlnPlot(HPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HPS <- subset(HPS, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

NAI[["percent.mt"]] <- PercentageFeatureSet(NAI, pattern = "^mt-")
VlnPlot(NAI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NAI <- subset(NAI, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)



#Normalization
HP <- NormalizeData(HP)
HPS<- NormalizeData(HPS)
NAI <- NormalizeData(NAI)

#Integration
HP <- FindVariableFeatures(HP, selection.method = "vst", nfeatures = 2000)
HPS <- FindVariableFeatures(HPS, selection.method = "vst", nfeatures = 2000)
NAI <- FindVariableFeatures(NAI, selection.method = "vst", nfeatures = 2000)

immune.list = list(Helped=HP, Helpless=HPS, Naive=NAI)
immune.features <- SelectIntegrationFeatures(object.list = immune.list)
immune.anchors <- FindIntegrationAnchors(object.list = immune.list, anchor.features = immune.features)

integrated <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(integrated) <- "integrated"


#scaling data
all.genes <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genes)


#PCA
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
print(integrated[["pca"]], dims = 1:30, nfeatures = 5)
VizDimLoadings(integrated, dims = 1:2, reduction = "pca")
ElbowPlot(integrated)
DimPlot(integrated, dim = c(1,2),reduction = "pca")


#clustring
integrated <- FindNeighbors(integrated, dims = 1:15)
integrated <- FindClusters(integrated, resolution = 0.5)
head(Idents(integrated), 5)
integrated <- RunUMAP(integrated, dims = 1:15)

DimPlot(integrated, reduction = "umap", label=TRUE, repel=TRUE)
saveRDS(integrated, "integ1.rds")

#findmarkers
integ.markers <- FindAllMarkers(integrated, logfc.threshold = 1, only.pos = T, test.use = "wilcox")
integ.markers=integ.markers%>%filter(p_val_adj < 0.01)
integ.markers=integ.markers%>%arrange(cluster,desc(avg_log2FC))
top5 <- integ.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(integ.markers,file = "integ.markerslist.csv",quote = F,row.names = F)



