############ XENIUM ANALYSIS WITH SEURAT ############
###############   Patient ID UCI-604 ################

install.packages("remotes")
remotes::install_github("satijalab/seurat")

library(Seurat)
library(future)
plan("multisession", workers = 24)

library(dplyr)

remotes::install_github("satijalab/seurat", "seurat5", force = TRUE)
library(Seurat, options(Seurat.object.assay.version = 'v5'))

remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
remotes::install_github("bnprks/BPCells", quiet = TRUE)

install.packages("rlang") ## Note: some error messages

devtools::install_github("dmcable/spacexr", build_vignettes = FALSE) ### Note: some error messages
library(spacexr) ### Error in library(spacexr) : there is no package called ‘spacexr’


############## Load the Xenium data one by one ############
# IPSILATERAL
ipsi.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/IPSILATERAL", fov = "ipsi")
# remove cells with 0 counts
ipsi.obj <- subset(ipsi.obj, subset = nCount_Xenium > 0)
VlnPlot(ipsi.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(ipsi.obj, fov = "ipsi", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P2-LIQ_C_renamed
P2.LIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P2-LIQ_C_renamed", fov = "P2.LIQ")
# remove cells with 0 counts
P2.LIQ.obj <- subset(P2.LIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P2.LIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P2.LIQ.obj, fov = "P2.LIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P2-LOQ_C
P2.LOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P2-LOQ_C", fov = "P2.LOQ")
# remove cells with 0 counts
P2.LOQ.obj <- subset(P2.LOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P2.LOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P2.LOQ.obj, fov = "P2.LOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P2-UIQ_C_renamed
P2.UIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P2-UIQ_C_renamed", fov = "P2.UIQ")
# remove cells with 0 counts
P2.UIQ.obj <- subset(P2.UIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P2.UIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P2.UIQ.obj, fov = "P2.UIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P2-UOQ_C_renamed
P2.UOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P2-UOQ_C_renamed", fov = "P2.UOQ")
# remove cells with 0 counts
P2.UOQ.obj <- subset(P2.UOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P2.UOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P2.UOQ.obj, fov = "P2.UOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)


# P3.A.LIQ_C
P3.A.LIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-A-LIQ", fov = "P3.A.LIQ")
# remove cells with 0 counts
P3.A.LIQ.obj <- subset(P3.A.LIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.A.LIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.A.LIQ.obj, fov = "P3.A.LIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.A.LOQ_C
P3.A.LOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-A-LOQ", fov = "P3.A.LOQ")
# remove cells with 0 counts
P3.A.LOQ.obj <- subset(P3.A.LOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.A.LOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.A.LOQ.obj, fov = "P3.A.LOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3-A-UIQ_C_renamed
P3.A.UIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-A-UIQ_C_renamed", fov = "P3.A.UIQ")
# remove cells with 0 counts
P3.A.UIQ.obj <- subset(P3.A.UIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.A.UIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.A.UIQ.obj, fov = "P3.A.UIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3-A-UOQ_C_renamed
P3.A.UOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-A-UOQ_C_renamed", fov = "P3.A.UOQ")
# remove cells with 0 counts
P3.A.UOQ.obj <- subset(P3.A.UOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.A.UOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.A.UOQ.obj, fov = "P3.A.UOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.M.LIQ_C
P3.M.LIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-M-LIQ", fov = "P3.M.LIQ")
# remove cells with 0 counts
P3.M.LIQ.obj <- subset(P3.M.LIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.M.LIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.M.LIQ.obj, fov = "P3.M.LIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.M.LOQ_C
P3.M.LOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-M-LOQ", fov = "P3.M.LOQ")
# remove cells with 0 counts
P3.M.LOQ.obj <- subset(P3.M.LOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.M.LOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.M.LOQ.obj, fov = "P3.M.LOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.M.UIQ_C
P3.M.UIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-M-UIQ", fov = "P3.M.UIQ")
# remove cells with 0 counts
P3.M.UIQ.obj <- subset(P3.M.UIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.M.UIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.M.UIQ.obj, fov = "P3.M.UIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.M.UOQ_C
P3.M.UOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-M-UOQ", fov = "P3.M.UOQ")
# remove cells with 0 counts
P3.M.UOQ.obj <- subset(P3.M.UOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.M.UOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.M.UOQ.obj, fov = "P3.M.UOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.P.LOQ_C
P3.P.LOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-P-LOQ", fov = "P3.P.LOQ")
# remove cells with 0 counts
P3.P.LOQ.obj <- subset(P3.P.LOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.P.LOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.P.LOQ.obj, fov = "P3.P.LOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.P.LIQ_C
P3.P.LIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-P-LIQ_D_renamed", fov = "P3.P.LIQ")
# remove cells with 0 counts
P3.P.LIQ.obj <- subset(P3.P.LIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.P.LIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.P.LIQ.obj, fov = "P3.P.LIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.P.UOQ_C
P3.P.UOQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-P-UOQ", fov = "P3.P.UOQ")
# remove cells with 0 counts
P3.P.UOQ.obj <- subset(P3.P.UOQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.P.UOQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.P.UOQ.obj, fov = "P3.P.UOQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# P3.P.UIQ_C
P3.P.UIQ.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/P3-P-UIQ", fov = "P3.P.UIQ")
# remove cells with 0 counts
P3.P.UIQ.obj <- subset(P3.P.UIQ.obj, subset = nCount_Xenium > 0)
VlnPlot(P3.P.UIQ.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(P3.P.UIQ.obj, fov = "P3.P.UIQ", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)

# TUMOR_C
TUMOR.obj <- LoadXenium("D:/XENIUM_DATA/xenium_prerelease_kessenbrock_hbreast_Jan05/TUMOR_D_renamed", fov = "TUMOR")
# remove cells with 0 counts
TUMOR.obj <- subset(TUMOR.obj, subset = nCount_Xenium > 0)
VlnPlot(TUMOR.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(TUMOR.obj, fov = "TUMOR", molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 20000)


################### MERGE DATA INTO ONE XENIUM OBJECT #######################
## OMIT SAMPLES IPSILATERAL & TUMOR FOR NOW TO FOCUS ONLY ON DEEP DIVE ##
UCI604.Xenium <- merge(P2.LIQ.obj, y = c(P2.LOQ.obj, P2.UIQ.obj, P2.UOQ.obj, P3.A.LIQ.obj, P3.A.LOQ.obj, 
                                         P3.A.UIQ.obj, P3.A.UOQ.obj, P3.M.LIQ.obj, P3.M.LOQ.obj, P3.M.UIQ.obj,
                                         P3.M.UOQ.obj, P3.P.LIQ.obj, P3.P.LOQ.obj, P3.P.UIQ.obj, P3.P.UOQ.obj),
                       add.cell.ids = c("P2LIQ", "P2LOQ", "P2UIQ","P2UOQ","P3ALIQ","P3ALOQ",
                                        "P3AUIQ","P3AUOQ","P3MLIQ","P3MLOQ","P3MUIQ",
                                        "P3MUOQ","P3PLIQ","P3PLOQ","P3PUIQ","P3PUOQ"))
UCI604.Xenium

# Save Seurat object #
saveRDS(UCI604.Xenium, file = "/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen/UCI604.Xenium.rds")

# Load in data #
UCI604.Xenium <- readRDS(file ="/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen/UCI604.Xenium.rds")

# visualize the expression level of some key layer marker genes at the per-cell level using ImageFeaturePlot(), 
# which is analogous to the FeaturePlot() function for visualizing expression on a 2D embedding. 
# We manually adjust the max.cutoff for each gene to roughly the 90th percentile (which can be specified with max.cutoff='q90')
# of it’s count distribution to improve contrast.
ImageFeaturePlot(UCI604.Xenium, fov = "P2.LIQ", features = c("KRT14", "MYLK", "KRT8", "LUM"), 
                 max.cutoff = c(3, 5, 3, 10), size = 0.75, cols = c("white", "red"), axes = TRUE)

ImageFeaturePlot(P2.LIQ.obj, fov = "P2.LIQ", features = c("KRT14", "MYLK", "KRT8", "LUM"), 
                 max.cutoff = c(3, 5, 3, 10), size = 0.75, cols = c("white", "red"), axes = TRUE) # just to test if cells are correct

ImageFeaturePlot(UCI604.Xenium, fov = "P2.LOQ", features = c("KRT14", "MYLK", "KRT8", "LUM"), 
                 max.cutoff = c(3, 5, 3, 10), size = 0.75, cols = c("white", "red"), axes = TRUE)
ImageFeaturePlot(UCI604.Xenium, fov = "P3.P.LIQ", features = c("KRT14", "MYLK", "KRT8", "LUM"), 
                 max.cutoff = c(3, 5, 3, 10), size = 0.75, cols = c("white", "red"), axes = TRUE)

#### We can zoom in on a chosen area with the Crop() function. 
#### Once zoomed-in, we can visualize cell segmentation boundaries along with individual molecules.
cropped.coords <- Crop(UCI604.Xenium[["P2.LIQ"]], x = c(1000, 2000), y = c(2000, 3000), coords = "plot") ## doesnt work


UCI604.Xenium[["P2.LIQ.zoom"]] <- cropped.coords ## Error --> doesnt work

#### visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(UCI604.Xenium[["P2.LIQ.zoom"]]) <- "segmentation" ## doesnt work
ImageDimPlot(UCI604.Xenium, fov = "P2.LIQ.zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("KRT14", "MYLK", "KRT8", "LUM"), nmols = 10000) ## Error --> doesnt work

## issue in code below was resolved by uninstalling package Matrix and reinstalling latest version
#remove.packages("Matrix")
#install.packages("Matrix")
#library("Matrix")

### Next, we use SCTransform for normalization followed by standard dimensionality reduction and clustering.
UCI604.Xenium <- SCTransform(UCI604.Xenium, assay = "Xenium")
UCI604.Xenium <- RunPCA(UCI604.Xenium, npcs = 30, features = rownames(UCI604.Xenium))
UCI604.Xenium <- RunUMAP(UCI604.Xenium, dims = 1:30)
UCI604.Xenium <- FindNeighbors(UCI604.Xenium, reduction = "pca", dims = 1:30)
UCI604.Xenium <- FindClusters(UCI604.Xenium, resolution = 0.3)

DimPlot(UCI604.Xenium)

#Highlight specific features
FeaturePlot(UCI604.Xenium, features = c("KRT14", "MYLK", "KRT8","ANKRD30A", "SLPI", "LUM"))

# Color cell positions colored by cluster labels
ImageDimPlot(UCI604.Xenium, fov="P2.LIQ", cols = "polychrome", size = 0.75)

UCI604.Xenium.res.0.3.markers <- FindAllMarkers(object = UCI604.Xenium, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UCI604.Xenium.res.0.3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10.UCI604.Xenium.res.0.3.markers
DoHeatmap(object = UCI604.Xenium, features = top10.UCI604.Xenium.res.0.3.markers$gene, raster=FALSE) + NoLegend()

write.csv(UCI604.Xenium.res.0.3.markers, "D:/Xenium_ANALYSIS_SEURAT/UCI604.Xenium.res.0.3.markers.csv")

VlnPlot(UCI604.Xenium, features = c('PTN','KRT5','KRT14'))

#set active ident to column that I want to rename --> in this case, numbered clusters
Idents(UCI604.Xenium) <- "seurat_clusters"

#Rename clusters to assign cell types
new_names <- (c("Fibroblasts", "Basal", "Vascular", "LumSec", "LumHR", "T cells", "T cells", "Myeloid", "Fibroblasts", "Adipocytes", "Pericytes", "Lymphatic" ,"B cells", "Mast cells", "T cells", "Vascular", "mixed epithelial"))
names(new_names) <- levels(UCI604.Xenium)
UCI604.Xenium <- RenameIdents(UCI604.Xenium, new_names)
UCI604.Xenium@meta.data$new_names <- Idents(UCI604.Xenium)

DimPlot(UCI604.Xenium, label = TRUE)
DimPlot(UCI604.Xenium, group.by = "seurat_clusters", label = TRUE)

### Add a new meta data column that contains information on which region the cells came from ######
######################################################################################################
# First, check which unique region names had been added in the first step when I merged the Xenium objects, using add.cell.ids
unique(sapply(X = strsplit(colnames(UCI604.Xenium), split = "_"), FUN = "[", 1))
# Then, create a vector that contains region information for all cells in the combined Xenium object
regions <- sapply(X = strsplit(colnames(UCI604.Xenium), split = "_"), FUN = "[", 1)
# Last, add a meta data column termed "region" with the information stored in the vector
UCI604.Xenium$region <- regions

# DimPlot to confirm meta data with region information was successfully added
DimPlot(UCI604.Xenium, label = TRUE, split.by = "region")

# How does cluster membership vary by condition?
table(Idents(UCI604.Xenium), UCI604.Xenium$region)
prop.table(table(Idents(UCI604.Xenium), UCI604.Xenium$region), margin = 2)

Fraction <- table(Idents(UCI604.Xenium), UCI604.Xenium$region)
Fraction <- as.data.frame(Fraction)
Fraction$Var1 <- as.character(Fraction$Var1)

colourcount = length(unique(UCI604.Xenium$region))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

### stacked barplot
library(ggplot2)
library(RColorBrewer)
ggplot(Fraction, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.7) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourcount)) +
  labs(fill="Cell Proportions")

UCI604.Xenium.Cell_type_fractions <- table(Idents(UCI604.Xenium), UCI604.Xenium$region)

write.csv(UCI604.Xenium.Cell_type_fractions, "D:/Xenium_ANALYSIS_SEURAT/UCI604.Xenium.cell_type_fractions.csv")

########### Label transfer from HBCA2 Patient UCI-604 object ########
### Note: this is the traditional Label transfer pipeline, not the updated one for Xenium specifically (RCTD, see below)
# Load in HBCA2  object
UCI_604.scRNASeq <- readRDS('/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen/UCI_604_scRNASeq.rds')

UCI_604.scRNASeq

# Define anchors from HBCA1 object and transfer them onto UCI604 data
# Set default assay to RNA for both objects
DefaultAssay(UCI_604.scRNASeq) <- "SCT"
DefaultAssay(UCI604.Xenium) <- "SCT"

Celltype_anchors <- FindTransferAnchors(reference = UCI_604.scRNASeq, query = UCI604.Xenium,
                                        dims = 1:30, reference.reduction = "pca")


predictions <- TransferData(anchorset= Celltype_anchors, refdata = UCI_604.scRNASeq$cell_type, dims = 1:30)
UCI604.Xenium <- AddMetaData(UCI604.Xenium, metadata = predictions)

DimPlot(UCI604.Xenium, reduction = "umap", group.by = "predicted.id", label = TRUE)
DimPlot(UCI604.Xenium, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(UCI604.Xenium, reduction = "umap", group.by = "new_names", label = TRUE)



### Run Robust Cell Type Decomposition, a computational approach to deconvolve spot-level data from spatial datasets, when provided with an scRNA-seq reference. RCTD has been shown to accurately annotate spatial data from a variety of technologies, including SLIDE-seq, Visium, and the 10x Xenium in-situ spatial platform.
#To run RCTD, we first install the spacexr package from GitHub which implements RCTD.

# install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

devtools::install_github("dmcable/spacexr", build_vignettes = FALSE) ## doesn't work somehow already here (and all below)

library(spacexr)

query.counts <- GetAssayData(UCI604.Xenium, assay = "Xenium", slot = "counts")[, Cells(UCI604.Xenium[["crop"]])]
coords <- GetTissueCoordinates(UCI604.Xenium[["crop"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# As reference, we want to use the scRNA-Seq fixed Seurat data from the same patient (Patient ID: UCI-604)
UCI_604.scRNASeq <- readRDS("/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen/UCI_604_scRNASeq.rds")
UCI_604.scRNASeq <- UpdateSeuratObject(UCI_604.scRNASeq)

Idents(UCI_604.scRNASeq) <- "subclass"
# remove CR cells because there aren't enough of them for annotation ### Note Maren: Not sure if I need this?
UCI_604.scRNASeq <- subset(UCI_604.scRNASeq, subset = subclass != "CR")
counts <- GetAssayData(UCI_604.scRNASeq, assay = "RNA", slot = "counts")
cluster <- as.factor(UCI_604.scRNASeq$subclass)
names(cluster) <- colnames(UCI_604.scRNASeq)
nUMI <- UCI_604.scRNASeq$nCount_RNA
names(nUMI) <- colnames(UCI_604.scRNASeq)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
UCI604.Xenium$predicted.celltype <- annotations
keep.cells <- Cells(UCI604.Xenium)[!is.na(UCI604.Xenium$predicted.celltype)]
UCI604.Xenium <- subset(UCI604.Xenium, cells = keep.cells)

##### In Seurat v5, we introduce support for ‘niche’ analysis of spatial data, which demarcates regions of tissue (‘niches’), each of which is defined by a different composition of spatially adjacent cell types. 
## Inspired by methods in Goltsev et al, Cell 2018 and He et al, NBT 2022, we consider the ‘local neighborhood’ for each cell - consisting of its k.neighbor spatially closest neighbors, and count the occurrences of each cell type present in this neighborhood. 
## We then use k-means clustering to group cells that have similar neighborhoods together, into spatial niches.

#We call the BuildNicheAssay function from within Seurat to construct a new assay called niche containing the cell type composition spatially neighboring each cell. A metadata column called niches is also returned, which contains cluster assignments based on the niche assay.
UCI604.Xenium <- BuildNicheAssay(object = UCI604.Xenium, fov = "crop", group.by = "new_names",
                                 niches.k = 5, neighbors.k = 30) ### Note: Error, fov appears to be indicated wrong (also when I try with an actual fov region, eg "P2.LIQ")

celltype.plot <- ImageDimPlot(UCI604.Xenium, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type") ## Note: While we haven't run RCTD, we can also group.by "new_names", the manually annotated cell type colum that I added

niche.plot <- ImageDimPlot(UCI604.Xenium, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

table(UCI604.Xenium$predicted.celltype, UCI604.Xenium$niches)