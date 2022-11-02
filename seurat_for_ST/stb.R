library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#library(hdf5r)
#library(SeuratDisk)
#library(anndata)

setwd("/data/stb/")
data_dir = "/data/stb/"

##Dataset
brain <- Load10X_Spatial(data_dir,
                         filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial",
                         slice = "tissue_hires_image.png",
                         filter.matrix = TRUE
)

##data_preprocessing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
svg("data_preprocessing.svg",width=16,height=8)
plot(plot2+plot1)
dev.off()

brain <- subset(brain, subset = nCount_Spatial > 0)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

##Gene expression visualization
svg("Gene_expression_visualization.svg",width=16,height=8)
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"), pt.size.factor = 1, alpha = c(0.1,1))
dev.off()

##Dimensionality reduction, clustering, and visualization
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
saveRDS(brain, file = "brain_PCA.RDS")

#find dims
#JackStraw cannot use in STC data.
#but we can use ElbowPlot to find the best dims
svg("brain_ElbowPlot.svg", width = 10, height = 6)
ElbowPlot(brain, ndims = 30)
dev.off()

brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
svg("data_clustering.svg",width=16,height=8)
plot(p1 + p2)
dev.off()

svg("data_clustering_facetwrap.svg",width=16,height=8)
SpatialDimPlot(brain,
               cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,5, 8)),
               facet.highlight = TRUE, ncol = 3)
dev.off()

##Interactive plotting
#SpatialDimPlot(brain, interactive = TRUE)
#SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
#LinkedDimPlot(brain)

##Identification of Spatially Variable Features
#way1
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
svg("Identification.svg",width=18,height=9)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
dev.off()

#way2
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
svg("Identification.svg",width=16,height=8)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
dev.off()

saveRDS(brain ,"./brain.rds")

brain <- readRDS("brain.rds")

##Subset out anatomical regions
Idents(brain)
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
#svg("choose_subset.svg",width=1660,height=800)
#SpatialDimPlot(cortex,
#               cells.highlight = WhichCells(cortex, expression = imagerow > 400 | imagecol < 150))
#dev.off()

#cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
#cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
#cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

#p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
#p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
#p1 + p2

##Integration with single-cell data
allen_reference <- readRDS("data_ref_umap_annotated.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50)

# After subsetting, we renormalize cortex
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
svg("allen_reference_annotation.svg",width=16,height=8)
DimPlot(allen_reference, group.by = "ident", label = TRUE)
dev.off()

anchors <- FindTransferAnchors(reference = allen_reference, query = brain, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference@active.ident, prediction.assay = TRUE,
                                  weight.reduction = brain[["pca"]], dims = 1:50)
brain[["predictions"]] <- predictions.assay

DefaultAssay(brain) <- "predictions"

svg("cluster_location.svg",width=16,height=8)
SpatialFeaturePlot(brain,
                   features = c("Oligodendrocytes", "Neurons", "Microglia", "Astrocytes"),
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)
dev.off()

brain <- FindSpatiallyVariableFeatures(brain, assay = "predictions", selection.method = "markvariogram",
                                       features = rownames(brain), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(brain), 4)

plot1 <- SpatialPlot(object = brain, features = top.clusters, ncol = 2)
plot2 <- SpatialFeaturePlot(brain,
                            features = c("Oligodendrocytes", "Nsc", "Neurons", "Microglia", "Astrocytes", "Neuroepithelials", "T cells"),
                            pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

svg("spatial_celltype.svg",width=16,height=8)
SpatialPlot(object = brain, features = top.clusters, ncol = 2)
dev.off()

svg("spatial_celltype_all.svg",width=16,height=33)
plot(plot2)
dev.off()

saveRDS(brain, "brain_transfer_annotated_final.rds")
