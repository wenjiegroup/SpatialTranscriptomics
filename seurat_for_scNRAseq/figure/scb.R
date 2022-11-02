library(dplyr)
library(Seurat)
library(SeuratObject)
library(SingleR)
library(celldex)

library(tidyverse)
library(patchwork)
library(ggplot2)

setwd("/data/scb/")
data_dir = "/data/scb/"

data_ref.data = Read10X_h5(filename = "filtered_feature_bc_matrix.h5")
data_ref = CreateSeuratObject(counts =data_ref.data, project="data_ref", min.cells =5, min.features = 200)
head(data_ref, 20)

###QC-线粒体
data_ref[["percent.mt"]] <- PercentageFeatureSet(data_ref, pattern = "^mt-")

data_ref <- subset(data_ref, subset = nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)

svg( paste0("data_ref","_VlnPlot2.svg"), width = 10, height = 6 )
plot(VlnPlot(data_ref, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
             ncol = 4,pt.size = 0))
dev.off()

##数据标准化
data_ref <- NormalizeData(data_ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
#鉴定细胞间表达量高变的基因
data_ref <- FindVariableFeatures(data_ref, selection.method = "vst", nfeatures = 2000, verbose = F)
#缩放数据(线性缩放：每个基因线性),PCA的输入数据要服从正态分布
data_ref <- ScaleData(data_ref, features = rownames(data_ref))
#PCA
data_ref <- RunPCA(data_ref, features = VariableFeatures(object = data_ref),npcs = 50)

# Examine and visualize PCA results a few different ways
print(data_ref[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(data_ref, dims = 1:2, reduction = "pca")
#DimPlot(data_ref, reduction = "pca")
#DimHeatmap(data_ref, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(data_ref, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
data_ref <- JackStraw(data_ref, dims = 50, num.replicate = 100)
data_ref <- ScoreJackStraw(data_ref, dims = 1:50)

svg( paste0("data_ref","_JackSPlot.svg"), width = 10, height = 6 )
JackStrawPlot(data_ref, dims = 1:50)
dev.off()

svg( paste0("data_ref","_ElbowPlot.svg"), width = 10, height = 6 )
ElbowPlot(data_ref,ndims = 50)
dev.off()

#cluster
data_ref <- FindNeighbors(data_ref, dims = 1:50)
data_ref <- FindClusters(data_ref, resolution =0.5 )
#umap
data_ref <- RunUMAP(data_ref, dims = 1:50 , label = T, reduction = "pca")

plot.umap1 <- DimPlot(data_ref, reduction = "umap", group.by = "orig.ident", pt.size=0.02, label = TRUE,repel = TRUE)
plot.umap2 <- DimPlot(data_ref, reduction = "umap", group.by = "ident",   pt.size=0.02, label = TRUE,repel = TRUE)

svg(paste0("data_ref","_umap.svg"),width=16,height=8)
plot(plot.umap2+plot.umap1)
dev.off()

saveRDS(data_ref ,"data_ref_umap_fine.rds")

data_ref.marker = FindAllMarkers(data_ref, min.pct = 0.25)

######################################### SingleR ###########################################
library(SingleR)
library(celldex)

immgen <- ImmGenData()
MouseRNAseq = MouseRNAseqData()
pred2 <- SingleR(test = data_ref@assays[["RNA"]]@data, ref = immgen, labels = immgen$label.fine, clusters = data_ref@active.ident, assay.type.test=1)
pred2@listData[["labels"]]

pred4 <- SingleR(test = data_ref@assays[["RNA"]]@data, ref = MouseRNAseq, labels = MouseRNAseq$label.fine, clusters = data_ref@active.ident, assay.type.test=1)
pred4@listData[["labels"]]

pred=data.frame(cluster=c(0:20),pred2=pred2@listData[["labels"]],
                pred4=pred4@listData[["labels"]])
pred

data_ref <- readRDS("data_ref_umap_fine.rds")

data_ref <- RenameIdents(data_ref, 
                         `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Neurons",`3` = "Neurons", `4` = "Astrocytes", 
                         `5` = "Neurons",`6` = "Microglia", `7` = "Astrocytes", `8` = "Neurons", `9` = "Neurons" ,
                         `10` = "Neuroepithelials", `11` = "Neuroepithelials", `12` = "Oligodendrocytes", `13` = "Neurons", `14` = "Neuroepithelials",
                         `15` = "Nsc", `16` = "Neuroepithelials", `17` = "Neurons", `18` = "Neuroepithelials", `19` = "Neurons",
                         `20` = "T cells")

plot.umap <- DimPlot(data_ref, reduction = "umap", group.by = "ident", pt.size=0.02, label = TRUE, repel = TRUE)


svg( "annotated_umap.svg", width = 12, height = 6 )
plot(plot.umap)
dev.off()

saveRDS(data_ref ,"data_ref_umap_annotated.rds")
