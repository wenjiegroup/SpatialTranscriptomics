library(devtools)
install_github("navinlabcode/CellTrek")

options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")

setwd("/data1/yangpin/jbj/cellTrek/")
data_dir = "/data1/yangpin/jbj/cellTrek/"

brain_st <- readRDS("brain_transfer_annotated_final.rds")
brain_sc <- readRDS("data_ref_umap_annotated.rds")

## Rename the cells/spots with syntactically valid names
brain_st <- RenameCells(brain_st, new.names=make.names(Cells(brain_st)))
brain_sc <- RenameCells(brain_sc, new.names=make.names(Cells(brain_sc)))

## Visualize the ST data
svg("ST_data.svg",width=16,height=8)
SpatialDimPlot(brain_st)
dev.off()

## Visualize the scRNA-seq data
svg("SC_data.svg",width=16,height=8)
DimPlot(brain_sc, label = T, label.size = 4.5)
dev.off()

#co-embedding the sc and st
brain_sc@meta.data["cell_type"] <- brain_sc@active.ident
brain_traint <- CellTrek::traint(st_data=brain_st, sc_data=brain_sc, sc_assay='RNA', cell_names='cell_type')

##check the co-embedding result to see if there is overlap between these two data modalities
svg("co-embed_ST_sc.svg",width=16,height=8)
DimPlot(brain_traint, group.by = "type") 
dev.off()

#chart single cells to their spatial locations use the non-linear interpolation
brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

#interactively visualize the CellTrek result using celltrek_vis
brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=sort(unique(brain_celltrek$cell_type)))

CellTrek::celltrek_vis(brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
                       brain_celltrek@images$tissue_hires_image.png@image, brain_celltrek@images$tissue_hires_image.png@scale.factors$lowres)

##Cell colocalization analysis
#summarize the colocalization patterns between different cell types using SColoc module
brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=unique(brain_celltrek$cell_type))

brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek, col_cell='cell_type', use_method='KL', eps=1e-50)

## We extract the minimum spanning tree (MST) result from the graph
brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- Fibro_cell[colnames(brain_sgraph_KL_mst_cons)]
## We then extract the metadata (including cell types and their frequencies)
brain_cell_class <- brain_celltrek@meta.data %>% dplyr::select(id=cell_type) %>% unique
brain_celltrek_count <- data.frame(freq = table(brain_celltrek$cell_type))
brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")

#visualize the colocalization result
CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons, meta_data=brain_cell_class_new)


##Spatial-weighted gene co-expression analysis within the cell type of interest
#further investigate the co-expression patterns within the cell type of interest using SCoexp module

#this need cluster with gene
brain_celltrek_Ast <- subset(brain_celltrek, subset=cell_type=='Astrocytes')
brain_celltrek_Ast@assays$RNA@scale.data <- matrix(NA, 1, 1)
brain_celltrek_Ast$cluster <- gsub('Neurons VISp ', '', brain_celltrek_Ast$cell_type)

svg("Astrocytes_genes_cluster.svg",width=16,height=8)
DimPlot(brain_celltrek_Ast, group.by = 'cell_type')
dev.off()

#We select top 2000 variable genes (exclude mitochondrial, ribosomal and high-zero genes)
brain_celltrek_Ast <- FindVariableFeatures(brain_celltrek_Ast)
vst_df <- brain_celltrek_Ast@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_Ast[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(brain_celltrek_Ast), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_Ast), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]

#use scoexp to do the spatial-weighted gene co-expression analysis
brain_celltrek_Ast_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_Ast, assay='RNA', approach='cc', gene_select = feature_temp, sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)

#visualize the co-expression modules using heatmap
brain_celltrek_Ast_k <- rbind(data.frame(gene=c(brain_celltrek_Ast_scoexp_res_cc$gs[[1]]), G='K1'), 
                              data.frame(gene=c(brain_celltrek_Ast_scoexp_res_cc$gs[[2]]), G='K2')) %>% 
  magrittr::set_rownames(.$gene) %>% dplyr::select(-1)

svg("coexpression_heatmap.svg",width=16,height=8)
pheatmap::pheatmap(brain_celltrek_Ast_scoexp_res_cc$wcor[rownames(brain_celltrek_Ast_k), rownames(brain_celltrek_Ast_k)], 
                   clustering_method='ward.D2', annotation_row=brain_celltrek_Ast_k, show_rownames=F, show_colnames=F, 
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='Neurons spatial co-expression')
dev.off()

brain_celltrek_Ast <- AddModuleScore(brain_celltrek_Ast, features=brain_celltrek_Ast_scoexp_res_cc$gs, name='CC_', nbin=10, ctrl=50, seed=42)
## First we look into the coexpression module based on the scRNA-seq embedding

svg("coexpression_scRNA-seq_embedding.svg",width=16,height=10)
FeaturePlot(brain_celltrek_Ast, grep('CC_', colnames(brain_celltrek_Ast@meta.data), value=T))
dev.off()

svg("coexpression_scores_spatial.svg",width=16,height=8)
SpatialFeaturePlot(brain_celltrek_Ast, grep('CC_', colnames(brain_celltrek_Ast@meta.data), value=T))
dev.off()
