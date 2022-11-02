library(SingleCellExperiment)
library(ggplot2)
library(scHOT)
library(scater)
library(matrixStats)

setwd("/data/schot/")

data(MOB_subset)
sce_MOB_subset <- MOB_subset$sce_MOB_subset

sce_MOB_subset

scHOT_spatial <- scHOT_buildFromSCE(sce_MOB_subset,
                                    assayName = "logcounts",
                                    positionType = "spatial",
                                    positionColData = c("x", "y"))

pairs <- t(combn(rownames(sce_MOB_subset),2))
rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")
head(pairs)

set.seed(2020)
pairs <- pairs[sample(nrow(pairs), 20), ]
if (!"Arrb1_Mtor" %in% rownames(pairs)) {
  pairs <- rbind(pairs, "Arrb1_Mtor" = c("Arrb1", "Mtor"))
}
if (!"Dnm1l_Fam63b" %in% rownames(pairs)) {
  pairs <- rbind(pairs, "Dnm1l_Fam63b" = c("Dnm1l", "Fam63b"))
}

scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, pairs)
scHOT_spatial@testingScaffold

scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
                                       positionColData = c("x","y"),
                                       positionType = "spatial",
                                       nrow.out = NULL,
                                       span = 0.05)

dim(slot(scHOT_spatial, "weightMatrix"))

cellID = 75

svg("weighting_scheme.svg",width=16,height=8)
ggplot(as.data.frame(colData(scHOT_spatial)), aes(x = -x, y = y)) +
  geom_point(aes(colour = slot(scHOT_spatial, "weightMatrix")[cellID,],
                 size = slot(scHOT_spatial, "weightMatrix")[cellID,])) +
  scale_colour_gradient(low = "black", high = "purple") +
  scale_size_continuous(range = c(1,5)) +
  theme_classic() +
  guides(colour = guide_legend(title = "Spatial Weight"),
         size = guide_legend(title = "Spatial Weight")) +
  ggtitle(paste0("Central cell: ", cellID))
dev.off()

scHOT_spatial <- scHOT_calculateGlobalHigherOrderFunction(
  scHOT_spatial,
  higherOrderFunction = weightedSpearman,
  higherOrderFunctionType = "weighted")

slot(scHOT_spatial, "scHOT_output")

head(diag(cor(t(assay(scHOT_spatial, "expression")[pairs[,1],]),
              t(assay(scHOT_spatial, "expression")[pairs[,2],]),
              method = "spearman")))

scHOT_spatial <- scHOT_setPermutationScaffold(scHOT_spatial,
                                              numberPermutations = 50,
                                              numberScaffold = 10)

slot(scHOT_spatial, "scHOT_output")

scHOT_spatial <- scHOT_calculateHigherOrderTestStatistics(scHOT_spatial)

slot(scHOT_spatial, "scHOT_output")

system.time(scHOT_spatial <- scHOT_performPermutationTest(
  scHOT_spatial,
  verbose = TRUE,
  parallel = FALSE))

slot(scHOT_spatial, "scHOT_output")

scHOT_spatial <- scHOT_estimatePvalues(scHOT_spatial,
                                       nperm_estimate = 100,
                                       maxDist = 0.1)
slot(scHOT_spatial, "scHOT_output")

colData(scHOT_spatial)[, "-x"] <- -colData(scHOT_spatial)[, "x"]

svg("GG_interaction.svg",width=16,height=8)
plotHigherOrderSequence(scHOT_spatial, c("Dnm1l_Fam63b"),
                       positionColData = c("-x", "y"))
dev.off()

svg("weighting_GG.svg",width=16,height=8)
plotOrderedExpression(scHOT_spatial, c("Dnm1l", "Fam63b"),
                      positionColData = c("-x", "y"),
                      assayName = "expression")
dev.off()

scHOT_spatial <- scHOT_stripOutput(scHOT_spatial, force = TRUE)

scHOT_spatial

slot(scHOT_spatial, "scHOT_output")

scHOT_spatial <- scHOT(scHOT_spatial,
                       testingScaffold = pairs,
                       positionType = "spatial",
                       positionColData = c("x", "y"),
                       nrow.out = NULL,
                       higherOrderFunction = weightedSpearman,
                       higherOrderFunctionType = "weighted",
                       numberPermutations = 50,
                       numberScaffold = 10,
                       higherOrderSummaryFunction = sd,
                       parallel = FALSE,
                       verbose = FALSE,
                       span = 0.05)

slot(scHOT_spatial, "scHOT_output")

svg("weighting_interaction.svg",width=16,height=8)
plotHigherOrderSequence(scHOT_spatial, "Arrb1_Mtor",
                        positionColData = c("-x", "y"))
dev.off()
