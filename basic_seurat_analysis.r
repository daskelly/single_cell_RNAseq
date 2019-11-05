library(Seurat)    # Seurat v3
library(magrittr)
library(assertthat)
# approx current as of Nov 4 2019

umi <- Read10X("path/to/10X_library/outs/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(umi)
# obj$batch <- "name"
mt_threshold <- 25
obj %<>% PercentageFeatureSet(pattern="^mt-", col.name="percent.mt") %>%
	subset(subset=percent.mt < mt_threshold) %>%
    NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(vars.to.regress=c("percent.mt", "nFeature_RNA")) %>%
    RunPCA(verbose=FALSE, npcs=80)
ElbowPlot(obj, ndims=50)
num_pc <- NA
obj %<>% FindNeighbors(reduction='pca', dims=1:num_pc, verbose=FALSE) %>%
	FindClusters(verbose=FALSE) %>%
    RunUMAP(reduction='pca', dims=1:num_pc, verbose=FALSE) %>% 
