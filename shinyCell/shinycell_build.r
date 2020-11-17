library(Seurat)
library(ShinyCell)

setwd('~/fastscratch/heart MI scRNAseq/7TimePointsSeurat-selected')
obj <- readRDS("AggregatePC24res0.5withNames.rds")
obj <- UpdateSeuratObject(obj)
scConf <- createConfig(obj)
makeShinyApp(obj, scConf, gene.mapping = TRUE,
             shiny.title = "Forte et al 2020 heart MI")
