library(Seurat)
library(ShinyCell)

setwd('~/fastscratch/heart MI scRNAseq/7TimePointsSeurat-selected')
obj <- readRDS("AggregatePC24res0.5withNames.rds")
obj <- UpdateSeuratObject(obj)
sobj <- SplitObject(obj, split.by='stim')
footnote = paste0(
  'strong("Reference: "), "Forte E., et al. ",',
  'em("Cell Reports "), strong("30,"), "3149-3163 (2020) ",',
  'a("doi:10.1016/j.celrep.2020.02.008",',
  'href = "https://doi.org/10.1016/j.celrep.2020.02.008",',
  'target="_blank"), style = "font-size: 125%;"'
)
daycodes <- names(sobj)
dayfull <- gsub("d", "day ", daycodes)
for (i in 1:length(sobj)) {
  message(i)
  daycode <- names(sobj)[i]
  a <- sobj[[daycode]]
  todel <- c('orig.ident', 'S.Score', 'G2M.Score', 'Phase', 'res.0.5')
  a@meta.data <- a@meta.data[, setdiff(colnames(a@meta.data), todel)]
  conf <- createConfig(a)
  conf <- modMetaName(conf, meta.to.mod=c('nCount_RNA', 'nFeature_RNA', 'CellType_0.5', 'stim'),
                      new.name=c("No. UMIs", "No. detected genes", "Cell type", 'Day post-MI'))
  conf <- reorderMeta(conf, conf$ID[c(3, 2, 1, 4:length(conf$ID))])
  makeShinyFiles(a, conf, gex.assay = "RNA", gex.slot = "data",
                 gene.mapping = TRUE, shiny.prefix = daycode,
                 shiny.dir = "shinyAppMulti/",
                 default.gene1 = "Adgre1", default.gene2 = "Pdgfra",
                 default.multigene = c("Ms4a1", "Ncr1", "Cd3e", "S100a9", 
                    "Cdh5", "Acta2", "Adgre1", "Col1a1", "Pdgfra"),
                 default.dimred=c('tSNE_1', 'tSNE_2'))
  
  # makeShinyCodes(shiny.title="Forte et al. 2020 heart MI", 
  #                shiny.footnotes = footnote,
  #                shiny.prefix=daycode, shiny.dir="shinyAppMulti/")
}
makeShinyCodesMulti(
  shiny.title = "Forte et al. 2020 heart MI", shiny.footnotes = footnote,
  shiny.prefix = daycodes,
  shiny.headers = dayfull, 
  shiny.dir = "shinyAppMulti/") 
