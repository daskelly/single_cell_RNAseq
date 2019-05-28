# Given a set of genes in mouse,
# return a logical TRUE/FALSE vector
# indicating whether each gene belongs to
# a gene ontology category of interest.

# IMPORTANT NOTE: THIS FUNCTION DOES NOT ALWAYS 
# WORK CORRECTLY!!!
# For example, seeing whether Adgre1 (F4/80) falls
# within the cell surface GO category (see below)
# should return TRUE but instead returns FALSE.
# This appears to be due to a mistake in the 
# biomaRt database and not something I can control.

go_match <- function(genes, go_id_string) {
    require(tidyverse)
    require(biomaRt)
    if (! exists("ensembl")) {
        ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }
    gene.data <- getBM(attributes=c("mgi_symbol", "ensembl_gene_id", "go_id"),
        filters="go", values=go_id_string, mart=ensembl) %>%
        filter(go_id == go_id_string)
    result <- tibble(mgi_symbol=genes)
    result[[go_id_string]] <- genes %in% gene.data$mgi_symbol
    result
}
    
# Example:
# GO:0009986 is cell surface (1101 genes as of 5/28/2019)
# see http://www.informatics.jax.org/vocab/gene_ontology/GO:0009986
genes <- c("Adgre1", "Cd8a", "Csf1", "Csf1r", "Pdgfra", "Ptprc", "Tcf21")
go_match(genes, "GO:0009986")
