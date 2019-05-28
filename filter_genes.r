# Given a set of genes in mouse,
# return a logical TRUE/FALSE vector
# indicating whether each gene belongs to
# a gene ontology category of interest.

load_db <- function(dbfile) {
    require(tidyverse)
    data.table::fread(cmd=paste0("gunzip -c ", dbfile),
        sep="\t", skip=24,
        col.names=c("database_designation", "mgi_id", "mgi_symbol", 
        "not_designation", "go_id", "mgi_ref_id", "go_evidence_code", 
        "inferred_from", "ontology", "mgi_name", "mgi_synonyms", "type", 
        "taxon", "modification_date", "assigned_by", "annotation_extension", 
        "gene_product")) %>% as_tibble()
}

go_match <- function(genes, go_id_string) {
    require(tidyverse)
    if (! exists("go_db")) go_db <- load_db('gene_association.mgi.gz')
    result <- tibble(mgi_symbol=genes)
    this_db <- filter(go_db, go_id == go_id_string)
    result[[go_id_string]] <- genes %in% this_db$mgi_symbol
    result
}
    
# Example:
# GO:0009986 is cell surface (1101 genes as of 5/28/2019)
# see http://www.informatics.jax.org/vocab/gene_ontology/GO:0009986
genes <- c("Adgre1", "Cd8a", "Csf1", "Csf1r", "Pdgfra", "Ptprc", "Tcf21")
go_match(genes, "GO:0009986")
