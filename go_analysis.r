library(tidyverse)
source('~/gbdb/go_analysis_find_children.r')
# find_children() function

gene_assoc <- data.table::fread("~/gbdb/gene_association_2019_10_22.mgi", skip=24,
	col.names=c("db", "db_obj_id", "db_obj_symbol", "qualifier", "go_id", "db_ref",
		"evidence_code", "with_or_from", "aspect", "db_obj_name", "db_obj_synonym",
		"db_obj_type", "taxon", "date", "assigned_by", "annot_ext", "gene_prod_form_id")) %>%
	as_tibble()

get_genes <- function(go_term) {
    kids <- get_children(go_term, name=TRUE)
    rows <- filter(gene_assoc, go_id %in% c(go_term, kids$go_id))
    unique(rows$db_obj_symbol) %>% sort()
}
# Note: the function call below will return slightly different results than 
# http://www.informatics.jax.org/go/term/GO:0009986
# Notably missing features include
# polymorphic pseudogenes (e.g. Clec7a)
# and unclassified other genome feature (e.g. Ly76)
# These are not excluded through anything under my control, rather they are
# apparently excluded from gene_association.mgi
genes <- get_genes("GO:0009986")  # cell surface GO term

