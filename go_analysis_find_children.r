library(ontologyIndex)
go_basic <- get_ontology("~/gbdb/go-basic.obo", propagate_relationships=c("is_a", "part_of"))

get_children <- function(go_term, name=TRUE) {
    require(assertthat)
    require(magrittr)
    require(tidyverse)
    assert_that(go_term %in% names(go_basic$id))
    kids <- go_basic$children[[go_term]] %>% unique() %>% sort()
    kids2 <- sapply(kids, function(x) go_basic$children[[x]])
    kids2 <- c(names(kids2), unname(unlist(kids2))) %>% unique() %>% sort()
    while(! are_equal(kids, kids2)) {
        kids <- kids2
        kids2 <- sapply(kids, function(x) go_basic$children[[x]])
        kids2 <- c(names(kids2), unname(unlist(kids2))) %>% unique() %>% sort()
    }
    if (name) {
        name_df <- go_basic$name[kids] %>% enframe("go_id", "name")
    } else {
        kids
    }
}

# Example GO term for cell surface
# get_children("GO:0009986", name=TRUE)
# agrees with results at
# http://www.informatics.jax.org/vocab/gene_ontology/GO:0009986
# http://www.mousemine.org/mousemine/template.do?name=Term_Descendants&scope=all
