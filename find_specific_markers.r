# Based on https://gist.github.com/daskelly/1ad3ab8a840ec5b297938a537dec5531


find_specific_markers <- function(obj, celltypes, n.cores,
	celltypes_to_find_markers=celltypes) {
    require(assertthat)
    assert_that(tolower(class(obj)) == "seurat")
    require(magrittr)
    require(Seurat)
    specific_markers <- NULL
    # Find specific markers for each cell type
    
    for (celltype1 in celltypes) {
        other_types <- setdiff(celltypes, celltype1)
        celltype_markers <- NULL
        for (celltype2 in other_types) {
            markers <- FindMarkers(obj, ident.1=celltype1, ident.2=celltype2,
                only.pos=TRUE, min.diff.pct=0.05, test.use="roc",
                max.cells.per.ident=1000, verbose=FALSE) %>%
                rownames_to_column("gene") %>%
                mutate(ident.2=celltype2)
            celltype_markers <- rbind(celltype_markers, markers)
        }
        celltype_markers$ident.1 <- celltype1
        specific_markers <- rbind(specific_markers, celltype_markers)
    } 
 
    # In specific_markers data.frame, ident.1 is which cell pop the marker is tagging
    # and ident.2 is which other cell pop we are contrasting with
    #
    # We want markers that are considerably higher in ident.1 than in any other pop.
    # As a basic filter we require any marker for ident.1 to be expressed in <50% of 
    # cells of each of the other cell types.
    #
    # Could filter this data using group_by(cluster) %>% top_n(10)
    npop <- length(celltypes)
    specific_markers %>% mutate(pct_diff=pct.1 - pct.2) %>%
        group_by(ident.1, gene) %>%
        filter(n() == npop - 1, all(pct.2 < 0.6)) %>%
        summarize(median_AUC=median(myAUC), pct.1=mean(pct.1), pct.2=mean(pct.2)) %>%
        mutate(cluster=ident.1) %>%
        filter(median_AUC > 0.65) %>%
        group_by(cluster) %>%
        arrange(desc(median_AUC))
}
