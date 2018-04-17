argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) != 2) {
    cat("USAGE --> knn_out_to_sparseMat.r infile.tsv outfile.mtx\n")
    quit(status=1)
}

library(Matrix)
library(tidyverse)
infile <- argv[1]
outfile <- argv[2]
x <- data.table::fread(infile) %>%
    as.data.frame() %>% remove_rownames() %>%
    column_to_rownames("V1") %>% as.matrix() %>%
    Matrix(sparse=TRUE)

writeMM(x, file=outfile)
