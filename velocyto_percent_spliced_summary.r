# Provide a summary of per-cell statistics from an RNA velocity
# analysis as calculated using my accessory script
# /home/skelld/repos/single_cell_RNAseq/velocyto_percent_spliced_per_cell.py
#
# Call as:
# singularity run /projects/skelld/sif/tidyverse-4.0.0.sif velocyto_percent_spliced_summary.r file.stats.tsv
#
library(tidyverse)
library(assertthat)
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) != 1) quit(status=1)

dat <- read_tsv(argv[1]) %>% 
    mutate(pct_spliced=spliced/(spliced + unspliced))

# Output the middle 50% of the distribution:
#summary(dat$pct_spliced)
z <- quantile(dat$pct_spliced, c(.25, .75)) %>% unname() %>% round(2)
cat(argv[1], " ", paste0(z, collapse="--"), "\n")

