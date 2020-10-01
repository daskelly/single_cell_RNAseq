# Starting with a velocyto loom file, write a text file with
# the count of spliced and unspliced reads for each cell.
#
# Script is based on code from
# https://cellgeni.readthedocs.io/en/0.0.1/visualisations.html
#
# Run using
# singularity run /projects/skelld/sif/scanpy_Sept2020.sif 
import sys, argparse
import loompy
import scanpy as sc
import pandas
import numpy
import scipy

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("loomfile", metavar="file.loom")
    parser.add_argument("outfile", metavar="file.tsv")
    args = parser.parse_args()
    assert args.loomfile.endswith('.loom')
    assert args.outfile.endswith('.tsv')
    
    adata = sc.read_loom(args.loomfile)
    vals = dict()
    # sanity checks:
    assert adata.obs.shape[0] == adata.X.shape[0] # n cells
    assert adata.var.shape[0] == adata.X.shape[1] # n features/genes
    cells = adata.obs.index.values
    genes = adata.var.index.values
    outfile = open(args.outfile, 'w')
    outfile.write('cell_id\tX\tmatrix\tambiguous\tspliced\tunspliced\n')
    for i in range(len(cells)):
        outfile.write('{}\t{}'.format(cells[i], adata.X[i, :].sum()))
        for key in ['matrix', 'ambiguous', 'spliced', 'unspliced']:
            outfile.write('\t{}'.format(adata.layers[key][i, :].sum()))
        outfile.write('\n')
    outfile.close()

if __name__ == "__main__":
    sys.exit(main())
