# Count reads for features in each cell
# Author: Yuanhua Huang
# Date: 27-04-2020

import os
import sys
import gzip
import time
import pysam
import subprocess
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from .version import __version__
from .utils.gtf_utils import load_genes
from .utils.count_utils import feature_count

FID = None
PROCESSED = 0
TOTAL_GENE = 0
START_TIME = time.time()

def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME, FID
    if RV is not None:
        FID.writelines(RV)
    
    PROCESSED += 1
    bar_len = 20
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[cf-count] [%s] %.1f%% done in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()
    return RV

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("Indexed & sorted sam/bam/cram file, possibly from CellRanger."))
    parser.add_option("--outDri", "-o", dest="out_dir", default=None,
        help=("Full path for output direcotry [default: $samFile_dir]"))
    parser.add_option("--barcodeFile", "-b", dest="barcode_file", default=None,
        help=("A plain file listing all effective cell barcodes."))
    parser.add_option("--GFF", "-G", dest="GFF_file", default=None,
        help=("A GFF/GTF file containing features annotation"))
    
    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group1.add_option("--featureType", "-t", dest="feature_type", default="gene",
        help=("Feature type to culculate [default: %default]."))
    # group1.add_option("--IDattribute", "-i", dest="id_attribute", default="ID",
    #     help=("Id tag for feature in GFF file [default: %default]."))
    group1.add_option("--cellTAG", dest="cell_tag", default="CB", 
        help="Tag for cell barcodes, turn off with None [default: %default]")
    group1.add_option("--UMItag", dest="UMI_tag", default="None", 
        help="Tag for UMI: UR (common for 10x RNA). None means no UMI but read "
        "counts [default: %default]")
    
    group2 = OptionGroup(parser, "Read filtering")
    group2.add_option("--minLEN", type="int", dest="min_LEN", default=30, 
        help="Minimum length mapped to the feature [default: %default]")
    group2.add_option("--minMAPQ", type="int", dest="min_MAPQ", default=20, 
        help="Minimum MAPQ for read filtering [default: %default]")
    group2.add_option("--maxFLAG", type="int", dest="max_FLAG", default=255, 
        help="Maximum FLAG for read filtering [default: %default]")
    
    parser.add_option_group(group1)
    parser.add_option_group(group2)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to cf-count in cellSNP v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.sam_file is None or options.barcode_file is None:
        print("Error: need samFile and barcodeFile.")
        sys.exit(1)

    sam_file = options.sam_file
    barcodes = np.genfromtxt(options.barcode_file, dtype="str", delimiter="\t")
    barcodes = sorted(list(barcodes))
    cell_tag = options.cell_tag
    UMI_tag  = options.UMI_tag if options.UMI_tag.upper() != "NONE" else None
        
    if options.out_dir is None:
        out_dir = os.path.dirname(options.sam_file)
    elif os.path.abspath(options.out_dir) == "":
        out_dir = "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if os.path.isdir(os.path.abspath(out_dir)) == False:
        print("Error: No such directory for file\n -- %s" %out_dir)
        sys.exit(1)        
    
    nproc = options.nproc
    min_LEN = options.min_LEN
    min_MAPQ = options.min_MAPQ
    max_FLAG = options.max_FLAG
    
    # load features from GTF file
    genes = load_genes(options.GFF_file, tranTag="", exonTag="")
    gene_ids = [g.geneID for g in genes]
    if len(np.unique(gene_ids)) != len(gene_ids):
        gene_ids = ["gene%d" %(x + 1) for x in range(len(genes))]

    # save cell barcodes and genes
    fid = open(out_dir + "/features.tsv", "w")
    for _gene_id in gene_ids:
        fid.writelines(_gene_id + "\n")
    fid.close()
    print("[cf-count] Count reads for %d features in %d cells with %d cores." 
          %(len(gene_ids), len(barcodes), nproc))

    fid = open(out_dir + "/barcodes.tsv", "w")
    for _cell_id in barcodes:
        fid.writelines(_cell_id + "\n")
    fid.close()

    global FID, TOTAL_GENE
    TOTAL_GENE = len(gene_ids)
    FID = open(out_dir + "/matrix.mtx", "w")

    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes=nproc)
        for ii in range(len(genes)):
            result.append(pool.apply_async(feature_count, (sam_file, 
                barcodes, genes[ii], ii, cell_tag, UMI_tag, min_MAPQ, max_FLAG, min_LEN), 
                callback=show_progress))
        pool.close()
        pool.join()
    else:
        for ii in range(len(genes)):
            RV = feature_count(sam_file, barcodes, genes[ii], ii, cell_tag, 
                               UMI_tag, min_MAPQ, max_FLAG, min_LEN)
            show_progress(RV)

    FID.close()

    data = np.genfromtxt(out_dir + "/matrix.mtx")
    fid = open(out_dir + "/matrix.mtx", "w")
    fid.writelines(
        "%" + "%MatrixMarket matrix coordinate integer general\n" + "%\n")
    fid.writelines("%d\t%d\t%d\n" %(len(genes), len(barcodes), data.shape[0]))
    
    sort_idx = np.argsort(data[:, 0] * len(barcodes) + data[:, 1])
    for ii in sort_idx:
        fid.writelines("%d\t%d\t%d\n" %(data[ii, 0], data[ii, 1], data[ii, 2]))
    fid.close()
    
    print("")
    

if __name__ == "__main__":
    main()
