# Utilility functions for count reads in each feature
# Author: Yuanhua Huang
# Date: 27/04/2020

## Note, this is a very basic read counting for features
 
import sys
import pysam
import numpy as np
from .base_utils import id_mapping, unique_list
from .pileup_utils import check_pysam_chrom
from ..version import __version__


global CACHE_CHROM
global CACHE_SAMFILE
CACHE_CHROM = None
CACHE_SAMFILE = None


def fetch_reads(sam_file, gene, cell_tag="CR", UMI_tag="UR", min_MAPQ=20, 
                max_FLAG=255, min_LEN=30):
    """ Fetch all reads mapped to a given gene.
    Filtering is also applied, including cell and UMI tags and read mapping 
    quality.
    """
    if sam_file is None or gene is None:
        if sam_file is None:
            print("Warning: samFile is None")
        if gene is None:
            print("Warning: gene is None")
        return np.array([]), np.array([])

    samFile, _chrom = check_pysam_chrom(sam_file, gene.chrom)

    UMIs_list, cell_list = [], []
    for _read in samFile.fetch(_chrom, gene.start - 1, gene.stop):
        ## filtering reads
        overhang = sum((_read.positions >= gene.start) * 
                       (_read.positions <= gene.stop))

        if _read.mapq < min_MAPQ or _read.flag > max_FLAG or overhang < min_LEN: 
            continue
        if cell_tag is not None and _read.has_tag(cell_tag) == False: 
            continue
        if UMI_tag is not None and _read.has_tag(UMI_tag) == False: 
            continue

        if UMI_tag is not None:
            UMIs_list.append(_read.get_tag(UMI_tag))
        if cell_tag is not None:
            cell_list.append(_read.get_tag(cell_tag))

    if len(cell_list) > 0 and len(cell_list) == len(UMIs_list):
        UMI_cell = [UMIs_list[x] + cell_list[x] for x in range(len(UMIs_list))]
        UMI_cell, idx, cnt = unique_list(UMI_cell)
        cell_list = [cell_list[x] for x in idx]
    
    cell_list_uniq, idx, read_count = unique_list(cell_list)

    return cell_list_uniq, read_count


def feature_count(sam_file, barcodes, gene, g_index, cell_tag, UMI_tag, 
    min_MAPQ, max_FLAG, min_LEN):
    """Fetch read count for a given feature.
    """
    cell_list_uniq, read_count = fetch_reads(sam_file, gene, cell_tag, UMI_tag, 
        min_MAPQ, max_FLAG, min_LEN)

    if len(cell_list_uniq) > 0:
        match_idx = id_mapping(cell_list_uniq, barcodes, uniq_ref_only=True, 
            IDs2_sorted=True)
        match_idx = match_idx.astype(float)

        idx1 = np.where(match_idx == match_idx)[0] #remove None
        idx2 = match_idx[idx1].astype(int)

        out_list = ["%d\t%d\t%d" %(g_index, idx1[x], read_count[idx2[x]])
            for x in range(len(idx1))]
        return "\t".join(out_list)
    else:
        return None
