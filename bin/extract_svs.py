#!/usr/bin/env python
"""Extract structural variations in regions of interest into portable feather format for later analysis
"""
import csv
import os
import glob
import subprocess
import itertools
import pandas as pd
import feather
import logging

import pysam

logging.basicConfig(
        format="%(asctime)s %(levelname)-5s %(name)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO)

log = logging.getLogger(__file__)

INPUT_DIR="/home/rvalls/dev/cancer-hiveplots/hyplot/tests/priv"

def _find_svcaller(in_file):
    """ Finds out the caller(s) (VEP, Snpeff, BASS...) used in a VCF.
        
        This would be easy if all SV callers agreed on things like
        ##SVCmd = ... or ##SVVersion in the headers.
    """ 
    callers = []
    
    log.debug("Identifying SV caller used in {}...".format(in_file))
    # Look at the VCF headers for hints...
    for rec in pysam.VariantFile(in_file).header.records:
        log.debug("Looking for SVcaller signature in the header line:\n {}\n".format(rec))
        if "##source=LUMPY" in str(rec):
            callers.append("LUMPY")
        elif "SVCLASS" in str(rec):
            # Used by ICGC SV
            callers.append("BASS")

    # ... try with an actual variant if not enough
    first_variant = pysam.VariantFile(in_file).next().id

    if "Manta" in first_variant:
        callers.append("MANTA")
    elif "SVCLASS" in first_variant:
        callers.append("BASS")
    #else:
    #    log.debug("Could not find SVcaller, the variant looks line like this: {}".format(first_variant))
    #    callers.append("UNKNOWN")

    log.debug("Found callers {}".format(callers))

    return callers

def extract_svs(in_file, depth, chroms):
    """Create CSV file of structural variants of interest, for Circos plots.
    """
    allowed_chroms = set([str(x) for x in range(1, 23)])
    # print(allowed_chroms + 'X')
    # allowed_chroms = set(allowed_chroms.append('X'))

    df = pd.DataFrame(columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "file", "caller", "svtype"])

    log.debug("Building dataframe from VCF file...")

    callers = _find_svcaller(in_file)
    for idx, caller in enumerate(callers):
        for p1, p2, svtype in parse_svs(in_file, depth):
            if len(chroms) == 0 or (p1[0] in chroms or p2[0] in chroms):
                if p1[0] in allowed_chroms and p2[0] in allowed_chroms:
                    row = pd.Series({"chrom1": p1[0], "start1": p1[1], "end1": p1[2], 
                                     "chrom2": p2[0], "start2": p2[1], "end2": p2[2], 
                                     "file": in_file, "caller": caller, "svtype": svtype})
                    df = df.append(row, ignore_index=True)

    try:
        log.info("Exporting to interoperable feather file {}.feather".format(in_file))
        feather.write_dataframe(df, "{}.feather".format(in_file))
    except feather.ext.FeatherError:
        log.error("Failed to serialize feather object (most likely empty source dataframe)")

def parse_svs(in_file, depth):
    bnds = {}
    min_su_bnd = depth
    min_su_other = depth
    try:
        for rec in pysam.VariantFile(in_file):
            if passes(rec):
                if rec.info["SVTYPE"] and rec.info["SVTYPE"] == "BND":
                    if rec.info.get("SU") is None or rec.info["SU"] > min_su_bnd:
                        p = [rec.chrom.replace("chr", ""), rec.start, rec.start]
                        if rec.id in bnds:
                            p_o = bnds[rec.id]
                            # ICGC samples encode SVTYPE (all appear to be BND)
                            # into SVCLASS in INFO (insertion, deletion)
                            if rec.info["SVCLASS"] is not None:
                                yield p_o, p, rec.info["SVCLASS"]
                            else:
                                yield p_o, p, rec.info["SVTYPE"]
                            del bnds[rec.id]
                        else:
                            mate_id = rec.info.get("MATEID")
                            bnds[mate_id] = p
                else:
                    if rec.info.get("SU") is None or rec.info["SU"] > min_su_other:
                        p1 = [rec.chrom.replace("chr", ""), rec.start, rec.start]
                        p2 = [rec.chrom.replace("chr", ""), rec.info["END"], rec.info["END"]]
                        yield p1, p2, rec.info["SVTYPE"]
    except ValueError:
        print("Skipping possibly empty file: {sample}".format(sample=in_file))

def passes(rec):
    try:
        if rec.filter.keys() == 0 or rec.filter.keys()[0] == "PASS":
            if not rec.samples[0].get("GT"):
                return True
            elif list(set(rec.samples[0].get("GT"))) != ["N"]:
                return True
    except IndexError:
        log.warn("Record ID {} does not seem to have PASS and/or GT information".format(rec.id))
    
    return True


def main(in_dir):

    chroms = [str(x) for x in range(1, 23)]
    #chroms = chroms.append(['X', 'Y'])
    depths = [10, 25]
    approaches = ["full", "subset"]

    for fname in in_dir:
        for depth in depths:
            for approach in approaches:
                if approach == "full":
                    cur_chroms = []
                else:
                    cur_chroms = chroms 
                extract_svs(fname, depth, cur_chroms)

if __name__ == "__main__":
    in_dir = glob.glob("{}/*.vcf".format(INPUT_DIR))

    main(in_dir)
