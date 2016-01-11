"""

"""
import os
import os.path as op

from argparse import ArgumentParser
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import gzip

from collections import Counter
import pandas as pd
# from ggplot import *

import vcf


def _sv_dist(fn_in):
    count = Counter()
    vcf_reader = vcf.Reader(open(fn_in, 'r'))
    samples = vcf_reader.samples
    remove_sr = remove_common = bdn = total = 0
    for record in vcf_reader:
        if record.genotype(samples[0])['GT'] == '0/0':
            continue # continue if genotype 0/0 in tumor
        if record.genotype(samples[1])['GT'] == '0/0':
            # count if genotype in control is 0/0 and
            # reads in tumor  > 5
            if record.genotype(samples[0])['SR'] > 5:
                count[(record.CHROM, record.var_subtype)] += 1
                if record.INFO["SVTYPE"] == "BND":
                    bdn += 1
                total += 1
            else:
                remove_sr += 1
        else:
            remove_common += 1
    print "Removed common %s, removed SR %s, total %s, BND %s" % (remove_sr, remove_common, total, bdn)
    return count

def _parse_count(count, sample):
    tab = []
    row = 0
    for k in count:
        row +=1
        tab.append([sample, k[0], k[1], count[k]])
    return  tab

def simple_report(args):
    summary = []
    for sample in args.files:
        name = os.path.basename(sample).replace("_Tumor-lumpy-pair.vcf.gz", "")
        dt = _sv_dist(sample)
        dt = pd.DataFrame(_parse_count(dt, name))
        dt.columns = ['sample', 'chrom', 'sv', 'counts']
        out_file = op.join(args.out, name + "_lumpy.tsv")
        dt.to_csv(out_file, sep='\t', index=False)
        summary.append(dt)

    out_file = op.join(args.out, "lumpy.tsv")
    pd.concat(summary).to_csv(out_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = ArgumentParser(description="clean SV VCF files.")
    parser.add_argument("--out",required=1)
    parser.add_argument("files", nargs="*", help="bcbio final yaml.")
    args = parser.parse_args()

    simple_report(args)
