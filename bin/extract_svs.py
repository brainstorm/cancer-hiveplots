#!/usr/bin/env python
"""Extract structural variations in regions of interest into tab delimited files for circos plots.
"""
import csv
import os
import glob
import subprocess

import pysam

def main():
    samples = {os.path.basename(k): ["TUMOR", "NORMAL"] for k in glob.glob("inputs/*.vcf.gz")}
    chroms = [str(x) for x in range(1, 23)]
    callers = ["BASS"]
    depths = [10, 25]
    approaches = ["full", "subset"]

    for fname, samples in samples.items():
        for depth in depths:
            for approach in approaches:
                if approach == "full":
                    cur_chroms = []
                else:
                    cur_chroms = chroms
                extract_svs(fname, samples, depth, cur_chroms, callers)

def extract_svs(project, samples, depth, chroms, callers):
    """Create CSV file of structural variants of interest, for Circos plots.
    """
    allowed_chroms = set([str(x) for x in range(1, 23)])
    out_dir = "sv_plots"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = os.path.join(out_dir, "%s-svs-%s-dp%s.csv" % (project, ("subset" if chroms else "full"), depth))
    region_bed = "%s-regions.bed" % os.path.splitext(out_file)[0]

    with open(out_file, "w") as out_handle:
        with open(region_bed, "w") as region_out:
            writer = csv.writer(out_handle)
            writer.writerow(["chrom1", "start1", "end1", "chrom2", "start2", "end2", "file", "sample", "caller", "svtype"])
            for sample in samples:
                for ci, caller in enumerate(callers):
                    sv_file = os.path.join("inputs", "%s" % project)
                    for p1, p2, svtype in parse_svs(sv_file, sample, depth):
                        if len(chroms) == 0 or (p1[0] in chroms or p2[0] in chroms):
                            if p1[0] in allowed_chroms and p2[0] in allowed_chroms:
                                writer.writerow(p1 + p2 + [project, sample, caller, svtype])
                                if ci == 0:
                                    for (chrom, start, end) in [p1, p2]:
                                        region_out.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, start - 1, end, project, sample))

def filter_shared(in_file, region_bed):
    cmd = "sort -k1,1 -k2,2n {region_bed} | bedtools merge -i - -d 10000 -c 4 -o count_distinct"
    overlaps = subprocess.check_output(cmd.format(**locals()), shell=True)
    shared = 0
    noshared = 0
    for line in overlaps.split("\n"):
        if line.strip():
            if int(line.strip().split()[-1]) > 1:
                shared += 1
            else:
                noshared += 1
    print shared, noshared
    return []

def parse_svs(in_file, sample, depth):
    bnds = {}
    min_su_bnd = depth
    min_su_other = depth
    for rec in pysam.VariantFile(in_file):
        if passes(rec):
            if rec.info["SVTYPE"] == "BND":
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

def passes(rec):
    #if rec.filter.keys() == 0 or rec.filter.keys()[0] == "PASS":
    #    if not rec.samples[0].get("GT"):
    #        return True
    #    elif list(set(rec.samples[0].get("GT"))) != ["N"]:
    #        return True
    return True

if __name__ == "__main__":
    main()
