#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import edlib
from pygz import PigzFile
from pyBioInfo.IO.File import FastqFile


DEBUG = False


def load_fastq(path):
    with FastqFile(path) as f:
        for i, read in enumerate(f):
            if DEBUG and i >= 100000:
                break
            yield read


def get_cell_name(path):
    cell = path.split("/")[-1]
    if cell.endswith(".gz"):
        cell = cell[:-3]
    if cell.endswith(".fastq"):
        cell = cell[:-6]
    if cell.endswith(".fq"):
        cell = cell[:-3]
    return cell


def edlib_align(que, ref):
    r = edlib.align(que, ref, task="locations", mode="HW")
    ed = r["editDistance"]
    x, y = r["locations"][0]
    y += 1
    return x, y, ed


def get_perc(n1, n2):
    if n2 == 0:
        return 0
    else:
        return n1 * 100 / n2


def write_fastq(fw, name, sequence, quality):
    fw.write("@%s\n" % name)
    fw.write("%s\n" % sequence)
    fw.write("+\n")
    fw.write("%s\n" % quality)

MIN_LENGTH = 200
MAX_LINKER_ED = 8
MAX_CHIMIRIC_LINKER_ED = 8

def main():
    infile, outdir = sys.argv[1:]
    print("Infile: %s" % infile)
    print("Outdir: %s" % outdir)
    
    linker1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    linker2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
    print("Linker 1:", linker1)
    print("Linker 2:", linker2)
    print("Minimum length:", MIN_LENGTH)
    print("Threshold edit distance:", MAX_LINKER_ED)
    print("Threshold for chimeric:", MAX_CHIMIRIC_LINKER_ED)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = os.path.join(outdir, "trimmed.fastq.gz")
    
    n_total = 0
    n_too_short = 0
    n_no_linker = 0
    n_chimeric = 0
    n_pass = 0  
    
    counter_length = defaultdict(int)
    counter_linker = defaultdict(int)
    counter_chimeric = defaultdict(int)

    counter_linker_edge = defaultdict(int) # max linker edit distance at edge
    counter_linker_inner = defaultdict(int) # min linker edit distance at inner
        
    fw = PigzFile(outfile, "wt")
    
    for read in load_fastq(infile):
        n_total += 1
        name, seq, qua = read.name, read.sequence, read.quality
        counter_length[len(seq)] += 1
        if len(seq) >= MIN_LENGTH:
            w = 40
            s1, s2 = seq[:w], seq[-w:]
            offset = len(seq) - w
            x1, y1, ed1 = edlib_align(linker1, s1)
            x2, y2, ed2 = edlib_align(linker2, s2)
            x2, y2 = x2 + offset, y2 + offset
            
            counter_linker[(ed1, ed2)] += 1
            counter_linker_edge[max(ed1, ed2)] += 1
            
            if max(ed1, ed2) <= MAX_LINKER_ED:
                # trim
                seq, qua = seq[y1:x2], qua[y1:x2]
                
                # chimeric
                x1, y1, ed1 = edlib_align(linker1, seq)
                x2, y2, ed2 = edlib_align(linker2, seq)
                counter_chimeric[(ed1, ed2)] += 1
                counter_linker_inner[min(ed1, ed2)] += 1
                if min(ed1, ed2) <= MAX_CHIMIRIC_LINKER_ED:
                    n_chimeric += 1
                    continue
                else:
                    write_fastq(fw, name, seq, qua)
                    n_pass += 1
            else:
                n_no_linker += 1
                continue
        else:
            n_too_short += 1
            continue
            
    cell_name = get_cell_name(infile)
    print("Cell: %s" % cell_name)
        
    fw.close()
        
    print("Total: %d (%.2f%%)" % (n_total, get_perc(n_total, n_total)))
    print("Too short: %d (%.2f%%)" % (n_too_short, get_perc(n_too_short, n_total)))
    print("No linker: %d (%.2f%%)" % (n_no_linker, get_perc(n_no_linker, n_total)))
    print("Is chimeric: %d (%.2f%%)" % (n_chimeric, get_perc(n_chimeric, n_total)))
    print("Pass: %d (%.2f%%)" % (n_pass, get_perc(n_pass, n_total)))
    
    with open(os.path.join(outdir, "stats.tsv"), "w+") as fw:
        fw.write("Total\tTooShort\tNoLiner\tIsChimeric\tPass\n")
        fw.write("\t".join(map(str, [n_total, n_too_short, n_no_linker, n_chimeric, n_pass])) + "\n")
    
    with open(os.path.join(outdir, "length.tsv"), "w+") as fw:
        fw.write("Length\tCount\tRatio\n")
        total = sum(counter_length.values())
        for length, count in sorted(counter_length.items()):
            fw.write("\t".join(map(str, [length, count, get_perc(count, total)])) + "\n")
            
    with open(os.path.join(outdir, "linker_at_edge.tsv"), "w+") as fw:
        fw.write("ED1\tED2\tCount\tRatio\n")
        total = sum(counter_linker.values())
        for (ed1, ed2), count in sorted(counter_linker.items()):
            fw.write("\t".join(map(str, [ed1, ed2, count, get_perc(count, total)])) + "\n")
            
    with open(os.path.join(outdir, "linker_at_inner.tsv"), "w+") as fw:
        fw.write("ED1\tED2\tCount\tRatio\n")
        total = sum(counter_chimeric.values())
        for (ed1, ed2), count in sorted(counter_chimeric.items()):
            fw.write("\t".join(map(str, [ed1, ed2, count, get_perc(count, total)])) + "\n")

    # Maximum edit distance at read edge
    with open(os.path.join(outdir, "linker_at_edge_max.tsv"), "w+") as fw:
        fw.write("MaxED\tCount\tRatio\tCumulative\n")
        r0 = 0
        t = sum(counter_linker_edge.values())
        for k, v in sorted(counter_linker_edge.items()):
            r = v / t
            r0 += r
            fw.write("\t".join(map(str, [k, v, r, r0])) + "\n")

    # Minimum edit distance at read inner
    with open(os.path.join(outdir, "linker_at_inner_min.tsv"), "w+") as fw:
        fw.write("MinED\tCount\tRatio\tCumulative\n")
        r0 = 0
        t = sum(counter_linker_inner.values())
        for k, v in sorted(counter_linker_inner.items()):
            r = v / t
            r0 += r
            fw.write("\t".join(map(str, [k, v, r, r0])) + "\n")


    # print("Finished!")
    

if __name__ == '__main__':
    main()
