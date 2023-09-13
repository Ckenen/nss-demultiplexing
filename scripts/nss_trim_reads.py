#!/usr/bin/env python
import sys
from collections import defaultdict
import edlib
from pygz import PigzFile
from pyBioInfo.IO.File import FastqFile


DEBUG = False


def load_fastq(path):
    with FastqFile(path) as f:
        for i, read in enumerate(f):
            if DEBUG and i >= 10000:
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


def cal_perc(n1, n2):
    if n2 == 0:
        return 0.0
    else:
        return n1 * 100 / n2


def write_fastq(fw, name, sequence, quality):
    fw.write("@%s\n" % name)
    fw.write("%s\n" % sequence)
    fw.write("+\n")
    fw.write("%s\n" % quality)

MIN_LENGTH = 200
MAX_ED = 8
MAX_CHIMIRIC = 6

def main():
    infile, outfile = sys.argv[1:]
    print("Infile: %s" % infile)
    print("Outfile: %s" % outfile)
    
    linker1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    linker2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
    print("Linker 1:", linker1)
    print("Linker 2:", linker2)
    print("Minimum length:", MIN_LENGTH)
    print("Threshold edit distance:", MAX_ED)
    print("Threshold for chimeric:", MAX_CHIMIRIC)
    
    n_total = 0
    n_too_short = 0
    n_no_linker = 0
    n_chimeric = 0
    n_pass = 0  
    
    counter_length = defaultdict(int)
    counter_linker = defaultdict(int)
    counter_chimeric = defaultdict(int)
    
    fw = PigzFile(outfile, "wt")
    
    for read in load_fastq(infile):
        n_total += 1
        
        seq = read.sequence
        qua = read.quality
        
        counter_length[len(seq)] += 1
        
        if len(seq) >= MIN_LENGTH:
            w = 40
            s1 = seq[:w]
            s2 = seq[-w:]
            offset = len(seq) - w
            x1, y1, ed1 = edlib_align(linker1, s1)
            x2, y2, ed2 = edlib_align(linker2, s2)
            x2, y2 = x2 + offset, y2 + offset
            
            counter_linker[(ed1, ed2)] += 1
            
            if ed1 <= MAX_ED and ed2 <= MAX_ED:
                # trim
                seq = seq[y1:x2]
                qua = qua[y1:x2]
                
                # chimeric
                x1, y1, ed1 = edlib_align(linker1, seq)
                x2, y2, ed2 = edlib_align(linker2, seq)
                counter_chimeric[(ed1, ed2)] += 1
                if ed1 <= MAX_CHIMIRIC or ed2 <= MAX_CHIMIRIC:
                    n_chimeric += 1
                    continue
                else:
                    write_fastq(fw, read.name, seq, qua)
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
        
    print("Total: %d (%.2f%%)" % (n_total, cal_perc(n_total, n_total)))
    print("Too short: %d (%.2f%%)" % (n_too_short, cal_perc(n_too_short, n_total)))
    print("No linker: %d (%.2f%%)" % (n_no_linker, cal_perc(n_no_linker, n_total)))
    print("Is chimeric: %d (%.2f%%)" % (n_chimeric, cal_perc(n_chimeric, n_total)))
    print("Pass: %d (%.2f%%)" % (n_pass, cal_perc(n_pass, n_total)))
    
    with open(outfile + ".length.txt", "w+") as fw:
        total = sum(counter_length.values())
        for length, count in sorted(counter_length.items()):
            fw.write("\t".join(map(str, [length, count, cal_perc(count, total)])) + "\n")
            
    with open(outfile + ".linker.txt", "w+") as fw:
        total = sum(counter_linker.values())
        for (ed1, ed2), count in sorted(counter_linker.items()):
            fw.write("\t".join(map(str, [ed1, ed2, count, cal_perc(count, total)])) + "\n")
            
    with open(outfile + ".chimeric.txt", "w+") as fw:
        total = sum(counter_chimeric.values())
        for (ed1, ed2), count in sorted(counter_chimeric.items()):
            fw.write("\t".join(map(str, [ed1, ed2, count, cal_perc(count, total)])) + "\n")
        
    print("Finished!")
    

if __name__ == '__main__':
    main()
