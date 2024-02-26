#!/usr/bin/env python
import os
import optparse
from collections import defaultdict
from pygz import PigzFile
from pyBioInfo.IO.File import FastqFile


DEBUG = False


def load_fastq(path):
    with FastqFile(path) as f:
        for i, read in enumerate(f):
            if DEBUG and i >= 100000:
                break
            yield read                                
                
def load_matrix(path):
    with PigzFile(path) as f:
        for line in f:
            row = line.strip("\n").split("\t")
            for i in [1, 5, 6, 7, 11, 12, 13]:
                row[i] = int(row[i])
            yield row


def get_perc(n1, n2):
    if n2 == 0:
        return 0
    else:
        return n1 * 100 / n2
    

def demultiplexing(f_fastq, 
                   f_matrix, 
                   f_config, 
                   outdir, 
                   max_edit_distance=5, 
                   min_length=400, 
                   keep_unclassified=False, 
                   trim_outer=True):
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fqdir = os.path.join(outdir, "fastqs")
    if not os.path.exists(fqdir):
        os.mkdir(fqdir)
        
    # dual barcode pair
    cell_list = []
    cell_reads_counter = dict()
    ed_counter = defaultdict(int)
    bc2cell = dict()  # barcode to cell
    fws = dict()
    
    with open(f_config) as f:
        for line in f:
            cell, bc1, bc2 = line.strip("\n").split("\t")[:3]
            cell_list.append(cell)
            bc2cell[(bc1, bc2)] = cell  # (1st, 2nd)
            fw_f = open(os.path.join(fqdir, "%s_F.fastq" % cell), "w+")
            fw_r = open(os.path.join(fqdir, "%s_R.fastq" % cell), "w+")
            fws[cell] = [fw_f, fw_r]
            cell_reads_counter[cell] = [0, 0, 0]  # total forward reverse
    cell_reads_counter["unclassified"] = [0, 0, 0]
    unfw = None
    if keep_unclassified:
        unfw = open(os.path.join(fqdir, "unclassified.fastq") , "w+")
        
    n_total = 0
    n_too_short = 0
    n_no_direction = 0
    n_large_ed = 0
    n_invalid_cell = 0
    n_pass = 0
    for read, row in zip(load_fastq(f_fastq), load_matrix(f_matrix)):
        name, seq, qua = read.name.split()[0], read.sequence, read.quality
        assert name == row[0]
        bc1, direct1, loc1, x1, y1, ed1 = row[2:8]
        bc2, direct2, loc2, x2, y2, ed2 = row[8:14]   
        fw = None
        n_total += 1
        cell = None
        if len(seq) <= min_length:
            n_too_short += 1
        else:
            direction = None
            if direct1 == "F" and loc1 == "H" and direct2 == "R" and loc2 == "T":
                direction = "F"
            elif direct1 == "R" and loc1 == "T" and direct2 == "F" and loc2 == "H":
                direction = "R"   
            if direction is None:
                n_no_direction += 1
            else:
                ed = max(ed1, ed2)
                ed_counter[ed] += 1
                if ed > max_edit_distance:
                    n_large_ed += 1
                else:                    
                    cell = bc2cell.get((bc1, bc2))
                    if cell is None:
                        n_invalid_cell += 1
                    else:
                        if trim_outer:
                            x, y = min(y1, y2), max(x1, x2)
                            seq = seq[x:y]
                            qua = qua[x:y]
                        cell_reads_counter[cell][0] += 1
                        if direction == "F":
                            fw = fws[cell][0]
                            cell_reads_counter[cell][1] += 1
                        else:
                            fw = fws[cell][1]
                            cell_reads_counter[cell][2] += 1
                        n_pass += 1      
        if fw is None:
            fw = unfw
            cell_reads_counter["unclassified"][0] += 1 
        if fw:
            fw.write("@%s\n%s\n+\n%s\n" % (name, seq, qua))            
    for fw1, fw2 in fws.values():
        fw1.close()
        fw2.close()
    if unfw:
        unfw.close()
        
    if True:
        print("Total reads: %d (%.2f%%)" % (n_total, get_perc(n_total, n_total)))
        print("Too short: %d (%.2f%%)" % (n_too_short, get_perc(n_too_short, n_total)))
        print("No direction: %d (%.2f%%)" % (n_no_direction, get_perc(n_no_direction, n_total)))
        print("Large edit distance: %d (%.2f%%)" % (n_large_ed, get_perc(n_large_ed, n_total)))
        print("Invalid cell: %d (%.2f%%)" % (n_invalid_cell, get_perc(n_invalid_cell, n_total)))
        print("Pass: %d (%.2f%%)" % (n_pass, get_perc(n_pass, n_total)))
        
    if True:
        with open(os.path.join(outdir, "stats.tsv"), "w+") as fw:
            fw.write("Total\tTooShort\tNoDirection\tLargeED\tInvalidCell\tPass\n")
            fw.write("\t".join(map(str, [n_total, n_too_short, n_no_direction, n_large_ed, n_invalid_cell, n_pass])) + "\n")
    
    if True:  
        with open(os.path.join(outdir, "reads.tsv"), "w+") as fw:
            fw.write("\t".join(["Cell", "Reads", "Forward", "Reverse"]) + "\n")
            for cell in cell_list:
                v1, v2, v3 = cell_reads_counter[cell]
                fw.write("\t".join(map(str, [cell, v1, v2, v3])) + "\n")
            v1, v2, v3 = cell_reads_counter["unclassified"]
            fw.write("\t".join(map(str, ["unclassified", v1, v2, v3])) + "\n")
          
    if True:  
        with open(os.path.join(outdir, "edit_distance.tsv"), "w+") as fw:
            fw.write("ED\tCount\tRatio\tCumulative\n")
            r0 = 0
            for k, v in sorted(ed_counter.items()):
                r = v / sum(ed_counter.values())
                r0 += r
                fw.write("%d\t%d\t%f\t%f\n" % (k, v, r, r0))

USAGE = """
    %prog [options] <input.fq.gz> <matrix.gz> <config.tsv> <outdir>
    
<input.fq.gz>: PATH of fastq file to be splitted.
<matrix.gz>: PATH of matrix file generated by FBILR.
<config.tsv>: Dual barcodes config file. Tab-delimited for each line: cell, 1st barcode, 2nd barcode.
<outdir>: Output directory.
"""

def main():    
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-k", "--keep", dest="keep", action="store_true", default=False, 
                      help="Output unclassified reads. [default: %default]")
    parser.add_option("-e", "--edit-distance", dest="ed", type="int", default=5, metavar="INT",
                      help="Max edit distance. [default: %default]")
    parser.add_option("-l", "--length", dest="length", type="int", default=400, metavar="INT",
                      help="Minimum length. [%default]")
    options, args = parser.parse_args()    
    keep_unclassfied = options.keep
    max_edit_distance = options.ed
    min_length = options.length
    fastq, matrix, config, outdir = args

    demultiplexing(f_fastq=fastq, 
                   f_matrix=matrix, 
                   f_config=config,
                   outdir=outdir,
                   max_edit_distance=max_edit_distance,
                   keep_unclassified=keep_unclassfied,
                   min_length=min_length)
    

if __name__ == '__main__':
    main()
