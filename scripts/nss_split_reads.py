#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time : May 2, 2022
# @Author : Zonggui Chen
# @Email : chenzonggui@stu.pku.edu.cn

# This script splits reads into single-cell reads by the metrics of FBILR. 
# This script is designed for Nanopore-based Strand-seq which contains two barcodes each reads.
# The NanoNASCseq.xls in the working directory are required to fetch the cell and correspond barcode.

import sys
import os
import optparse
from collections import defaultdict
import pandas as pd
from pygz import PigzFile
from fbilr.reader import Matrix2Reader

# Spliting Nanopore-based Strand-seq reads to different files.


DEBUG = False


def load_fastq(path):
    assert path.endswith(".gz")
    with PigzFile(path, "rt") as f:
        name = None
        seq = None
        qua = None
        for i, line in enumerate(f):
            if DEBUG and i >= 400000:
                break
            j = i % 4
            if j == 0:
                name = line[1:-1].split()[0]
            elif j == 1:
                seq = line[:-1]
            elif j == 3:
                qua = line[:-1]
                yield name, seq, qua
                                
                
def load_matrix(path):
    with PigzFile(path) as f:
        for line in f:
            row = line.strip("\n").split("\t")
            for i in [1, 5, 6, 7, 11, 12, 13]:
                row[i] = int(row[i])
            yield row


def calculate_percentage(n1, n2):
    if n2 == 0:
        return 0
    else:
        return n1 * 100 / n2
    

def demultiplexing(f_fastq, f_matrix, f_barcode, 
                   f_outdir, f_reads=None, 
                   max_ed=5, min_len=400, 
                   keep_unclassified=False, trim_outer=True):
    
    # output directory
    
    if not os.path.exists(f_outdir):
        os.mkdir(f_outdir)
        
    # dual barcode pair
    
    cell_list = []
    cell_reads_counter = dict()
    bc2cell = dict()  # barcode to cell
    fws = dict()
    
    with open(f_barcode) as f:
        for line in f:
            cell, bc1, bc2 = line.strip("\n").split("\t")[:3]
            cell_list.append(cell)
            
            bc2cell[(bc1, bc2)] = cell  # (p7, p5)
            fw_f = open(os.path.join(f_outdir, "%s_F.fastq" % cell), "w+")
            fw_r = open(os.path.join(f_outdir, "%s_R.fastq" % cell), "w+")
            fws[cell] = [fw_f, fw_r]
            cell_reads_counter[cell] = [0, 0, 0]  # total forward reverse
            
    cell_reads_counter["unclassified"] = [0, 0, 0]
    if keep_unclassified:
        unfw = open(os.path.join(f_outdir, "unclassified.fastq") , "w+")
    else:
        unfw = None
    
    barcode_counter = defaultdict(int)
    
    n_total = 0
    n_too_short = 0
    n_no_barcode = 0
    n_invalid_distance = 0
    n_invalid_cell = 0
    n_pass = 0
    
    for read, row in zip(load_fastq(f_fastq), load_matrix(f_matrix)):
        name, seq, qua = read
        
        assert name == row[0]
        
        bc1, direction1, location1, x1, y1, ed1 = row[2:8]
        bc2, direction2, location2, x2, y2, ed2 = row[8:14]
        
        barcode_counter[tuple(row[2:])] += 1
        
        fw = None
        n_total += 1
        
        cell = None
        
        if len(seq) > min_len:
            if ed1 <= max_ed and ed2 <= max_ed:
                read_direction = None
                if direction1 == "F" and location1 == "H" and direction2 == "R" and location2 == "T":
                    read_direction = "F"
                elif direction1 == "R" and location1 == "T" and direction2 == "F" and location2 == "H":
                    read_direction = "R"
                if read_direction is None:
                    n_no_barcode += 1
                else:
                    # check distance
                    x, y = min(x1, x2), max(y1, y2)
                    dis1 = x
                    dis2 = len(seq) - y
                    if dis1 < 70 and dis2 < 70:
                        cell = bc2cell.get((bc1, bc2))
                        if cell is None:
                            n_invalid_cell += 1
                        else:
                            if trim_outer:
                                x, y = min(y1, y2), max(x1, x2)
                                seq = seq[x:y]
                                qua = qua[x:y]
                            cell_reads_counter[cell][0] += 1
                            if read_direction == "F":
                                fw = fws[cell][0]
                                cell_reads_counter[cell][1] += 1
                            else:
                                fw = fws[cell][1]
                                cell_reads_counter[cell][2] += 1
                            n_pass += 1
                    else:
                        n_invalid_distance += 1
            else:
                n_no_barcode += 1
        else:
            n_too_short += 1

 
        if fw is None:
            fw = unfw
            cell_reads_counter["unclassified"][0] += 1
            
        # reduce file size by simplify read name
        if fw:
            fw.write("@%s\n%s\n+\n%s\n" % (name, seq, qua))            
    
    for fw1, fw2 in fws.values():
        fw1.close()
        fw2.close()
    if unfw:
        unfw.close()
        
    print("Total reads: %d (%.2f%%)" % (n_total, calculate_percentage(n_total, n_total)))
    print("Too short: %d (%.2f%%)" % (n_too_short, calculate_percentage(n_too_short, n_total)))
    print("No barcode: %d (%.2f%%)" % (n_no_barcode, calculate_percentage(n_no_barcode, n_total)))
    print("Invalid distance: %d (%.2f%%)" % (n_invalid_distance, calculate_percentage(n_invalid_distance, n_total)))
    print("Invalid cell: %d (%.2f%%)" % (n_invalid_cell, calculate_percentage(n_invalid_cell, n_total)))
    print("Pass: %d (%.2f%%)" % (n_pass, calculate_percentage(n_pass, n_total)))
                
    if f_reads:
        with open(f_reads, "w+") as fw:
            fw.write("\t".join(["Cell", "Reads", "Forward", "Reverse"]) + "\n")
            for cell in cell_list:
                v1, v2, v3 = cell_reads_counter[cell]
                fw.write("\t".join(map(str, [cell, v1, v2, v3])) + "\n")
            v1, v2, v3 = cell_reads_counter["unclassified"]
            fw.write("\t".join(map(str, ["unclassified", v1, v2, v3])) + "\n")
    

def main():
    
    parser = optparse.OptionParser()

    parser.add_option("-m", "--matrix", dest="matrix", 
                      help="PATH of matrix2 file.")
    parser.add_option("-b", "--barcode-cofnig", dest="barcode", 
                      help="Dual barcodes config file. Tab-delimited for each line: cell, p7 barcode, p5 barcode.")
    parser.add_option("-f", "--fastq", dest="fastq", 
                      help="PATH of fastq file to be splitted.")
    parser.add_option("-k", "--keep", dest="keep", action="store_true", default=False, 
                      help="Output unclassified reads. [default: %default]")
    parser.add_option("-e", "--edit-distance", dest="ed", type="int", default=5, 
                      help="Max edit distance. [default: %default]")
    parser.add_option("-l", "--length", dest="length", type="int", default=400, 
                      help="Minimum length. [%default]")
    parser.add_option("-o", "--outdir", dest="outdir", default="./", 
                      help="Output directory. [default: %default]")
    # parser.add_option("-s", "--stats", dest="stats", 
    #                   help="PATH to statistic. [default: %default]")
    parser.add_option("-r", "--reads", dest="reads", 
                      help="PATH to cell reads summary. [default: %default]")
    options, args = parser.parse_args()
    # print(options)
    
    if len(args) != 0:
        parser.print_usage()
        exit(1)

    demultiplexing(f_fastq=options.fastq, 
                   f_matrix=options.matrix, 
                   f_barcode=options.barcode,
                   f_outdir=options.outdir,
                   f_reads=options.reads, 
                   max_ed=options.ed, 
                   min_len=options.length)
    

if __name__ == '__main__':
    main()
