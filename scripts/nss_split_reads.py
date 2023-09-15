#!/usr/bin/env python
import os
import optparse
from collections import defaultdict
from pygz import PigzFile


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
    ed_counter = defaultdict(int)
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
    n_no_direction = 0
    n_large_ed = 0
    n_invalid_cell = 0
    n_pass = 0
    
    for read, row in zip(load_fastq(f_fastq), load_matrix(f_matrix)):
        name, seq, qua = read
        
        assert name == row[0]
        
        bc1, direct1, loc1, x1, y1, ed1 = row[2:8]
        bc2, direct2, loc2, x2, y2, ed2 = row[8:14]
        
        barcode_counter[tuple(row[2:])] += 1
        
        fw = None
        n_total += 1
        
        cell = None
        
        if len(seq) <= min_len:
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
                if ed > max_ed:
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
        
    print("Total reads: %d (%.2f%%)" % (n_total, calculate_percentage(n_total, n_total)))
    print("Too short: %d (%.2f%%)" % (n_too_short, calculate_percentage(n_too_short, n_total)))
    print("No direction: %d (%.2f%%)" % (n_no_direction, calculate_percentage(n_no_direction, n_total)))
    print("Large edit distance: %d (%.2f%%)" % (n_large_ed, calculate_percentage(n_large_ed, n_total)))
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
            
    print("-" * 80)
    print("ED\tCount\tRatio\tCumulativeRatio")
    print("-" * 80)
    r0 = 0
    for k, v in sorted(ed_counter.items()):
        r = v / sum(ed_counter.values())
        r0 += r
        print("%d\t%d\t%f\t%f" % (k, v, r, r0))
    

def main():
    
    parser = optparse.OptionParser()

    parser.add_option("-m", "--matrix", dest="matrix", metavar="PATH",
                      help="PATH of matrix2 file.")
    parser.add_option("-b", "--barcode-config", dest="barcode", metavar="PATH",
                      help="Dual barcodes config file. Tab-delimited for each line: cell, p7 barcode, p5 barcode.")
    parser.add_option("-f", "--fastq", dest="fastq", metavar="PATH",
                      help="PATH of fastq file to be splitted.")
    parser.add_option("-k", "--keep", dest="keep", action="store_true", default=False, 
                      help="Output unclassified reads. [default: %default]")
    parser.add_option("-e", "--edit-distance", dest="ed", type="int", default=5, metavar="INT",
                      help="Max edit distance. [default: %default]")
    parser.add_option("-l", "--length", dest="length", type="int", default=400, metavar="INT",
                      help="Minimum length. [%default]")
    parser.add_option("-o", "--outdir", dest="outdir", default="./", metavar="DIR",
                      help="Output directory. [default: %default]")
    # parser.add_option("-s", "--stats", dest="stats", 
    #                   help="PATH to statistic. [default: %default]")
    parser.add_option("-r", "--reads", dest="reads", metavar="PATH",
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
