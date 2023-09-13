#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    infile, outfile = sys.argv[1:]
    dat = pd.read_csv(infile, sep="\t", index_col=0)
    name = infile.split("/")[-1]
    total = dat["Reads"].sum()
    uncls = dat.loc["unclassified"]["Reads"]
    ratio = uncls / total
    dat = dat.loc[list(filter(lambda item: item != "unclassified", dat.index))]
    mean = dat["Reads"].mean()
    median = dat["Reads"].median()
    std = dat["Reads"].std()
    xs = np.arange(len(dat))
    xticks = [cell.split(".")[-1] for cell in dat.index]
    ys1 = dat["Forward"]
    ys2 = dat["Reverse"]
    ylim = max(ys1 + ys2) * 1.5
    plt.figure(figsize=(max(10, 1 + len(xs) * 0.2), 4))
    plt.bar(xs, ys1, label="Forward", color="C0")
    plt.bar(xs, ys2, bottom=ys1, label="Reverse", color="C1")
    plt.axhline(mean, color="grey", ls="--")
    plt.xlim(min(xs) - 1, max(xs) + 1)
    plt.xticks(xs, xticks, rotation=90)
    plt.ylim(0, ylim)
    plt.text(0, ylim * 0.95, name)
    plt.text(0, ylim * 0.90, "Mean = %s" % format(int(mean), ","))
    plt.text(0, ylim * 0.85, "Median = %s" % format(int(median), ","))
    plt.text(0, ylim * 0.80, "Std = %s" % format(int(std), ","))
    plt.text(0, ylim * 0.75, "Unclassified = %.2f%%" % (ratio * 100))
    plt.ylabel("Number of read")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    
    
if __name__ == '__main__':
    main()
    