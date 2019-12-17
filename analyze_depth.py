import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('overlaps', type = str, metavar = 'FILE', help = 'overlap.bed.')
    args = parser.parse_args()
    
    bed = pd.read_csv( args.overlaps, sep = "\t", header=None, names=['contig', 'start', 'end',"coverage", "cr"])

    # I want to eliminte the really low or really high coverage things because they are probably
    # not assembled correctly and then assecess what the mean and standard deviation is
    top = bed.coverage.quantile(.90)
    bot = bed.coverage.quantile(.10)
    
    # save stats like mean coverage 
    stats = bed["coverage"][ (bed.coverage < top) & ( bed.coverage > bot) ].describe()
    
    # filter for high coverage regsion
    mincov = stats["mean"]  + 3 * np.sqrt(stats["mean"])

    print("Mean coverage:", stats["mean"])
    print("Standard deviation:", np.sqrt(stats["mean"]))
    print("Mean + 3 * standard dev.:", mincov)

if __name__ == '__main__':
    main()
