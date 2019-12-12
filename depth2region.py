import argparse
from math import sqrt
import sys

def readBins(bin_file):
    with open(bin_file, 'r') as b:
        for line in b:
            fields = line.rstrip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            depth = float(fields[3])
            yield (chrom, start, end, depth)

def compute_depth_mean(bin_file):
    total_bases = 0
    total_depth = 0
    for chrom, start, end, depth in readBins(bin_file):
        total_bases += end - start
        total_depth += (end - start) * depth
    return total_depth, total_bases

def compute_depth_sd(bin_file, mean):
    total_bases = 0
    total_squared_deviation = 0
    for chrom, start, end, depth in readBins(bin_file):
        total_bases += end - start
        total_squared_deviation += (depth - mean) * (depth - mean) * (end - start)
    sd = sqrt(total_squared_deviation / (total_bases - 1))
    return sd

def outputRegion(region, form = 'bed'):
    if form == 'bed':
        print(region[0]+'\t'+str(region[1])+'\t'+str(region[2]))
    elif form == 'region':
        print(region[0]+':'+str(region[1])+'-'+str(region[2]))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bins', type = str, metavar = 'FILE', help = 'Output file from "bin_depths.py".')
    parser.add_argument('-l', '--min_length', type = int, metavar = 'INT', default = 10000, help = 'Minimum region length [10000] to consider')
    parser.add_argument('-t', '--tolerance', type = int, default = 2, help = 'While going through each bin, allow [5] bins that have coverage lower than threshold. (i.e. A short drop in coverage')
    parser.add_argument('-f', '--format', type = str, metavar = '[BED|region]', default = 'bed', help = 'Output format. bed refers to "chr\\tfrom\\tto" for each line; region refers to "chr:from-to" for each line')
    args = parser.parse_args()
    
    total_depth, total_bases = compute_depth_mean(args.bins)
    mean_depth = total_depth / total_bases
    print("Mean:", mean_depth, file=sys.stderr)
    
    sd = compute_depth_sd(args.bins, mean_depth)
    print("SD:", sd, file=sys.stderr)
    
    min_cov = mean_depth + 3*sd
    print("Min coverage:", min_cov, file=sys.stderr)
    
    in_region = False
    region_chrom = ''
    region_start = -1
    region_end = -1
    regions = []
    
    for chrom, start, end, depth in readBins(args.bins):
        if depth >= min_cov:
            if in_region:
                if chrom == region_chrom:
                    region_end = end
                else:
                    regions.append((region_chrom, region_start, region_end))
                    region_chrom = chrom
                    region_start = start
                    region_end = end
            else:
                in_region = True
                region_chrom = chrom
                region_start = start
                region_end = end
        else:
            if in_region:
                if chrom == region_chrom:
                    if (end - region_end) > args.tolerance:
                        regions.append((region_chrom, region_start, region_end))
                        in_region = False
                        region_chrom = ''
                        region_start = -1
                        region_end = -1
                else:
                    regions.append((region_chrom, region_start, region_end))
                    in_region = False
                    region_chrom = ''
                    region_start = -1
                    region_end = -1
    if in_region:
        regions.append((region_chrom, region_start, region_end))
    
    print("Found %d regions (%d larger than threshold)" % (len(regions), len([r for r in regions if (r[2] - r[1]) >= args.min_length])), file=sys.stderr)
    region_chrom = ''
    region_start = -1
    region_end = -1
    for chrom, start, end in regions:
        if (end - start) >= args.min_length:
            if region_chrom == '':
                region_chrom = chrom
                region_start = max(0, start - 15000)
                region_end = end + 15000
            else:
                if region_chrom == chrom and ((start - 15000) - region_end) < 10000:
                    region_end = end + 15000
                else:
                    outputRegion((region_chrom, region_start, region_end), form = args.format)
                    region_chrom = chrom
                    region_start = max(0, start - 15000)
                    region_end = end + 15000
    if region_chrom != '':
        outputRegion((region_chrom, region_start, region_end), form = args.format)

if __name__ == '__main__':
    main()
