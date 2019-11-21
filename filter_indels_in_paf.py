import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('paf', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help = 'Input paf file (default: stdin')
    parser.add_argument('--max_indel', type = int, default = 20, help = 'Maximum allowed indel size (default: 20)')
    args = parser.parse_args()

    num_filtered = 0
    for line in args.paf:
        fields = line.split('\t')
        query_name = fields[0]
        query_length = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        query_dir = fields[4]
        target_name = fields[5]
        target_start = int(fields[7])
        target_end = int(fields[8])
        aln_block_length = int(fields[10])
        tags = fields[12:]
        for field in tags:
            if field[:5] == 'cs:Z:':
                cigar2 = field[5:]
            if field[:5] == 'cg:Z:':
                cigar1 = field[5:]
        matches = re.findall(r'[0-9]+[A-Z]', cigar1)
        longest_indel = -1
        for match in matches:
        	length = int(match[:-1])
        	operation = match[-1]
        	if (operation == "I" or operation == "D") and length > longest_indel:
        		longest_indel = length
        if longest_indel > args.max_indel:
        	num_filtered += 1
        else:
        	print(line.strip())
    print("Done. Filtered %d paf lines." % (num_filtered), file=sys.stderr)


if __name__ == '__main__':
    main()