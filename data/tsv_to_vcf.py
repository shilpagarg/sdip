import sys
import argparse

def parseArguments(args):
    """Set up the command-line parser and call it on the command-line arguments to the program.

    arguments:
    args -- command-line arguments to the program"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='Reads BED file of SVs from VISOR and converts it to VCF')
    parser.add_argument('bed', type=argparse.FileType('r'), help='BED file from VISOR')
    parser.add_argument('contigs', type=argparse.FileType('r'), help='contig definition in reference genome (generate with "awk \'{print "##contig=<ID="$1",length="$2">"}\' genome.fa.fai")')
    return parser.parse_args()

def print_vcf_header(contigs_file):
    print("##fileformat=VCFv4.2")
    for line in contigs_file:
        print(line.strip())
    print("""
##INFO=<ID=END,Number=1,Type=Integer,Description="end position of the variant described in this record">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE""")

def parse_records(bed):
    num_del = 0
    num_ins = 0
    for line in bed:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        variant_type = fields[3]
        info = fields[4]
        genotype = fields[6]
        if variant_type == "deletion":
            num_del += 1
            print("%s\t%d\t%s\tN\t<DEL>\t.\tPASS\tEND=%d;SVTYPE=DEL;SVLEN=%d\tGT\t%s" % (chrom, start, "DEL_"+str(num_del), end, -(end-start), genotype))
        elif variant_type == "insertion":
            num_ins += 1
            print("%s\t%d\t%s\tN\t<INS>\t.\tPASS\tEND=%d;SVTYPE=INS;SVLEN=%d\tGT\t%s" % (chrom, start, "INS_"+str(num_ins), end, len(info), genotype))

def main(args):
    options = parseArguments(args)

    print_vcf_header(options.contigs)
    parse_records(options.bed)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
