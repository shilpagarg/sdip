#PACBIO
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/pacbiosv/data/HG002_pacbio_ccs/fastq/allreads.fastq.fai
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/heller_d-data/shilpa/NA12878/reads_pacbio/na12878.pacbiocss.fastq.fai
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/heller_d-data/shilpa/CHM13/reads_ccs/allreads.fastq.fai

cat /project/pacbiosv/data/HG002_pacbio_ccs/fastq/allreads.fastq.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'
cat /project/heller_d-data/shilpa/NA12878/reads_pacbio/na12878.pacbiocss.fastq.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'
cat /project/heller_d-data/shilpa/CHM13/reads_ccs/allreads.fastq.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'

#NANOPORE
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/pacbiosv/data/HG002_nanopore_ultralong/fastq/reads.fastq.fai
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/pacbiosv/data/NA12878_reference/raw/rel5-guppy-0.3.0-chunk10k.fastq.fai
awk '{SUM+=$2} END {print "Bases:", SUM; print "Coverage:", SUM/3000000000; print "Reads:", NR; print "Mean length", SUM/NR}' /project/heller_d-data/shilpa/CHM13/reads_nano/rel3.fastq.gz.fai

cat /project/pacbiosv/data/HG002_nanopore_ultralong/fastq/reads.fastq.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'
cat /project/pacbiosv/data/NA12878_reference/raw/rel5-guppy-0.3.0-chunk10k.fastq.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'
cat /project/heller_d-data/shilpa/CHM13/reads_nano/rel3.fastq.gz.fai | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print lastcontig } { lastsize=$2; lastcontig=$1 }'
