#Total size and mean size of collapsed regions
cat /project/heller_d-data/shilpa/NA12878/sda/coverage/na12878.collapses.bed | awk '{SUM+=$3-$2} END {print SUM, NR, SUM/NR}'
#Total size of aligned collapsed regions (alignment to hg38)
cat /project/heller_d-data/shilpa/NA12878/regions/segdups/aln_to_hg38/*.bed | bedtools sort -i - | bedtools merge | awk '{SUM+=$3-$2} END {print SUM}'
#Number of collapsed regions
cat *.bed | wc -l
#Number of collapsed regions overlapping annotated SDs
cat *.bed | bedtools sort -i - | bedtools intersect -a - -b /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -u | wc -l
#Total size of collapsed regions overlapping annotated SDs
cat *.bed | bedtools sort -i - | bedtools intersect -a - -b /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -u | bedtools merge | awk '{SUM+=$3-$2} END {print SUM}'
