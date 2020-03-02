bedtools bamtobed -i /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bam > /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bed

##########
#Unmerged#
##########

#Retrieve segdup regions that are covered or uncovered by base assembly
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bed -f 1.0 -v > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bed -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.ref.bed

#Retrieve segdup regions that are covered by contigs
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5000.d2/bams/polished.diploid.to.hg38.bam -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sdip.bed
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/sda/bams/polished.haploid.to.hg38.bam -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sda.bed

bedtools intersect -a /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed -b /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sdip.bed -f 1.0 -r -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.covered.sdip.bed
bedtools intersect -a /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed -b /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sda.bed -f 1.0 -r -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.covered.sda.bed

#Total number of SDs
wc -l /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.bed
#Number of SDs uncovered by base assembly
wc -l /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed
#Span of SDs uncovered by base assembly
cat /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed | bedtools sort | bedtools merge | awk '{SUM+=$3-$2} END {print SUM}'

#Number of SDs uncovered by base assembly but covered by SDip
wc -l /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.covered.sdip.bed
#Span of SDs uncovered by base assembly but covered by SDip
cat /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.covered.sdip.bed | bedtools sort | bedtools merge | awk '{SUM+=$3-$2} END {print SUM}'

#Number of SDs covered by SDip
wc -l /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sdip.bed
#Span of SDs covered by SDip
cat /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sdip.bed | bedtools sort | bedtools merge | awk '{SUM+=$3-$2} END {print SUM}'

#Size range of SDs covered by SDip and SDs uncovered by base assembly
cat /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.sdip.bed | awk '{print $0, $3-$2}' | sort -k7,7nr | less
cat /project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed | awk '{print $0, $3-$2}' | sort -k7,7nr | less

########
#Merged#
########

#Retrieve segdup regions that are covered or uncovered by base assembly
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.merged.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bed -f 1.0 -v > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.uncovered.ref.bed
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.merged.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/bams/ref.to.hg38.bed -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.covered.ref.bed

#Retrieve segdup regions that are covered by contigs
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.merged.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5000.d2/bams/polished.diploid.to.hg38.bam -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.covered.sdip.bed
bedtools intersect -a /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.merged.bed -b /project/heller_d-data/shilpa/NA12878/regions/eval/sda/bams/polished.haploid.to.hg38.bam -f 1.0 -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.covered.sda.bed

bedtools intersect -a /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.uncovered.ref.bed -b /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.covered.sdip.bed -f 1.0 -r -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.uncovered.covered.sdip.bed
bedtools intersect -a /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.uncovered.ref.bed -b /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.covered.sda.bed -f 1.0 -r -u > /project/heller_d-data/shilpa/NA12878/resolved_ucsc/merged/segdups.uncovered.covered.sda.bed

