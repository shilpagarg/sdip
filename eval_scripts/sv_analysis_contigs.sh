cd /project/heller_d-data/shilpa/NA12878/regions/svim/t5.b5000.d2/

#View lengths and genotypes
bcftools query -i 'SVTYPE=="DEL"' -f '[ %GT]\n' contigs.diploid/variants.norm.vcf | sort -k 1,1n | uniq -c
bcftools query -i 'SVTYPE=="DEL"' -f '%SVLEN\n' contigs.diploid/variants.norm.vcf | sort -k 1,1n | less
bcftools query -i 'SVTYPE=="INS"' -f '[ %GT]\n' contigs.diploid/variants.norm.vcf | sort -k 1,1n | uniq -c
bcftools query -i 'SVTYPE=="INS"' -f '%SVLEN\n' contigs.diploid/variants.norm.vcf | sort -k 1,1n | less


cd /project/heller_d-data/shilpa/NA12878/regions/svim
#Find overlap regions
bedtools intersect -a <(bedtools bamtobed -i ../eval/bams/bacs.to.hg38.bam) -b <(bedtools bamtobed -i ../eval/t5.b5000.d2/bams/polished.diploid.to.hg38.bam) | bedtools sort | bedtools merge > overlap.sdip.bacs.bed
bedtools intersect -a <(bedtools bamtobed -i ../eval/bams/bacs.to.hg38.bam) -b <(bedtools bamtobed -i ../eval/sda/bams/polished.haploid.to.hg38.bam) | bedtools sort | bedtools merge > overlap.sda.bacs.bed
bedtools intersect -a <(bedtools bamtobed -i ../eval/bams/bacs.to.hg38.bam) -b <(bedtools bamtobed -i ../eval/t5.b5000.d2/bams/polished.diploid.to.hg38.bam) | bedtools intersect -a - -b <(bedtools bamtobed -i ../eval/sda/bams/polished.haploid.to.hg38.bam) | bedtools sort | bedtools merge > overlap.sdip.sda.bacs.bed

#Compute total span of overlap
awk '{SUM+=$3-$2} END {print SUM}' overlap.sdip.sda.bacs.bed

#Subset calls in overlap regions
bcftools view -R overlap.sdip.bacs.bed t5.b5000.d2/contigs.diploid/variants.norm.vcf.gz > t5.b5000.d2/contigs.diploid/variants.norm.overlap.sdip.bacs.vcf
bcftools view -R overlap.sdip.bacs.bed t5.b5000.d2/bacs.haploid/variants.norm.vcf.gz > t5.b5000.d2/bacs.haploid/variants.norm.overlap.sdip.bacs.vcf

bcftools view -R overlap.sda.bacs.bed sda/contigs.haploid/variants.norm.vcf.gz > sda/contigs.haploid/variants.norm.overlap.sda.bacs.vcf
bcftools view -R overlap.sda.bacs.bed sda/bacs.haploid/variants.norm.vcf.gz > sda/bacs.haploid/variants.norm.overlap.sda.bacs.vcf

bcftools view -R overlap.sdip.sda.bacs.bed t5.b5000.d2/contigs.diploid/variants.norm.vcf.gz > t5.b5000.d2/contigs.diploid/variants.norm.overlap.sdip.sda.bacs.vcf
bcftools view -R overlap.sdip.sda.bacs.bed sda/contigs.haploid/variants.norm.vcf.gz > sda/contigs.haploid/variants.norm.overlap.sdip.sda.bacs.vcf
bcftools view -R overlap.sdip.sda.bacs.bed t5.b5000.d2/bacs.haploid/variants.norm.vcf.gz > t5.b5000.d2/bacs.haploid/variants.norm.overlap.sdip.sda.bacs.vcf

#Prepare VCF lists
echo /project/heller_d-data/shilpa/NA12878/regions/svim/t5.b5000.d2/contigs.diploid/variants.norm.overlap.sdip.bacs.vcf > vcfs.sdip.bacs.fofn
echo /project/heller_d-data/shilpa/NA12878/regions/svim/t5.b5000.d2/bacs.haploid/variants.norm.overlap.sdip.bacs.vcf >> vcfs.sdip.bacs.fofn

echo /project/heller_d-data/shilpa/NA12878/regions/svim/sda/contigs.haploid/variants.norm.overlap.sda.bacs.vcf > vcfs.sda.bacs.fofn
echo /project/heller_d-data/shilpa/NA12878/regions/svim/sda/bacs.haploid/variants.norm.overlap.sda.bacs.vcf >> vcfs.sda.bacs.fofn

echo /project/heller_d-data/shilpa/NA12878/regions/svim/t5.b5000.d2/contigs.diploid/variants.norm.overlap.sdip.sda.bacs.vcf > vcfs.sdip.sda.bacs.fofn
echo /project/heller_d-data/shilpa/NA12878/regions/svim/sda/contigs.haploid/variants.norm.overlap.sdip.sda.bacs.vcf >> vcfs.sdip.sda.bacs.fofn
echo /project/heller_d-data/shilpa/NA12878/regions/svim/t5.b5000.d2/bacs.haploid/variants.norm.overlap.sdip.sda.bacs.vcf >> vcfs.sdip.sda.bacs.fofn

#Merge calls
/project/pacbiosv/bin/SURVIVOR/Debug/SURVIVOR merge vcfs.sdip.bacs.fofn 0.5 1 1 0 1 40 merged.sdip.bacs.vcf
/project/pacbiosv/bin/SURVIVOR/Debug/SURVIVOR merge vcfs.sda.bacs.fofn 0.5 1 1 0 1 40 merged.sda.bacs.vcf
/project/pacbiosv/bin/SURVIVOR/Debug/SURVIVOR merge vcfs.sdip.sda.bacs.fofn 0.5 1 1 0 1 40 merged.sdip.sda.bacs.vcf

bcftools sort merged.sdip.bacs.vcf | bgzip -c > merged.sdip.bacs.vcf.gz
bcftools sort merged.sda.bacs.vcf | bgzip -c > merged.sda.bacs.vcf.gz
bcftools sort merged.sdip.sda.bacs.vcf | bgzip -c > merged.sdip.sda.bacs.vcf.gz
tabix -p vcf merged.sdip.bacs.vcf.gz
tabix -p vcf merged.sda.bacs.vcf.gz
tabix -p vcf merged.sdip.sda.bacs.vcf.gz

#Print SV sets
bcftools query -f '%SUPP_VEC\n' merged.sdip.bacs.vcf.gz | sort -k 1,1n | uniq -c
bcftools query -f '%SUPP_VEC\n' merged.sda.bacs.vcf.gz | sort -k 1,1n | uniq -c
bcftools query -f '%SUPP_VEC\n' merged.sdip.sda.bacs.vcf.gz | sort -k 1,1n | uniq -c

#Convert set sizes for Rscript
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' merged.sdip.sda.bacs.vcf | sed -e 's/\(.\)/\1 /g' > merged.sdip.sda.bacs.txt

#Produce Venn plot
Rscript /project/pacbiosv/code/sdip/eval_scripts/plot_sv_venn.R
