#generate genome from segdup regions
cd /project/heller_d-data/shilpa/simulated/data/

mkdir genome
mkdir svs
mkdir reads

cd /project/heller_d-data/shilpa/simulated/data/genome
bedtools slop -i /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/segdups.hg38.sorted.merged.head100.bed -g /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai -b 50000 | bedtools merge | awk '{print $1":"$2"-"$3}' > segdups.hg38.merged.slop50k.txt
samtools faidx -r segdups.hg38.merged.slop50k.txt /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa > genome.simulated.fa
#produce N regions
seqtk cutN -n 200 -gp100 genome.simulated.fa > genome.simulated.n.bed
#remove N regions
seqtk cutN -n 200 -p100 genome.simulated.fa | seqtk seq -L 1000 > genome.simulated.no-n.fa
samtools faidx genome.simulated.no-n.fa
cut -f1,2 genome.simulated.no-n.fa.fai > genome.simulated.no-n.genome

#simulate SVs
cd /project/heller_d-data/shilpa/simulated/data
#script source: https://raw.githubusercontent.com/davidebolo1993/VISOR/master/scripts/randomregion.r
Rscript randomregion.r -d genome/genome.simulated.no-n.genome -n '25:25:25:25:7362:7362' --hap1 svs/svs.simulated.hap1.bed --hap2 svs/svs.simulated.hap2.bed
sortBed -i svs/svs.simulated.hap1.bed > svs/svs.simulated.hap1.sorted.bed
sortBed -i svs/svs.simulated.hap2.bed > svs/svs.simulated.hap2.sorted.bed
conda activate visorenv
VISOR HACk -g genome/genome.simulated.no-n.fa -bed svs/svs.simulated.hap1.sorted.bed svs/svs.simulated.hap2.sorted.bed -o genome/hack.out

#simulate reads
cd /project/heller_d-data/shilpa/simulated/data/reads
pbsim --seed 10000 --prefix sim.hap1 --data-type CCS --depth 20 --length-min 500 --length-max 20000 --sample-fastq /project/heller_d-data/shilpa/CHM13/reads_ccs/allreads.fastq ../genome/hack.out/h1.fa
pbsim --seed 10000 --prefix sim.hap2 --data-type CCS --depth 20 --length-min 500 --length-max 20000 --sample-fastq /project/heller_d-data/shilpa/CHM13/reads_ccs/allreads.fastq ../genome/hack.out/h2.fa
cat sim.hap1_*.fastq | awk 'NR%4==1 {printf("%s_HAP1\n",$0)} NR%4!=1 {print}' > sim.hap1.pooled.fastq
cat sim.hap2_*.fastq | awk 'NR%4==1 {printf("%s_HAP2\n",$0)} NR%4!=1 {print}' > sim.hap2.pooled.fastq
cat sim.hap1.pooled.fastq sim.hap2.pooled.fastq > sim.pooled.fastq
samtools faidx -f sim.pooled.fastq

#convert bed to vcf
cd /project/heller_d-data/shilpa/simulated/data/svs
awk '{print "##contig=<ID="$1",length="$2">"}' ../genome/genome.simulated.no-n.fa.fai > svs.simulated.sorted.contigs.txt
python ../tsv_to_vcf.py svs.simulated.hap1.sorted.bed svs.simulated.sorted.contigs.txt | bcftools sort -Oz > svs.simulated.sorted.vcf.gz
tabix -p vcf svs.simulated.sorted.vcf.gz
