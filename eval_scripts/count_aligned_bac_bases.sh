cat /project/heller_d-data/shilpa/CHM13/regions/contigs/pooled.t5.b5.d2.polished.fa /project/heller_d-data/shilpa/CHM13/ref/corrected.fasta > contigs.plus.assm.fa
minimap2 -I 8G -t 10 --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 contigs.plus.assm.fa /project/heller_d-data/shilpa/CHM13/bacs/bacs.fa | samtools view -F 2308 -u - | samtools sort - > bacs.to.combined.bam
python3 /project/pacbiosv/code/WHdenovo/paftest/samIdentity.py --header --mask /project/heller_d-data/shilpa/CHM13/regions/eval/bams/bacs.in.sd.names bacs.to.combined.bam > bacs.to.combined.tbl
awk '{SUM+=$3-$2} END {print SUM}' bacs.to.combined.tbl 
/project/heller_d-data/shilpa/CHM13/combined_test> awk '{SUM+=$2} END {print SUM}' ../bacs/bacs.fa.fai

minimap2 -I 8G -t 10 --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 /project/heller_d-data/shilpa/CHM13/ref/corrected.fasta /project/heller_d-data/shilpa/CHM13/bacs/bacs.fa | samtools view -F 2308 -u - | samtools sort - > bacs.to.corrected.bam
python3 /project/pacbiosv/code/WHdenovo/paftest/samIdentity.py --header --mask /project/heller_d-data/shilpa/CHM13/regions/eval/bams/bacs.in.sd.names bacs.to.corrected.bam > bacs.to.corrected.tbl
awk '{SUM+=$3-$2} END {print SUM}' bacs.to.corrected.tbl

minimap2 -I 8G -t 10 --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 /project/heller_d-data/genomes/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /project/heller_d-data/shilpa/CHM13/bacs/bacs.fa | samtools view -F 2308 -u - | samtools sort - > bacs.to.hg38.bam
python3 /project/pacbiosv/code/WHdenovo/paftest/samIdentity.py --header --mask /project/heller_d-data/shilpa/CHM13/regions/eval/bams/bacs.in.sd.names bacs.to.hg38.bam > bacs.to.hg38.tbl
awk '{SUM+=$3-$2} END {print SUM}' bacs.to.hg38.tbl
