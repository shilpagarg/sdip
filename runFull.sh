r=$1
#samtools view r${r}.bam | cut -d$'\t' -f1 | sort | uniq > r${r}.reads
#python ../tmp/hg002-asm-r3-pg0.1.5.3/getOv.py ../tmp/hg002-asm-r3-pg0.1.5.3/0-seqdb/seq_dataset.idx ../tmp/hg002-asm-r3-pg0.1.5.3/3-asm/split/ r${r}.reads | sort | uniq > r${r}.full.reads
cd split
#rm *
#cp ../r${r}.full.reads ./
#split -l 200 r${r}.full.reads
#rm r${r}.full.reads
ls * | parallel 'seqtk subseq ../../hg002.pacbioccs.fasta {} > {}.fasta'
cat *.fasta > ../r${r}.full.2.fasta
cd ../
python paftest/haplotype.py r${r}.full.2.fasta -o r${r}.full.2.gfa -t30 2> r${r}.full.2.log
