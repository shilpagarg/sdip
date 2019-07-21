i=11
while read bed; do 
    ((i=$i + 1));
    chr=`echo $bed | awk '{print $1}'`;
    start=`echo $bed | awk '{print $2}'`;
    ((start=$start - 50000));
    end=`echo $bed | awk '{print $3}'`;
    ((end=$end + 50000));
    samtools view -b -@32 ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/GRCh38_no_alt_analysis/HG002.15kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam ${chr}:${start}-${end} >  r${i}.bam;
    bamToFastq -i r${i}.bam -fq r${i}.fastq;
    seqtk seq -A r${i}.fastq > r${i}.fasta;
    awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' r${i}.fasta | awk '!seen[$1]++' | sed "s/\t/\n/g" > r${i}.rd.fasta
    #cat unmapped.fasta r${i}.rd.fasta > r${i}.full.fasta
    python paftest/haplotype.py r${i}.rd.fasta -o r${i}.gfa -t32 -l r${i}.list 2> r${i}.log;
done < $1
