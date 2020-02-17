import networkx as nx
import pandas as pd
import numpy as np
import pysam
from collections import defaultdict, Counter

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

wildcard_constraints:
    region="\d+",
    tip_max_size="\d+",
    bubble_max_size="\d+",
    degree_max_size="\d+"

#############################
#Find similar segdup regions#
#############################

rule fetch_segdup_sequences:
    input:
        regions = config["regions"]["original"]["sda"],
        reference = config["genome"]
    output:
        "regions/segdups/r{region}.fa"
    shell:
        "reg=`awk '{{if (NR=={wildcards.region}) {{ print $1\":\"$2\"-\"$3 }} }}' {input.regions}`; samtools faidx {input.reference} ${{reg}} > {output}"

rule map_segdup_sequences:
    input:
        fasta = "regions/segdups/r{region}.fa",
        reference = config["genome"]
    output:
        "regions/segdups/r{region}.bam"
    threads: 10
    shell:
        "minimap2 -ax asm10 -t {threads} --secondary=yes -N 100 -p 0.80 {input.reference} {input.fasta} | samtools view -b | samtools sort > {output}"

rule bam_to_tsv:
    input:
        "regions/segdups/r{region}.bam"
    output:
        "regions/segdups/r{region}.tsv"
    shell:
        "bedtools bamtobed -i {input} | awk 'OFS=\"\\t\" {{ print {wildcards.region}, $0 }}' > {output}"

rule faidx:
    input:
        "{name}.fa"
    output:
        "{name}.fa.fai"
    shell:
        "samtools faidx {input}"

rule concat_regions:
    input:
        expand("regions/segdups/r{region}.tsv", region=all_segdup_regions)
    output:
        "segdups.similar.tsv"
    run:
        o = open(output[0], 'w')
        index = 0
        for p in input:
            with open(p, 'r') as f:
                lines = f.readlines()
            if len(lines) < 5:
                index += 1
                for l in lines:
                    fields = l.strip().split()                    
                    print("%d\t%s\t%s\t%s" % (index, fields[1], fields[2], fields[3]), file = o)
        o.close()

######################
#Recruit reads#
######################

rule recruit_pacbio_reads:
    input:
        regions = "segdups.similar.tsv",
        bam = "alignment_pacbio/pooled.sorted.bam",
        index = "alignment_pacbio/pooled.sorted.bam.bai"
    output:
        "regions/bams_pacbio/r{region}.bam"
    params:
        primary_mapq_threshold = 1 
    threads: 5
    shell:
        "reg=`awk '{{if ($1=={wildcards.region}) {{ print $2\":\"$3\"-\"$4 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{reg}} > {output}"

rule recruit_nanopore_reads:
    input:
        regions = "segdups.similar.tsv",
        bam = "alignment_nanopore/pooled.sorted.bam",
        index = "alignment_nanopore/pooled.sorted.bam.bai"
    output:
        "regions/bams_nanopore/r{region}.bam"
    params:
        primary_mapq_threshold = 30
    threads: 5
    shell:
        "reg=`awk '{{if ($1=={wildcards.region}) {{ print $2\":\"$3\"-\"$4 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{reg}} > {output}"

rule get_read_names:
    input:
        "regions/bams_{technology}/r{region}.bam"
    output:
        "regions/reads_{technology}/r{region}.reads"
    shell:
        "samtools view {input} | cut -d$'\t' -f1 | sort | uniq > {output}"

rule get_pacbio_fastq:
    input:
        allreads = config["reads"]["pacbio"],
        names = "regions/reads_pacbio/r{region}.reads"
    output:
        temp("regions/fastas_pacbio/r{region}.fastq")
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ] && [ $LINES -gt 0 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

rule get_nanopore_fastq:
    input:
        allreads = config["reads"]["nanopore"],
        names = "regions/reads_nanopore/r{region}.reads"
    output:
        temp("regions/fastas_nanopore/r{region}.fastq")
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ] && [ $LINES -gt 0 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

rule convert_fastq_to_fasta:
    input:
        "regions/fastas_{technology}/r{region}.fastq"
    output:
        "regions/fastas_{technology}/r{region}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule all_vs_all_alignment:
    input:
        fasta = "regions/fastas_pacbio/r{region}.fasta"
    output:
        paf = "regions/pafs/r{region}.fasta.paf"
    threads: 10
    shell:
        "minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t {threads} {input.fasta} {input.fasta} | sort -k8n -k6 > {output.paf}"

rule filter_paf_for_long_indels:
    input:
        paf = "regions/pafs/r{region}.fasta.paf"
    output:
        paf = "regions/pafs/r{region}.fasta.filtered.paf"
    shell:
        "python3 %s/filter_indels_in_paf.py {input.paf} --max_indel 50 > {output.paf}" % (config["tools"]["paftest"])

rule generate_graph_use_paf:
    input:
        fasta = "regions/fastas_pacbio/r{region}.fasta",
        paf = "regions/pafs/r{region}.fasta.filtered.paf"
    output:
        gfa = "regions/gfas/r{region}.reducted.gfa",
        contained = "regions/gfas/r{region}.contained.reads"
    threads: 3
    params:
        output = "regions/gfas/r{region}.gfa"
    #log:
    #    "regions/logs/r{region}.log.gz"
    shell:
        "python3 %s/haplotype.py {input.fasta} -o {params.output} -t {threads} -p {input.paf} 2>/dev/null" % (config["tools"]["paftest"])

rule prune_graph:
    input:
        gfa = "regions/gfas/r{region}.reducted.gfa"
    output:        
        gfa = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_min_length}.d{degree_max_size}.gfa"
    shell:
        """python3 %s/prune_tips.py {input.gfa} remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python3 %s/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python3 %s/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_degree3.py --max_size {wildcards.degree_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python3 %s/prune_tips.py remove --max_size {wildcards.tip_max_size} > {output.gfa}""" % (config["tools"]["paftest"], config["tools"]["paftest"],
                                                                                                    config["tools"]["paftest"], config["tools"]["paftest"],
                                                                                                    config["tools"]["paftest"], config["tools"]["paftest"],
                                                                                                    config["tools"]["paftest"], config["tools"]["paftest"])

#################
#Extract contigs#
#################

rule convert_lemon:
    input:
        gfa = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa"
    output:
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon",
        table = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.tbl"
    shell:
        "python3 %s/convert_to_lemon.py {input.gfa} {output.table} > {output.lemon}" % (config["tools"]["paftest"])

rule compute_path_cover:
    input:
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon"
    output:
        cover = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    shell:
        "%s {input.lemon} {output.cover}" % (config["tools"]["mc-mpc"])

rule cover_statistics:
    input:
        covers = expand("regions/gfas/pruned/r{region}.reducted.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.cover", region = good_regions),
        gfas = expand("regions/gfas/pruned/r{region}.reducted.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.gfa", region = good_regions)
    output:
        "regions/stats/stats.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.txt"
    run:
        assert len(input["covers"]) == len(input["gfas"])
        max_covers = []
        for i in input["covers"]:
            with open(i, "r") as cover:
                max_cover = -1
                for line in cover:
                    fields = line.strip().split()
                    if int(fields[1]) > max_cover:
                        max_cover = int(fields[1])
                max_covers.append(max_cover)

        component_numbers = []
        for j in input["gfas"]:
            with open(j, "r") as gfa:
                G = nx.Graph()
                for line in gfa:
                    if line[0] == 'L':
                        fields = line.split('\t')
                        node1 = fields[1]
                        node2 = fields[3]
                        G.add_edge(node1, node2)
                component_numbers.append(nx.number_connected_components(G))

        with open(output[0], "w") as output_file:
            print("#region\tmax_cover\tnum_components", file=output_file)
            for k in range(len(input["covers"])):
                print("%d\t%d\t%d" % (k+1, max_covers[k], component_numbers[k]), file=output_file)

rule ul_align_to_graph:
    input:
        graph = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        nano = "regions/fastas_nanopore/r{region}.fasta"
    output:
        aln = "regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json"
    log: "regions/logs/nanopore_alignment/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.log"
    threads: 5
    shell:
        "%s -g {input.graph} -f {input.nano} -t {threads} --seeds-mum-count 30000 -a {output.aln} --seeds-mxm-length 10 -b 35 1>{log} 2>&1" % (config["tools"]["graphaligner"])

rule extract_contigs_from_cover:
    input:
        gfa = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        json = "regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json",
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon",
        table = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.tbl",
        cover = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    output:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa",
        reads = "regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/reads.tsv"
    params:
        prefix = "r{region}"
    shell:
        "python3 %s/contigs_from_nanopore_cover.py --prefix {params.prefix} {input.gfa} {input.lemon} {input.table} {input.cover} --json {input.json} --paths {output.reads} > {output.contigs}" % (config["tools"]["paftest"])

rule polish_contigs:
    input:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa",
        reads = "regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/reads.tsv",
        contained_reads = "regions/gfas/r{region}.contained.reads",
        allreads = config["reads"]["pacbio"]
    output:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    params:
        wd = "regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/"
    threads: 4
    run:
        with open(input["contigs"], 'r') as contigs_file:
            content = contigs_file.read()
        if len(content) < 1:
            shell("touch {output.contigs}")
        else:
            #Extract contigs into separate files
            shell("samtools faidx {input.contigs}")
            fai_path = input["contigs"] + ".fai"
            with open(fai_path, 'r') as index_file:
                num = 0
                for line in index_file:
                    fields = line.strip().split()
                    name = fields[0]
                    num += 1
                    path = params["wd"] + "contig%d.fa" % (num)
                    shell("samtools faidx {input.contigs} {name} > {path}")
            num_haps = num
            
            #Store reads for each contig
            with open(input["reads"], 'r') as read_file:
                reads = defaultdict(list)
                for line in read_file:
                    fields = line.strip().split()
                    num = int(fields[0])
                    assert num <= num_haps
                    reads[num].append(fields[1])
            
            #Add contained reads
            with open(input["contained_reads"], 'r') as contained_file:
                for num in range(1, num_haps + 1):
                    for line in contained_file:
                        fields = line.strip().split()
                        if fields[1] in reads[num]:
                            reads[num].append(fields[0])
            
            #Polish each contig
            for num in range(1, num_haps + 1):
                reads_path = params["wd"] + "reads%d.txt" % (num)
                with open(reads_path, 'w') as read_file:
                    for r in reads[num]:
                        print(r, file=read_file)
                fastq_path = params["wd"] + "reads%d.fq" % (num)
                fasta_path = params["wd"] + "reads%d.fa" % (num)
                shell("samtools faidx -f -r {reads_path} {input.allreads} > {fastq_path}")
                shell("seqtk seq -A {fastq_path} > {fasta_path}")
                contig_path = params["wd"] + "contig%d.fa" % (num)
                bam_path = params["wd"] + "aligned%d.bam" % (num)
                shell("minimap2 -ax map-pb -t4 {contig_path} {fasta_path} | samtools view -b | samtools sort > {bam_path}")
                shell("samtools index {bam_path}")
                pileup_path = params["wd"] + "pileup%d.vcf" % (num)
                snps_path = params["wd"] + "snps%d.vcf.gz" % (num)
                indels_path = params["wd"] + "indels%d.vcf.gz" % (num)
                calls_path = params["wd"] + "calls%d.vcf.gz" % (num)
                shell("bcftools mpileup -Q0 -o20 -e10 -Ov -f {contig_path} {bam_path} > {pileup_path}")
                shell("bcftools view -m2 -M3 -v snps -i 'DP > 3 & (I16[2] + I16[3]) / (I16[0] + I16[1]) > 1' -Oz {pileup_path} > {snps_path}")
                shell("bcftools index {snps_path}")
                shell("bcftools view -m2 -M2 -v indels -i 'DP > 3 & (I16[2] + I16[3]) / (I16[0] + I16[1]) > 1 & IMF>0.5' -Oz {pileup_path} > {indels_path}")
                shell("bcftools index {indels_path}")
                shell("bcftools concat {snps_path} {indels_path} | bcftools sort -Oz > {calls_path}")
                shell("bcftools index {calls_path}")
                norm_path = params["wd"] + "calls%d.norm.bcf" % (num)
                consensus_path = params["wd"] + "consensus%d.fa" % (num)
                shell("bcftools norm -f {contig_path} {calls_path} -Ob -o {norm_path}")
                shell("bcftools index {norm_path}")
                shell("cat {contig_path} | bcftools consensus {norm_path} > {consensus_path}")
            
            #Concat polished contigs
            shell("cat %s > {output.contigs}" % (" ".join([(params["wd"] + "consensus%d.fa" % (num)) for num in range(1, num_haps + 1)])))

rule self_align_contigs:
    input:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    output:
        bam = "regions/contigs/alignments.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/r{region}.bam"
    threads: 4
    shell:
        "minimap2 -x asm20 -Y -a --eqx -t {threads} {input.contigs} {input.contigs} | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule compute_identity_table:
    input:
        bam = "regions/contigs/alignments.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/r{region}.bam"
    output:
        tbl = "regions/contigs/alignments.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/r{region}.tbl"
    shell:
        "python3 %s/samIdentity.py --header {input.bam} | awk '$1 != $5' > {output.tbl}" % (config["tools"]["paftest"])

rule separate_paralogs:
    input:
        tbl = "regions/contigs/alignments.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/r{region}.tbl",
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    output:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.diploid.fa"
    log:
        "regions/logs/merge_haplotypes/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.log"
    shell:
        "python3 %s/merge_haplotypes.py {input.contigs} {input.tbl} {wildcards.region} 2> {log} > {output.contigs}" % (config["tools"]["paftest"])

rule remove_duplicates_haploid:
    input:
        tbl = "regions/contigs/alignments.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/r{region}.tbl",
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    output:
        contigs = "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haploid.fa"
    log:
        "regions/logs/merge_haplotypes/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.log"
    shell:
        "python3 %s/merge_haplotypes.py {input.contigs} {input.tbl} {wildcards.region} 2> {log} > {output.contigs}" % (config["tools"]["paftest"])

#############
#Plot graphs#
#############

rule plot_bandage_raw:
    input:
        graph = "regions/gfas/r{region}.reducted.gfa"
    output:
        png = "regions/pngs/r{region}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then %s image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi" % (config["tools"]["bandage"])

rule plot_bandage_pruned:
    input:
        graph = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa"
    output:
        png = "regions/pngs/r{region}.pruned.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then %s image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi" % (config["tools"]["bandage"])

#############
#Map contigs#
#############


rule pool_contigs_diploid:
    input:
        expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.diploid.fa", i=good_regions)
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa"
    shell:
        "cat {input} > {output}"

rule pool_contigs_haploid:
    input:
        expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.haploid.fa", i=good_regions)
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haploid.fa"
    shell:
        "cat {input} > {output}"

rule split_to_diploid:
    input:
        fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa",
        fai = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa.fai"
    output:
        h0 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype0.fa",
        h1 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype1.fa",
        h2 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype2.fa"
    run:
        shell("samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype0\") {input.fa} > {output.h0}")
        shell("samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype1\") {input.fa} > {output.h1}")
        shell("samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype2\") {input.fa} > {output.h2}")

rule filter_haplotype0:
    input:
        fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa",
        fai = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa.fai"
    output:
        fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.diploid.fa",
    run:
        shell("samtools faidx -r <(cat {input.fai} | cut -f 1 | grep -v \"_haplotype0\") {input.fa} > {output.fa}")


rule map_contigs_to_assembly:
    input:
        fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        ref = config["genome"]
    output:
        "regions/eval/{parameters}/bams/polished.{ploidy}.to.assembly.bam"
    threads: 10
    shell:
        "minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
        {input.ref} {input.fa} | samtools view -F 260 -u - | samtools sort - > {output}"

rule map_contigs_to_bacs:
    input:
        asm = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        bacs = config["bacs"]["fasta"]
    output:
        bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.bacs.bam"
    threads: 10
    shell:
        "minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
        {input.bacs} {input.asm} | samtools view -F 2308 -u - | samtools sort - > {output.bam}"

rule map_bacs_to_contigs:
    input:
        asm = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        bacs = config["bacs"]["fasta"]
    output:
        bam = "regions/eval/{parameters}/bams/bacs.to.polished.{ploidy}.bam"
    threads: 10
    shell:
        "minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
        {input.asm} {input.bacs} | samtools view -F 2308 -u - | samtools sort - > {output.bam}"

rule index_bam:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    shell:
        "samtools index {input}"

rule map_contigs_to_hg38:
    input:
        asm = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        hg38 = config["hg38"]
    output:
        bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.hg38.bam"
    threads: 10
    shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
    -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
     {input.hg38} {input.asm} | samtools view -F 260 -u - | samtools sort -@ {threads} - > {output.bam}
"""

rule map_bacs_to_assembly:
    input:
        bacs = config["bacs"]["fasta"],
        ref = config["genome"]
    output:
        "regions/eval/bams/bacs.to.assembly.bam"
    threads: 10
    shell:
        "minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
        {input.ref} {input.bacs} | samtools view -F 260 -u - | samtools sort - > {output}"

rule map_bacs_to_hg38:
    input:
        hg38 = config["hg38"],
        sd = config["segdups"], 
        bacs = config["bacs"]["fasta"],
    output:
        bam = "regions/eval/bams/bacs.to.hg38.bam",
        names = "regions/eval/bams/bacs.in.sd.names",
    threads: 10
    shell:"""
minimap2 -t {threads} --secondary=no -a -Y -x asm20 \
    -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
     {input.hg38} {input.bacs} | samtools view -F 260 -u - | samtools sort -@ {threads} - > {output.bam}

bedtools intersect -a {output.bam} -b {input.sd} | samtools view - | awk '{{print $1}}' > {output.names}
"""

##########
#Evaluate#
##########

rule quast_to_assembly:
    input:
        fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        ref = config["genome"]
    output:
        "regions/eval/{parameters}/quast_to_assembly/{ploidy}/report.html"
    params:
        wd = "regions/eval/{parameters}/quast_to_assembly/{ploidy}"
    shell:
        "%s --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 95.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}" % (config["tools"]["quast"])

rule quast_to_bacs:
    input:
        fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        ref = config["bacs"]["fasta"]
    output:
        "regions/eval/{parameters}/quast_to_bacs/{ploidy}/report.html"
    params:
        wd = "regions/eval/{parameters}/quast_to_bacs/{ploidy}"
    shell:
        "%s --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 98.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}" % (config["tools"]["quast"])

rule quast_to_hg38:
    input:
        fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        ref = config["hg38"]
    output:
        "regions/eval/{parameters}/quast_to_hg38/{ploidy}/report.html"
    params:
        wd = "regions/eval/{parameters}/quast_to_hg38/{ploidy}"
    shell:
        "%s --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 95.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}" % (config["tools"]["quast"])

rule bam_to_bed:
    input:
        bam = "{name}.bam"
    output:
        bed = "{name}.bed"
    shell:
        "bedtools bamtobed -i {input.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}"

rule count_resolved_segdup_regions_assembly:
    input:
        bed = "regions/eval/{parameters}/bams/polished.{ploidy}.to.assembly.bed",
        regions = config["regions"]["original"]["sda"],
        index = config["genome"] + ".fai"
    output:
        inter = "regions/eval/{parameters}/resolved.{ploidy}/inter.bed",
        stats = "regions/eval/{parameters}/resolved.{ploidy}/stats.txt",
        status_list = "regions/eval/{parameters}/resolved.{ploidy}/list.tbl"
    params:
        min_mapq = 30,
        extra = 0
    run:
        shell("bedtools intersect -a {input.bed} -b {input.regions} -wa -wb > {output.inter}")
        #Read segdup regions
        with open(input["regions"]) as segdupFile:
            segdupLines = segdupFile.readlines()
        segdups = {}
        for line in segdupLines:
            vals = line.split()[0:3]
            start = int(vals[1]) 
            end = int(vals[2])
            segdups["_".join(vals[0:3])] = []
        #Read fasta index
        contig_lengths = dict()
        with open(input["index"], 'r') as fai_file:
            for line in fai_file:
                fields = line.strip().split()
                contig_lengths[fields[0]] = int(fields[1])
        with open(output["inter"]) as tabFile:
            for line in tabFile:
                vals = line.split()
                mapq = int(vals[4])
                if (mapq < params["min_mapq"]):
                    continue
                (ctgChrom, ctgStart, ctgEnd, ctgName) = (vals[0], int(vals[1]), int(vals[2]), vals[3])
                (segChrom, segStart, segEnd) = (vals[5], int(vals[6]), int(vals[7]))
                pref = segStart - ctgStart
                suff = ctgEnd - segEnd
                segdup = "_".join(vals[5:8])
                segdups[segdup].append((pref, suff, ctgStart < 10, contig_lengths[ctgChrom] - ctgEnd < 10))
        counter = Counter()
        stats = open(output["stats"], 'w')
        status_list = open(output["status_list"], 'w')
        for line in segdupLines:
            vals = line.split()
            segdup = "_".join(vals[0:3])
            resolving_contigs = 0
            for pref, suff, atStart, atEnd in segdups[segdup]:
                if (atStart or pref > params["extra"]) and (atEnd or suff > params["extra"]):
                    resolving_contigs += 1
            counter[resolving_contigs] += 1
            # write all entries
            vals.append( "|".join(["%d,%d" % (pref, suff) for pref, suff, atStart, atEnd in segdups[segdup]]))
            line = "\t".join(vals)
            print(line, file=status_list)
        for i in range(max(counter.keys()) + 1):
            print("%d: %d" % (i, counter[i]), file=stats)
        print("Total: %d" % (sum(counter.values())), file=stats)
        stats.close()
        status_list.close()

ruleorder: count_resolved_segdup_regions_assembly > bam_to_bed


rule bacs_to_contigs_tbl:
    input:
        bam = "regions/eval/{parameters}/bams/bacs.to.polished.{ploidy}.bam",
        names = "regions/eval/bams/bacs.in.sd.names"
    output:
        tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
    shell:
        "python3 %s/samIdentity.py --header --mask {input.names} {input.bam} > {output.tbl}" % (config["tools"]["paftest"])

rule contigs_to_hg38_tbl:
    input:
        bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.hg38.bam"
    output:
        tbl = "regions/eval/{parameters}/tables/polished.{ploidy}.to.hg38.tbl"
    shell:
        "python3 %s/samIdentity.py --header <(samtools view -h -F 2308 {input.bam}) > {output.tbl}" % (config["tools"]["paftest"])

rule contigs_to_bacs_tbl:
    input:
        bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.bacs.bam"
    output:
        tbl = "regions/eval/{parameters}/tables/polished.{ploidy}.to.bacs.tbl"
    shell:
        "python3 %s/samIdentity.py --header {input.bam} > {output.tbl}" % (config["tools"]["paftest"])

rule make_qv_sum:
    input:
        bac_tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
    output:
        qv_sum = "regions/eval/{parameters}/tables/qv_sum.{ploidy}.txt"
    run:
        pd.options.mode.use_inf_as_na = True
        out = ""
        sys.stderr.write(input["bac_tbl"] + "\n")
        df = pd.read_csv(input["bac_tbl"], sep="\t")
        val_all = 1 - df["perID_by_all"]/100
        val_matches = 1 - df["perID_by_matches"]/100
        df["qv_all"] = -10 * np.log10( val_all )
        df["qv_matches"] = -10 * np.log10( val_matches )
        for mask in [True, False]:
            tmp = df[df["mask"] == mask]
            mean_all = tmp["perID_by_all"].mean()
            mean_all_qv = -10 * np.log10(1 - mean_all/100)
            mean_matches = tmp["perID_by_matches"].mean()
            mean_matches_qv = -10 * np.log10(1 - mean_matches/100)
            perfect_all = tmp["qv_all"].isna()
            perfect_matches = tmp["qv_matches"].isna()
            out += "Segdup BACs? = {}\n".format(mask)
            out += "\tMean identity\tQV of mean identity\n"
            out += "All\t{}\t{}\n".format(mean_all, mean_all_qv)
            out += "Matches\t{}\t{}\n\n".format(mean_matches, mean_matches_qv)
            out += "All (SDA method)\nPerfect\t{}\n{}\n\n".format(sum(perfect_all), tmp.qv_all.describe())            
            out += "Matches (SDA method)\nPerfect\t{}\n{}\n\n".format(sum(perfect_matches), tmp.qv_matches.describe())            
        open(output["qv_sum"], "w+").write(out)

rule count_misassemblies:
    input:
        bac_tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
    output:
        confirmed = "regions/eval/{parameters}/misassemblies.{ploidy}/confirmed.txt",
        misassembled = "regions/eval/{parameters}/misassemblies.{ploidy}/misassembled.txt",
        unclear = "regions/eval/{parameters}/misassemblies.{ploidy}/unclear.txt"
    params:
        threshold = 5000
    run:
        shell("awk -v t=\"{params.threshold}\" '{{ if (($6 < t) && (($8-$7) < t)) {{ print $5 }} }}' {input.bac_tbl} | uniq | sort -k 1,1 > {output.confirmed}")
        shell("join -t$'\\t' -v 1 -1 5 -2 1 <(tail -n+2 {input.bac_tbl} | sort -k 5,5) {output.confirmed} | awk -v t=\"{params.threshold}\" '{{ if (($3>t && $6>t) || ($5-$4>t && $8-$7>t)) {{ print $1 }} }}' | uniq > {output.misassembled}")
        shell("join -t$'\\t' -v 1 -1 5 -2 1 <(tail -n+2 {input.bac_tbl} | sort -k 5,5) <(cat {output.confirmed} {output.misassembled} | sort -k1,1) | cut -f 1 | uniq > {output.unclear}")

#########
#SV Call#
#########

rule svim_call_contigs:
    input:
        bam =  "regions/eval/{parameters}/bams/polished.{ploidy}.to.hg38.bam",
        bai =  "regions/eval/{parameters}/bams/polished.{ploidy}.to.hg38.bam.bai",
        genome = config["hg38"]
    output:
        "regions/svim/{parameters}/contigs.haploid/variants.vcf"
    params:
        working_dir = "regions/svim/{parameters}/contigs.haploid/"
    shell:
        "%s haploid {params.working_dir} {input.bam} {input.genome}" % (config["tools"]["svim-asm"])

rule svim_call_bacs:
    input:
        bam =  "regions/eval/bams/bacs.to.hg38.bam",
        bai =  "regions/eval/bams/bacs.to.hg38.bam.bai",
        genome = config["hg38"]
    output:
        "regions/svim/{parameters}/bacs.haploid/variants.vcf"
    params:
        working_dir = "regions/svim/{parameters}/bacs.haploid/"
    shell:
        "%s haploid {params.working_dir} {input.bam} {input.genome}" % (config["tools"]["svim-asm"])


