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


########################
# Align BACs and Reads #
########################

rule align_pacbio_reads:
    input:
        ref = config["hg38"],
        fasta = config["pacbio"]
    output:
        bam = "pipeline/alignment_pacbio/pooled.sorted.bam"
    threads: 20
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -t {threads} --eqx --MD --secondary=no {input.ref} {input.fasta} | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule align_nanopore_reads:
    input:
        ref = config["hg38"],
        fasta = config["nanopore"]
    output:
        bam = "pipeline/alignment_nanopore/pooled.sorted.bam"
    threads: 20
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -aL -z 600,200 -x map-ont -t {threads} --eqx --secondary=yes -N 100 -p 0.97 {input.ref} {input.fasta} | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule align_long_nanopore_reads:
    input:
        ref = config["hg38"],
        fasta = config["nanopore"]
    output:
        bam = "pipeline/alignment_nanopore/pooled.long.sorted.bam"
    threads: 20
    conda:
        "../envs/minimap.yaml"
    shell:
        "/project/pacbiosv/bin/bioawk/bioawk -c fastx  '{{ if(length($seq) > 100000){{ print \"@\"$name; print $seq; print \"+\"$comment; print $qual }} }}' {input.fasta} | \
         minimap2 -aL -z 600,200 -x map-ont -t {threads} --eqx --MD --secondary=yes -N 100 -p 0.97 {input.ref} - | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule map_bacs_to_hg38:
    input:
        hg38 = config["hg38"], 
        bacs = config["bacs"],
    output:
        bam = "pipeline/alignment_bacs/bacs.sorted.bam"
    threads: 10
    conda:
        "../envs/minimap.yaml"
    shell:"""
minimap2 -t {threads} --secondary=no -a -Y -x asm20 \
    -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
     {input.hg38} {input.bacs} | samtools view -F 260 -u - | samtools sort -@ {threads} - > {output.bam}
"""

rule bam_to_bed:
    input:
        bam = "{name}.bam"
    output:
        bed = "{name}.bed"
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed}"

rule bedtools_merge:
    input:
        bed = "{name}.bed"
    output:
        bed = "{name}.merged.bed"
    shell:
        "bedtools merge -i {input.bed} > {output.bed}"

rule index_bam:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

############
# Call SVs #
############

#Run Sniffles
rule sniffles_call_pacbio_reads:
    input:
        bam = "pipeline/alignment_pacbio/pooled.sorted.bam",
        bai = "pipeline/alignment_pacbio/pooled.sorted.bam.bai"
    output:
        "pipeline/calls/sniffles/raw_{min_support}.vcf"
    conda:
        "../envs/sniffles.yaml"
    shell:
        "sniffles --mapped_reads {input.bam} --min_length 40 --min_support {wildcards.min_support} --vcf {output} --threads 1"

#see https://github.com/spiralgenetics/truvari/issues/43
rule fix_sniffles_pacbio_reads:
    input:
        "pipeline/calls/sniffles/raw_{min_support}.vcf"
    output:
        "pipeline/calls/sniffles/min_{min_support,[0-9]+}.vcf"
    shell:
        "sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}"

#Run svim-asm on Bacs
rule svim_call_bacs:
    input:
        bam =  "pipeline/alignment_bacs/bacs.sorted.bam",
        bai =  "pipeline/alignment_bacs/bacs.sorted.bam.bai",
        genome = config["hg38"]
    output:
        "pipeline/calls/bacs/variants.vcf"
    params:
        working_dir = "pipeline/calls/bacs/"
    shell:
        "/home/heller_d/bin/anaconda3/bin/svim-asm haploid --min_sv_size 40 --min_mapq 0 {params.working_dir} {input.bam} {input.genome}"

rule sniffles_call_nanopore_reads:
    input:
        bam = "pipeline/alignment_nanopore/pooled.long.sorted.bam",
        bai = "pipeline/alignment_nanopore/pooled.long.sorted.bam.bai"
    output:
        "pipeline/calls/ul_ont/raw_{min_support}.vcf"
    conda:
        "../envs/sniffles.yaml"
    shell:
        "sniffles --mapped_reads {input.bam} --min_length 40 --min_support {wildcards.min_support} --vcf {output} --threads 1"

#see https://github.com/spiralgenetics/truvari/issues/43
rule fix_sniffles_ont_reads:
    input:
        "pipeline/calls/ul_ont/raw_{min_support}.vcf"
    output:
        "pipeline/calls/ul_ont/min_{min_support,[0-9]+}.vcf"
    shell:
        "sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}"

###############
#Recruit reads#
###############

rule recruit_pacbio_reads:
    input:
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        bam = "pipeline/alignment_pacbio/pooled.sorted.bam",
        index = "pipeline/alignment_pacbio/pooled.sorted.bam.bai"
    output:
        "pipeline/regions/bams_pacbio/r{region}.bam"
    params:
        primary_mapq_threshold = 1 
    threads: 5
    conda:
        "../envs/samtools.yaml"
    shell:
        "reg=`awk '{{if (NR=={wildcards.region}) {{ print $1\":\"$2\"-\"$3 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{reg}} > {output}"


rule recruit_nanopore_reads:
    input:
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        bam = "pipeline/alignment_nanopore/pooled.sorted.bam",
        index = "pipeline/alignment_nanopore/pooled.sorted.bam.bai"
    output:
        "pipeline/regions/bams_nanopore/r{region}.bam"
    params:
        primary_mapq_threshold = 30
    threads: 5
    conda:
        "../envs/samtools.yaml"
    shell:
        "reg=`awk '{{if (NR=={wildcards.region}) {{ print $1\":\"$2\"-\"$3 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{reg}} > {output}"

rule get_read_names:
    input:
        "pipeline/regions/bams_{technology}/r{region}.bam"
    output:
        "pipeline/regions/reads_{technology}/r{region}.reads"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {input} | cut -d$'\t' -f1 | sort | uniq > {output}"

rule get_pacbio_fastq:
    input:
        allreads = config["pacbio"],
        names = "pipeline/regions/reads_pacbio/r{region}.reads"
    output:
        temp("pipeline/regions/fastas_pacbio/r{region}.fastq")
    conda:
        "../envs/samtools.yaml"
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ] && [ $LINES -gt 0 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

rule get_nanopore_fastq:
    input:
        allreads = config["nanopore"],
        names = "pipeline/regions/reads_nanopore/r{region}.reads"
    output:
        temp("pipeline/regions/fastas_nanopore/r{region}.fastq")
    conda:
        "../envs/samtools.yaml"
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ] && [ $LINES -gt 0 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

rule convert_fastq_to_fasta:
    input:
        "{name}.fastq"
    output:
        "{name}.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -A {input} > {output}"

rule all_vs_all_alignment:
    input:
        fasta = "pipeline/regions/fastas_pacbio/r{region}.fasta"
    output:
        paf = "pipeline/regions/pafs/r{region}.fasta.paf"
    threads: 10
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t {threads} {input.fasta} {input.fasta} | sort -k8n -k6 > {output.paf}"

rule filter_paf_for_long_indels:
    input:
        paf = "pipeline/regions/pafs/r{region}.fasta.paf"
    output:
        paf = "pipeline/regions/pafs/r{region}.fasta.filtered.paf"
    shell:
        "python workflow/scripts/filter_indels_in_paf.py {input.paf} --max_indel 50 > {output.paf}"

rule generate_graph_use_paf:
    input:
        fasta = "pipeline/regions/fastas_pacbio/r{region}.fasta",
        paf = "pipeline/regions/pafs/r{region}.fasta.filtered.paf"
    output:
        gfa = "pipeline/regions/gfas/r{region}.reducted.gfa",
        contained = "pipeline/regions/gfas/r{region}.contained.reads"
    threads: 3
    params:
        output = "pipeline/regions/gfas/r{region}.gfa"
    #log:
    #    "regions/logs/r{region}.log.gz"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/haplotype.py {input.fasta} -o {params.output} -t {threads} -p {input.paf} 2>/dev/null"

rule prune_graph:
    input:
        gfa = "pipeline/regions/gfas/r{region}.reducted.gfa"
    output:        
        gfa = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_min_length}.d{degree_max_size}.gfa"
    shell:
        """python workflow/scripts/prune_tips.py {input.gfa} remove --max_size {wildcards.tip_max_size} | \
           python workflow/scripts/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python workflow/scripts/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python workflow/scripts/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python workflow/scripts/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python workflow/scripts/prune_degree3.py --max_size {wildcards.degree_max_size} | \
           python workflow/scripts/prune_ultrabubbles.py remove --min_length {wildcards.bubble_min_length} | \
           python workflow/scripts/prune_tips.py remove --max_size {wildcards.tip_max_size} > {output.gfa}"""

#################
#Extract contigs#
#################

rule convert_lemon:
    input:
        gfa = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa"
    output:
        lemon = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon",
        table = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.tbl"
    shell:
        "python workflow/scripts/convert_to_lemon.py {input.gfa} {output.table} > {output.lemon}"

rule compute_path_cover:
    input:
        lemon = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon"
    output:
        cover = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    shell:
        "workflow/bin/mc-mpc {input.lemon} {output.cover}"

rule cover_statistics:
    input:
        covers = expand("pipeline/regions/gfas/pruned/r{region}.reducted.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.cover", region = good_regions),
        gfas = expand("pipeline/regions/gfas/pruned/r{region}.reducted.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.gfa", region = good_regions)
    output:
        "pipeline/regions/stats/stats.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.txt"
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
        graph = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        nano = "pipeline/regions/fastas_nanopore/r{region}.fasta"
    output:
        aln = "pipeline/regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json"
    log: "pipeline/regions/logs/nanopore_alignment/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.log"
    threads: 5
    conda:
        "../envs/graphaligner.yaml"
    shell:
        "GraphAligner -g {input.graph} -f {input.nano} -t {threads} --seeds-mum-count 30000 -a {output.aln} --seeds-mxm-length 10 -b 35 1>{log} 2>&1"

rule extract_contigs_from_cover:
    input:
        gfa = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        json = "pipeline/regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json",
        lemon = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon",
        table = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.tbl",
        cover = "pipeline/regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    output:
        contigs = "pipeline/regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa",
        reads = "pipeline/regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/reads.tsv"
    params:
        prefix = "r{region}"
    shell:
        "python workflow/scripts/contigs_from_nanopore_cover.py --prefix {params.prefix} {input.gfa} {input.lemon} {input.table} {input.cover} --json {input.json} --paths {output.reads} > {output.contigs}"

rule polish_contigs:
    input:
        contigs = "pipeline/regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa",
        reads = "pipeline/regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/reads.tsv",
        contained_reads = "pipeline/regions/gfas/r{region}.contained.reads",
        sequences = config["pacbio"]
    output:
        contigs = "pipeline/regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    params:
        wd = "pipeline/regions/polishing/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/"
    threads: 4
    conda:
        "../envs/polishing.yaml"
    script:
        "../scripts/polishing.py"

rule self_align_contigs:
    input:
        contigs = "pipeline/regions/contigs/r{region}.{parameters}.polished.fa"
    output:
        bam = "pipeline/regions/contigs/alignments.{parameters}/r{region}.bam"
    threads: 4
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -x asm20 -Y -a --eqx -t {threads} {input.contigs} {input.contigs} | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule compute_identity_table:
    input:
        bam = "pipeline/regions/contigs/alignments.{parameters}/r{region}.bam"
    output:
        tbl = "pipeline/regions/contigs/alignments.{parameters}/r{region}.tbl"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/samIdentity.py --header {input.bam} | awk '$1 != $5' > {output.tbl}"

rule separate_paralogs:
    input:
        tbl = "pipeline/regions/contigs/alignments.{parameters}/r{region}.tbl",
        contigs = "pipeline/regions/contigs/r{region}.{parameters}.polished.fa",
        fai = "pipeline/regions/contigs/r{region}.{parameters}.polished.fa.fai"
    output:
        contigs = "pipeline/regions/contigs/r{region}.{parameters}.polished.diploid.fa"
    log:
        "pipeline/regions/logs/merge_haplotypes/r{region}.{parameters}.log"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/merge_haplotypes.py --max_similarity 99.95 {input.contigs} {input.tbl} {wildcards.region} 2> {log} > {output.contigs}"

rule faidx:
    input:
        "{name}.fa"
    output:
        "{name}.fa.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

rule remove_duplicates_haploid:
    input:
        tbl = "pipeline/regions/contigs/alignments.{parameters}/r{region}.tbl",
        contigs = "pipeline/regions/contigs/r{region}.{parameters}.polished.fa",
        fai = "pipeline/regions/contigs/r{region}.{parameters}.polished.fa.fai"
    output:
        contigs = "pipeline/regions/contigs/r{region}.{parameters}.polished.haploid.fa"
    log:
        "pipeline/regions/logs/merge_haplotypes/r{region}.{parameters}.log"
    shell:
        "python workflow/scripts/merge_haplotypes.py --haploid --max_similarity 99.95 {input.contigs} {input.tbl} {wildcards.region} 2> {log} > {output.contigs}"

rule pool_contigs_haploid:
    input:
        expand("pipeline/regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.haploid.fa", i=good_regions)
    output:
        "pipeline/regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haploid.fa"
    shell:
        "cat {input} > {output}"

rule map_contigs_to_hg38:
    input:
        asm = "pipeline/regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
        hg38 = config["hg38"]
    output:
        bam = "pipeline/alignment_contigs/contigs.{parameters}.{ploidy}.sorted.bam"
    threads: 10
    conda:
        "../envs/minimap.yaml"
    shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
    -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
     {input.hg38} {input.asm} | samtools view -F 2308 -u - | samtools sort -@ {threads} - > {output.bam}
"""

#Run svim-asm on contigs
rule svim_call_sdip:
    input:
        bam =  "pipeline/alignment_contigs/contigs.{parameters}.{ploidy}.sorted.bam",
        bai =  "pipeline/alignment_contigs/contigs.{parameters}.{ploidy}.sorted.bam.bai",
        genome = config["hg38"]
    output:
        "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/variants.vcf"
    params:
        working_dir = "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/"
    shell:
        "/home/heller_d/bin/anaconda3/bin/svim-asm haploid --min_sv_size 40 --min_mapq {wildcards.mapq} {params.working_dir} {input.bam} {input.genome}"


############
# Eval SVs #
############

rule filter_on_regions_and_type:
    input:
        calls = "pipeline/calls/{tool}/{name}.vcf"
    output:
        "pipeline/calls/{tool}/{name}.filtered.vcf.gz"
    shell:
        "bcftools view -i \"SVTYPE = 'INS' | SVTYPE = 'DEL'\" {input.calls} | bcftools sort -Oz > {output}"

rule deduplicate:
    input:
        "pipeline/calls/{tool}/{name}.filtered.vcf.gz"
    output:
        "pipeline/calls/{tool}/{name}.filtered.dedup.vcf.gz"
    shell:
        "bcftools norm -d all -Oz {input} > {output}"

rule tabix_calls:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "tabix -p vcf {input}"

rule truvari_to_bacs_sniffles:
    input:
        reference = config["hg38"],
        truth = "pipeline/calls/bacs/variants.filtered.dedup.vcf.gz",
        truth_index = "pipeline/calls/bacs/variants.filtered.dedup.vcf.gz.tbi",
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        calls = "pipeline/calls/sniffles/min_{min_support}.filtered.vcf.gz",
        calls_index = "pipeline/calls/sniffles/min_{min_support}.filtered.vcf.gz.tbi"
    output:
        "pipeline/truvari_to_bacs/sniffles/{min_support}/summary.txt",
        "pipeline/truvari_to_bacs/sniffles/{min_support}/tp-call.vcf",
        "pipeline/truvari_to_bacs/sniffles/{min_support}/tp-base.vcf",
        "pipeline/truvari_to_bacs/sniffles/{min_support}/fp.vcf",
        "pipeline/truvari_to_bacs/sniffles/{min_support}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_to_bacs/sniffles/{min_support}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --passonly \
                -r 1000 \
                -p 0.00 \
                -c {input.calls}"""

rule truvari_to_ont_sniffles:
    input:
        reference = config["hg38"],
        truth = "pipeline/calls/ul_ont/min_10.filtered.dedup.vcf.gz",
        truth_index = "pipeline/calls/ul_ont/min_10.filtered.dedup.vcf.gz.tbi",
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        calls = "pipeline/calls/sniffles/min_{min_support}.filtered.vcf.gz",
        calls_index = "pipeline/calls/sniffles/min_{min_support}.filtered.vcf.gz.tbi"
    output:
        "pipeline/truvari_to_ont/sniffles/{min_support}/summary.txt",
        "pipeline/truvari_to_ont/sniffles/{min_support}/tp-call.vcf",
        "pipeline/truvari_to_ont/sniffles/{min_support}/tp-base.vcf",
        "pipeline/truvari_to_ont/sniffles/{min_support}/fp.vcf",
        "pipeline/truvari_to_ont/sniffles/{min_support}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_to_ont/sniffles/{min_support}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --passonly \
                -r 1000 \
                -p 0.00 \
                -c {input.calls}"""

rule truvari_to_bacs_sdip:
    input:
        reference = config["hg38"],
        truth = "pipeline/calls/bacs/variants.filtered.dedup.vcf.gz",
        truth_index = "pipeline/calls/bacs/variants.filtered.dedup.vcf.gz.tbi",
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        calls = "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/variants.filtered.vcf.gz",
        calls_index = "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/variants.filtered.vcf.gz.tbi"
    output:
        "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}/summary.txt",
        "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}/tp-call.vcf",
        "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}/tp-base.vcf",
        "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}/fp.vcf",
        "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_to_bacs/sdip/{parameters}.{ploidy}.{mapq}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --passonly \
                -r 1000 \
                -p 0.00 \
                -c {input.calls}"""

rule truvari_to_ont_sdip:
    input:
        reference = config["hg38"],
        truth = "pipeline/calls/ul_ont/min_10.filtered.dedup.vcf.gz",
        truth_index = "pipeline/calls/ul_ont/min_10.filtered.dedup.vcf.gz.tbi",
        regions = "pipeline/alignment_bacs/bacs.sorted.merged.bed",
        calls = "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/variants.filtered.vcf.gz",
        calls_index = "pipeline/calls/contigs.{parameters}.{ploidy}.{mapq}/variants.filtered.vcf.gz.tbi"
    output:
        "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}/summary.txt",
        "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}/tp-call.vcf",
        "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}/tp-base.vcf",
        "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}/fp.vcf",
        "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_to_ont/sdip/{parameters}.{ploidy}.{mapq}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --passonly \
                -r 1000 \
                -p 0.00 \
                -c {input.calls}"""

rule reformat_truvari_results:
    input:
        "pipeline/{truvari}/{tool}/{parameters}/summary.txt"
    output:
        "pipeline/{truvari}/{tool}/{parameters}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | grep -v 'gt' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.tool}\", \"{wildcards.parameters}\", $1, $2 }}' > {output}"

rule cat_truvari_results:
    input:
        sniffles = expand("pipeline/{{truvari}}/sniffles/{threshold}/pr_rec.txt", 
                          threshold = [1, 4, 7, 10]),
        sdip = expand("pipeline/{{truvari}}/sdip/t5.b5000.d2.haploid.{mapq}/pr_rec.txt",
                       mapq = [0, 30, 50, 60])
    output:
        all = "pipeline/eval/{truvari}/all_results.txt"
    threads: 1
    shell:
        "cat {input.sniffles} {input.sdip} > {output.all}"

rule plot_pr_tools:
    input:
        "pipeline/eval/{truvari}/all_results.txt"
    output:
        "pipeline/eval/{truvari}/results.tools.png"
    threads: 1
    log:
        "pipeline/logs/rplot.{truvari}.tools.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools.R {input} {output} > {log}"

#############################
#Find similar segdup regions#
#############################

# rule fetch_segdup_sequences:
#     input:
#         regions = config["regions"],
#         reference = config["base"]
#     output:
#         "regions/segdups/r{region}.fa"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "reg=`awk '{{if (NR=={wildcards.region}) {{ print $1\":\"$2\"-\"$3 }} }}' {input.regions}`; samtools faidx {input.reference} ${{reg}} > {output}"

# rule map_segdup_sequences:
#     input:
#         fasta = "regions/segdups/r{region}.fa",
#         reference = config["base"]
#     output:
#         "regions/segdups/r{region}.bam"
#     threads: 10
#     conda:
#         "../envs/minimap.yaml"
#     shell:
#         "minimap2 -ax asm10 -t {threads} --secondary=yes -N 100 -p 0.80 {input.reference} {input.fasta} | samtools view -b | samtools sort > {output}"

# rule map_segdup_sequences_to_hg38:
#     input:
#         fasta = "regions/segdups/r{region}.fa",
#         reference = config["hg38"]
#     output:
#         "regions/segdups/aln_to_hg38/r{region}.bam"
#     threads: 2
#     conda:
#         "../envs/minimap.yaml"
#     shell:
#         "minimap2 -ax asm10 -t {threads} --secondary=yes -N 100 -p 0.80 {input.reference} {input.fasta} | samtools view -b | samtools sort > {output}"

# rule bam_to_tsv:
#     input:
#         "regions/segdups/r{region}.bam"
#     output:
#         "regions/segdups/r{region}.tsv"
#     conda:
#         "../envs/bedtools.yaml"
#     shell:
#         "bedtools bamtobed -i {input} | awk 'OFS=\"\\t\" {{ print {wildcards.region}, $0 }}' > {output}"



# rule faidx2:
#     input:
#         "{name}.fasta"
#     output:
#         "{name}.fasta.fai"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools faidx {input}"

# rule concat_regions:
#     input:
#         expand("regions/segdups/r{region}.tsv", region=all_segdup_regions)
#     output:
#         "segdups.similar.tsv"
#     run:
#         o = open(output[0], 'w')
#         index = 0
#         for p in input:
#             with open(p, 'r') as f:
#                 lines = f.readlines()
#             if len(lines) < 5:
#                 index += 1
#                 for l in lines:
#                     fields = l.strip().split()                    
#                     print("%d\t%s\t%s\t%s" % (index, fields[1], fields[2], fields[3]), file = o)
#         o.close()


#############
#Plot graphs#
#############

# rule plot_bandage_raw:
#     input:
#         graph = "regions/gfas/r{region}.reducted.gfa"
#     output:
#         png = "regions/pngs/r{region}.png"
#     conda:
#         "../envs/bandage.yaml"
#     shell:
#         "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"

# rule plot_bandage_pruned:
#     input:
#         graph = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa"
#     output:
#         png = "regions/pngs/r{region}.pruned.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.png"
#     conda:
#         "../envs/bandage.yaml"
#     shell:
#         "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"

#############
#Map contigs#
#############

# rule pool_contigs_diploid:
#     input:
#         expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.diploid.fa", i=good_regions)
#     output:
#         "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa"
#     shell:
#         "cat {input} > {output}"

# rule pool_contigs_haploid:
#     input:
#         expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.haploid.fa", i=good_regions)
#     output:
#         "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haploid.fa"
#     shell:
#         "cat {input} > {output}"

# rule split_to_diploid:
#     input:
#         fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa",
#         fai = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa.fai"
#     output:
#         h0 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype0.fa",
#         h1 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype1.fa",
#         h2 = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.haplotype2.fa"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype0\") {input.fa} > {output.h0} &&
#            samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype1\") {input.fa} > {output.h1} &&
#            samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype2\") {input.fa} > {output.h2}"""

# rule filter_haplotype0:
#     input:
#         fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa",
#         fai = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.grouped.fa.fai"
#     output:
#         fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.diploid.fa",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools faidx -r <(cat {input.fai} | cut -f 1 | grep -v \"_haplotype0\") {input.fa} > {output.fa}"

# rule map_contigs_to_bacs:
#     input:
#         asm = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
#         bacs = config["bacs"]
#     output:
#         bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.bacs.bam"
#     threads: 10
#     conda:
#         "../envs/minimap.yaml"
#     shell:
#         "minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 \
#         {input.bacs} {input.asm} | samtools view -F 2308 -u - | samtools sort - > {output.bam}"

# rule map_bacs_to_contigs:
#     input:
#         asm = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
#         bacs = config["bacs"]
#     output:
#         bam = "regions/eval/{parameters}/bams/bacs.to.polished.{ploidy}.bam"
#     threads: 10
#     conda:
#         "../envs/minimap.yaml"
#     shell:
#         "minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
#         {input.asm} {input.bacs} | samtools view -F 2308 -u - | samtools sort - > {output.bam}"


##########
#Evaluate#
##########

# rule quast_to_assembly:
#     input:
#         fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
#         ref = config["base"]
#     output:
#         "regions/eval/{parameters}/quast_to_assembly/{ploidy}/report.html"
#     params:
#         wd = "regions/eval/{parameters}/quast_to_assembly/{ploidy}"
#     conda:
#         "../envs/quast.yaml"
#     shell:
#         "quast.py --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 95.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}"

# rule quast_to_bacs:
#     input:
#         fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
#         ref = config["bacs"]
#     output:
#         "regions/eval/{parameters}/quast_to_bacs/{ploidy}/report.html"
#     params:
#         wd = "regions/eval/{parameters}/quast_to_bacs/{ploidy}"
#     conda:
#         "../envs/quast.yaml"
#     shell:
#         "quast.py --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 98.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}"

# rule quast_to_hg38:
#     input:
#         fa = "regions/contigs/pooled.{parameters}.polished.{ploidy}.fa",
#         ref = config["hg38"]
#     output:
#         "regions/eval/{parameters}/quast_to_hg38/{ploidy}/report.html"
#     params:
#         wd = "regions/eval/{parameters}/quast_to_hg38/{ploidy}"
#     conda:
#         "../envs/quast.yaml"
#     shell:
#         "quast.py --fragmented --min-contig 20000 --min-alignment 5000 --min-identity 95.0 --unaligned-part-size 2000 --no-icarus -o {params.wd} -r {input.ref} {input.fa}"

# rule bam_to_bed:
#     input:
#         bam = "{name}.bam"
#     output:
#         bed = "{name}.bed"
#     conda:
#         "../envs/bedtools.yaml"
#     shell:
#         "bedtools bamtobed -i {input.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}"

# rule overlap_alignments_with_regions:
#     input:
#         bed = "regions/eval/{parameters}/bams/polished.{ploidy}.to.assembly.bed",
#         regions = config["regions"]
#     output:
#         inter = "regions/eval/{parameters}/resolved.{ploidy}/inter.bed"
#     shell:
#         "bedtools intersect -a {input.bed} -b {input.regions} -wa -wb > {output.inter}"

# rule count_resolved_segdup_regions_assembly:
#     input:
#         regions = config["regions"],
#         inter = "regions/eval/{parameters}/resolved.{ploidy}/inter.bed",
#         index = config["base"] + ".fai"
#     output:
#         stats = "regions/eval/{parameters}/resolved.{ploidy}/stats.txt",
#         status_list = "regions/eval/{parameters}/resolved.{ploidy}/list.tbl",
#         resolved = "regions/eval/{parameters}/resolved.{ploidy}/resolved.txt"
#     params:
#         min_mapq = 30,
#         extra = 0
#     run:
#         #Read segdup regions
#         with open(input["regions"]) as segdupFile:
#             segdupLines = segdupFile.readlines()
#         segdups = {}
#         for line in segdupLines:
#             vals = line.split()[0:3]
#             start = int(vals[1]) 
#             end = int(vals[2])
#             segdups["_".join(vals[0:3])] = []
#         #Read fasta index
#         contig_lengths = dict()
#         with open(input["index"], 'r') as fai_file:
#             for line in fai_file:
#                 fields = line.strip().split()
#                 contig_lengths[fields[0]] = int(fields[1])
#         with open(input["inter"]) as tabFile:
#             for line in tabFile:
#                 vals = line.split()
#                 mapq = int(vals[4])
#                 if (mapq < params["min_mapq"]):
#                     continue
#                 (ctgChrom, ctgStart, ctgEnd, ctgName) = (vals[0], int(vals[1]), int(vals[2]), vals[3])
#                 (segChrom, segStart, segEnd) = (vals[5], int(vals[6]), int(vals[7]))
#                 pref = segStart - ctgStart
#                 suff = ctgEnd - segEnd
#                 segdup = "_".join(vals[5:8])
#                 segdups[segdup].append((pref, suff, ctgStart < 10, contig_lengths[ctgChrom] - ctgEnd < 10))
#         counter = Counter()
#         resolved_regions = set()
#         stats = open(output["stats"], 'w')
#         status_list = open(output["status_list"], 'w')
#         resolved_file = open(output["resolved"], 'w')
#         for line in segdupLines:
#             vals = line.split()
#             segdup = "_".join(vals[0:3])
#             resolving_contigs = 0
#             for pref, suff, atStart, atEnd in segdups[segdup]:
#                 if (atStart or pref > params["extra"]) and (atEnd or suff > params["extra"]):
#                     resolving_contigs += 1
#                     resolved_regions.add((vals[0], vals[1], vals[2]))
#             counter[resolving_contigs] += 1
#             # write all entries
#             vals.append( "|".join(["%d,%d" % (pref, suff) for pref, suff, atStart, atEnd in segdups[segdup]]))
#             line = "\t".join(vals)
#             print(line, file=status_list)
#         for i in range(max(counter.keys()) + 1):
#             print("%d: %d" % (i, counter[i]), file=stats)
#         for region in resolved_regions:
#             print("\t".join(region), file=resolved_file)
#         print("Total: %d" % (sum(counter.values())), file=stats)
#         stats.close()
#         status_list.close()

# rule bacs_to_contigs_tbl:
#     input:
#         bam = "regions/eval/{parameters}/bams/bacs.to.polished.{ploidy}.bam"
#     output:
#         tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
#     conda:
#         "../envs/python.yaml"
#     shell:
#         "python workflow/scripts/samIdentity.py --header {input.bam} > {output.tbl}"

# rule contigs_to_hg38_tbl:
#     input:
#         bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.hg38.bam"
#     output:
#         tbl = "regions/eval/{parameters}/tables/polished.{ploidy}.to.hg38.tbl"
#     conda:
#         "../envs/python.yaml"
#     shell:
#         "python workflow/scripts/samIdentity.py --header {input.bam} > {output.tbl}"

# rule contigs_to_bacs_tbl:
#     input:
#         bam = "regions/eval/{parameters}/bams/polished.{ploidy}.to.bacs.bam"
#     output:
#         tbl = "regions/eval/{parameters}/tables/polished.{ploidy}.to.bacs.tbl"
#     conda:
#         "../envs/python.yaml"
#     shell:
#         "python workflow/scripts/samIdentity.py --header {input.bam} > {output.tbl}"

# rule make_qv_sum:
#     input:
#         bac_tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
#     output:
#         qv_sum = "regions/eval/{parameters}/tables/qv_sum.{ploidy}.txt"
#     run:
#         pd.options.mode.use_inf_as_na = True
#         out = ""
#         sys.stderr.write(input["bac_tbl"] + "\n")
#         df = pd.read_csv(input["bac_tbl"], sep="\t")
#         val_all = 1 - df["perID_by_all"]/100
#         val_matches = 1 - df["perID_by_matches"]/100
#         df["qv_all"] = -10 * np.log10( val_all )
#         df["qv_matches"] = -10 * np.log10( val_matches )
#         mean_all = df["perID_by_all"].mean()
#         mean_all_qv = -10 * np.log10(1 - mean_all/100)
#         mean_matches = df["perID_by_matches"].mean()
#         mean_matches_qv = -10 * np.log10(1 - mean_matches/100)
#         perfect_all = df["qv_all"].isna()
#         perfect_matches = df["qv_matches"].isna()
#         out += "\tMean identity\tQV of mean identity\n"
#         out += "All\t{}\t{}\n".format(mean_all, mean_all_qv)
#         out += "Matches\t{}\t{}\n\n".format(mean_matches, mean_matches_qv)
#         out += "All (SDA method)\nPerfect\t{}\n{}\n\n".format(sum(perfect_all), df.qv_all.describe())            
#         out += "Matches (SDA method)\nPerfect\t{}\n{}\n\n".format(sum(perfect_matches), df.qv_matches.describe())            
#         open(output["qv_sum"], "w+").write(out)

# rule count_misassemblies:
#     input:
#         bac_tbl = "regions/eval/{parameters}/tables/bacs.to.polished.{ploidy}.tbl"
#     output:
#         confirmed = "regions/eval/{parameters}/misassemblies.{ploidy}/confirmed.txt",
#         misassembled = "regions/eval/{parameters}/misassemblies.{ploidy}/misassembled.txt",
#         unclear = "regions/eval/{parameters}/misassemblies.{ploidy}/unclear.txt"
#     params:
#         threshold = 5000
#     run:
#         shell("awk -v t=\"{params.threshold}\" '{{ if (($6 < t) && (($8-$7) < t)) {{ print $5 }} }}' {input.bac_tbl} | uniq | sort -k 1,1 > {output.confirmed}")
#         shell("join -t$'\\t' -v 1 -1 5 -2 1 <(tail -n+2 {input.bac_tbl} | sort -k 5,5) {output.confirmed} | awk -v t=\"{params.threshold}\" '{{ if (($3>t && $6>t) || ($5-$4>t && $8-$7>t)) {{ print $1 }} }}' | uniq > {output.misassembled}")
#         shell("join -t$'\\t' -v 1 -1 5 -2 1 <(tail -n+2 {input.bac_tbl} | sort -k 5,5) <(cat {output.confirmed} {output.misassembled} | sort -k1,1) | cut -f 1 | uniq > {output.unclear}")

#########
#SV Call#
#########

# rule svim_call_diploid:
#     input:
#         bam1 =  "regions/eval/{parameters}/bams/polished.haplotype1.to.hg38.bam",
#         bam2 =  "regions/eval/{parameters}/bams/polished.haplotype2.to.hg38.bam",
#         bai1 =  "regions/eval/{parameters}/bams/polished.haplotype1.to.hg38.bam.bai",
#         bai2 =  "regions/eval/{parameters}/bams/polished.haplotype2.to.hg38.bam.bai",
#         genome = config["hg38"]
#     output:
#         "regions/svim/{parameters}/contigs.diploid/variants.vcf"
#     params:
#         working_dir = "regions/svim/{parameters}/contigs.diploid/"
#     conda:
#         "../envs/svimasm.yaml"
#     shell:
#         "svim-asm diploid {params.working_dir} {input.bam1} {input.bam2} {input.genome}"

# rule svim_call_haploid:
#     input:
#         bam =  "regions/eval/{parameters}/bams/polished.haploid.to.hg38.bam",
#         bai =  "regions/eval/{parameters}/bams/polished.haploid.to.hg38.bam.bai",
#         genome = config["hg38"]
#     output:
#         "regions/svim/{parameters}/contigs.haploid/variants.vcf"
#     params:
#         working_dir = "regions/svim/{parameters}/contigs.haploid/"
#     conda:
#         "../envs/svimasm.yaml"
#     shell:
#         "svim-asm haploid --min_mapq 0 {params.working_dir} {input.bam} {input.genome}"

# rule normalize_calls:
#     input:
#         "{name}.vcf"
#     output:
#         "{name}.norm.vcf"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "bcftools norm -d all {input} > {output}"

# rule compress_calls:
#     input:
#         "{name}.vcf"
#     output:
#         gz = "{name}.vcf.gz",
#         tbi = "{name}.vcf.gz.tbi"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "bgzip -c {input} > {output.gz} && tabix -p vcf {output.gz}"
