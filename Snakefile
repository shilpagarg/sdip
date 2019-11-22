configfile: "config.yaml"
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/
big_regions = [17, 32, 96, 101, 102, 108, 162, 166, 216, 217, 218, 267, 268, 275, 278, 285, 296, 306, 307, 433, 467, 470, 471, 503, 504]
good_regions = [elem for elem in list(range(1, 505)) if not elem in big_regions]

rule all:
    input:
        expand("regions/pngs/r{i}.{subset}.png", i=good_regions, subset=["primary", "primary.full"]),
        expand("regions/pngs/r{i}.{subset}.pruned.notips{tip_max_size}.nobubbles{bubble_max_size}.png", i=good_regions,
                                                                                                        subset=["primary", "primary.full"],
                                                                                                        tip_max_size=[5],
                                                                                                        bubble_max_size=[5]),
        expand("regions/contigs/pooled.{subset}.notips{tip_max_size}.nobubbles{bubble_max_size}.fa", subset=["primary", "primary.full"],
                                                                                                     tip_max_size=[5],
                                                                                                     bubble_max_size=[5])

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

wildcard_constraints:
    region="\d+",
    tip_max_size="\d+",
    bubble_max_size="\d+"


rule minimap:
    input:
        fastq = get_samples,
        reference = config["genome"]
    output:
        temp("alignment/{sample}.bam")
    threads: 10
    shell:
        "minimap2 -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -t {threads} --eqx --secondary=yes -N 100 -p 0.9 {input.reference} {input.fastq} | samtools view -b | samtools sort > {output}"

rule pool_samples:
    input:
        expand("alignment/{sample}.bam", sample=config["samples"])
    output:
        "alignment/pooled.bam"
    threads: 4
    shell:
        "samtools merge -@ {threads} -r {output} {input}"

rule samtools_index:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    threads: 4
    shell:
        "samtools index -@ {threads} {input}"

rule recruit_primary_reads:
    input:
        regions = config["regions"],
        bam = "alignment/pooled.bam",
        index = "alignment/pooled.bam.bai"
    output:
        temp("regions/bams/r{region,\d+}.primary.bam")
    params:
        primary_mapq_threshold = 1 
    threads: 5
    shell:
        "chr=`awk '{{if (NR=={wildcards.region}) {{ print $1 }} }}' {input.regions}`; start=`awk '{{ if (NR=={wildcards.region}) {{ print $2 }} }}' {input.regions}`; end=`awk '{{ if (NR=={wildcards.region}) {{ print $3 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{chr}}:${{start}}-${{end}} > {output}"

rule recruit_secondary_reads:
    input:
        regions = config["regions"],
        bam = "alignment/pooled.bam",
        index = "alignment/pooled.bam.bai"
    output:
        temp("regions/bams/r{region,\d+}.secondary.bam")
    threads: 5
    shell:
        "chr=`awk '{{if (NR=={wildcards.region}) {{ print $1 }} }}' {input.regions}`; start=`awk '{{ if (NR=={wildcards.region}) {{ print $2 }} }}' {input.regions}`; end=`awk '{{ if (NR=={wildcards.region}) {{ print $3 }} }}' {input.regions}`; samtools view -b -@ {threads} -f 256 {input.bam} ${{chr}}:${{start}}-${{end}} > {output}"

rule get_read_names:
    input:
        "regions/bams/r{region}.{subset}.bam"
    output:
        "regions/reads/r{region,\d+}.{subset}.reads"
    shell:
        "samtools view {input} | cut -d$'\t' -f1 | sort | uniq > {output}"

rule grab_more_from_overlaps:
    input:
        "regions/reads/r{region}.{subset}.reads"
    output:
        "regions/reads/r{region}.{subset}.full.reads"
    threads: 10
    shell:
        "python3 ../WHdenovo/paftest/getOv.py --num_chunks 50 --num_threads 10 ../../../tmp/hg002-asm-r3-pg0.1.5.3/0-seqdb/seq_dataset.idx ../../../tmp/hg002-asm-r3-pg0.1.5.3/3-asm/split/ {input} | sort | uniq > {output}"

rule merge_primary_and_secondary:
    input:
        primary = "regions/reads/r{region}.primary.reads",
        secondary = "regions/reads/r{region}.secondary.reads"
    output:
        "regions/reads/r{region,\d+}.merged.reads"
    shell:
        "cat {input} > {output}"

rule get_fastq:
    input:
        allreads = config["reads"],
        names = "regions/reads/r{region}.{subset}.reads"
    output:
        temp("regions/fastas/r{region}.{subset}.fastq")
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ]; then seqtk subseq {input.allreads} {input.names} > {output}; else echo '' > {output}; fi"

rule convert_fastq_to_fasta:
    input:
        "regions/fastas/r{region}.{subset}.fastq"
    output:
        "regions/fastas/r{region}.{subset}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule all_vs_all_alignment:
    input:
        fasta = "regions/fastas/r{region}.{subset}.fasta"
    output:
        paf = "regions/pafs/r{region}.{subset}.fasta.paf"
    threads: 10
    shell:
        "minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t {threads} {input.fasta} {input.fasta} | sort -k8n -k6 > {output.paf}"

rule filter_paf_for_long_indels:
    input:
        paf = "regions/pafs/r{region}.{subset}.fasta.paf"
    output:
        paf = "regions/pafs/r{region}.{subset}.fasta.filtered.paf"
    shell:
        "python3 ../WHdenovo/paftest/filter_indels_in_paf.py {input.paf} --max_indel 50 > {output.paf}"

rule generate_graph_use_paf:
    input:
        fasta = "regions/fastas/r{region}.{subset}.fasta",
        paf = "regions/pafs/r{region}.{subset}.fasta.filtered.paf"
    output:
        gfa = "regions/gfas/r{region}.{subset}.reducted.gfa"
    threads: 3
    params:
        output = "regions/gfas/r{region}.{subset}.gfa"
    #log:
    #    "regions/logs/r{region}.log.gz"
    shell:
        "python3 ../WHdenovo/paftest/haplotype.py {input.fasta} -o {params.output} -t {threads} -p {input.paf} 2>/dev/null"

rule prune_graph:
    input:
        gfa = "regions/gfas/r{region}.{subset}.reducted.gfa"
    output:        
        gfa = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.gfa"
    shell:
        """python3 ../WHdenovo/paftest/remove_tips.py {input.gfa} remove --max_size {wildcards.tip_max_size} | \
           python3 ../WHdenovo/paftest/find_ultrabubbles.py remove --max_size {wildcards.bubble_max_size} | \
           python3 ../WHdenovo/paftest/remove_tips.py remove --max_size {wildcards.tip_max_size} | \
           python3 ../WHdenovo/paftest/find_ultrabubbles.py remove --max_size {wildcards.bubble_max_size} | \
           python3 ../WHdenovo/paftest/remove_tips.py remove --max_size {wildcards.tip_max_size} > {output.gfa}"""

rule convert_lemon:
    input:
        gfa = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.gfa"
    output:
        lemon = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.lemon"
    shell:
        "python3 ../WHdenovo/paftest/convert_to_lemon.py {input.gfa} > {output.lemon}"

rule compute_path_cover:
    input:
        lemon = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.lemon"
    output:
        cover = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.cover"
    shell:
        "../bin/mc-mpc {input.lemon} {output.cover}"

rule extract_contigs:
    input:
        gfa = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.gfa",
        lemon = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.lemon",
        cover = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.cover"
    output:
        "regions/contigs/r{region}.{subset}.notips{tip_max_size}.nobubbles{bubble_max_size}.fa"
    shell:
        "python3 ../WHdenovo/paftest/iteratePathsFromCover.py {input.gfa} {input.lemon} {input.cover} > {output}"

rule plot_bandage_raw:
    input:
        graph = "regions/gfas/r{region}.{subset}.reducted.gfa"
    output:
        png = "regions/pngs/r{region}.{subset}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then ../bin/Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"

rule plot_bandage_pruned:
    input:
        graph = "regions/gfas/pruned/r{region}.{subset}.reducted.notips{tip_max_size}.nobubbles{bubble_max_size}.gfa"
    output:
        png = "regions/pngs/r{region}.{subset}.pruned.notips{tip_max_size}.nobubbles{bubble_max_size}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then ../bin/Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"

rule pool_contigs:
    input:
        expand("regions/contigs/r{i}.{{subset}}.notips{{tip_max_size}}.nobubbles{{bubble_max_size}}.fa", i=good_regions)
    output:
        "regions/contigs/pooled.{subset}.notips{tip_max_size}.nobubbles{bubble_max_size}.fa"    
    shell:
        "cat {input} > {output}"

rule map_contigs:
    input:
        fa = "regions/contigs/pooled.fa",
        ref = "/project/heller_d-data/genomes/hg38/hg38.mmi"
    output:
        "regions/contigs/pooled.sorted.bam"
    shell:
        "minimap2 -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -t 10 --eqx --secondary=yes -N 100 -p 0.9 \
        {input.ref} {input.fa} | samtools view -b | samtools sort > {output}"

rule index_bam:
    input:
        "regions/contigs/pooled.sorted.bam"
    output:
        "regions/contigs/pooled.sorted.bam.bai"
    shell:
        "samtools index {input}"
