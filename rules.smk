def get_samples(wildcards):
    return config["samples"][wildcards.sample]

wildcard_constraints:
    region="\d+",
    tip_max_size="\d+",
    bubble_max_size="\d+"

#############################
#Find similar segdup regions#
#############################

rule fetch_segdup_sequences:
    input:
        regions = config["regions"],
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

rule bamtobed:
    input:
        "regions/segdups/r{region}.bam"
    output:
        "regions/segdups/r{region}.tsv"
    shell:
        "bedtools bamtobed -i {input} | awk 'OFS=\"\\t\" {{ print {wildcards.region}, $0 }}' > {output}"

rule concat_regions:
    input:
        expand("regions/segdups/r{region}.tsv", region=regions)
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
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

rule get_nanopore_fastq:
    input:
        allreads = config["reads"]["nanopore"],
        names = "regions/reads_nanopore/r{region}.reads"
    output:
        temp("regions/fastas_nanopore/r{region}.fastq")
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ]; then samtools faidx -f -r {input.names} {input.allreads} > {output}; else echo '' > {output}; fi"

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
        "python3 /project/pacbiosv/code/WHdenovo/paftest/filter_indels_in_paf.py {input.paf} --max_indel 50 > {output.paf}"

rule generate_graph_use_paf:
    input:
        fasta = "regions/fastas_pacbio/r{region}.fasta",
        paf = "regions/pafs/r{region}.fasta.filtered.paf"
    output:
        gfa = "regions/gfas/r{region}.reducted.gfa"
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
        gfa = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa"
    shell:
        """python3 %s/prune_tips.py {input.gfa} remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --max_size {wildcards.bubble_max_size} | \
           python3 %s/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --max_size {wildcards.bubble_max_size} | \
           python3 %s/prune_tips.py remove --max_size {wildcards.tip_max_size} | \
           python3 %s/prune_degree3.py --max_size {wildcards.degree_max_size} | \
           python3 %s/prune_ultrabubbles.py remove --max_size {wildcards.bubble_max_size} | \
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
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon"
    shell:
        "python3 %s/convert_to_lemon.py {input.gfa} > {output.lemon}" % (config["tools"]["paftest"])

rule compute_path_cover:
    input:
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon"
    output:
        cover = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    shell:
        "%s {input.lemon} {output.cover}" % (config["tools"]["mc-mpc"])

rule cover_statistics:
    input:
        expand("regions/gfas/pruned/r{region}.reducted.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.cover", region = regions)
    output:
        "regions/stats/cover.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.txt"
    run:
        max_covers = []
        for i in input:
            input_file = open(i, "r")
            max_cover = -1
            for line in input_file:
                fields = line.strip().split()
                if int(fields[1]) > max_cover:
                    max_cover = int(fields[1])
            input_file.close()
            max_covers.append(max_cover)
        output_file = open(output[0], "w")
        for i, m in enumerate(max_covers):
            print("%d\t%d" % (i+1, m), file=output_file)
        output_file.close()

rule ul_align_to_graph:
    input:
        graph = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        nano = "regions/fastas_nanopore/r{region}.fasta"
    output:
        aln = "regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json"
    log: "regions/logs/nanopore_alignment/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.log"
    threads: 5
    shell:
        "GraphAligner -g {input.graph} -f {input.nano} -t {threads} --seeds-mum-count 30000 -a {output.aln} --seeds-mxm-length 10 -b 35 1>{log} 2>&1"

rule extract_contigs_from_cover:
    input:
        gfa = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.gfa",
        #json = "regions/jsons/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.json",
        lemon = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.lemon",
        cover = "regions/gfas/pruned/r{region}.reducted.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.cover"
    output:
        "regions/contigs/r{region}.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa"
    params:
        prefix = "r{region}"
    shell:
        "python3 %s/contigs_from_nanopore_cover.py --prefix {params.prefix} {input.gfa} {input.lemon} {input.cover} > {output}" % (config["tools"]["paftest"])


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

rule pool_contigs:
    input:
        expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.fa", i=regions)
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa"    
    shell:
        "cat {input} > {output}"

rule map_contigs:
    input:
        fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa",
        ref = config["genome"]
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.sorted.bam"
    threads: 10
    shell:
        "minimap2 -ax asm5 -t {threads} --eqx \
        {input.ref} {input.fa} | samtools view -b | samtools sort > {output}"

rule index_bam:
    input:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.sorted.bam"
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.sorted.bam.bai"
    shell:
        "samtools index {input}"