# Collection of Snakemake rules imported in the main Snakefile.

def get_pacbio_reads(wildcards):
    return config["reads"][wildcards.sample]["pacbio"]

def get_pacbio_index(wildcards):
    return config["reads"][wildcards.sample]["pacbio"] + ".fai"

def get_nanopore_reads(wildcards):
    return config["reads"][wildcards.sample]["nanopore"]

rule faidx:
    input:
        "{name}.fa"
    output:
        "{name}.fa.fai"
    shell:
        "samtools faidx {input}"

rule faidx2:
    input:
        "{name}.fasta"
    output:
        "{name}.fasta.fai"
    shell:
        "samtools faidx {input}"

#########################
#Construct overlap graph#
#########################

rule all_vs_all_alignment:
    input:
        fasta = get_pacbio_reads
    output:
        paf = "{sample}/{sample}.paf"
    threads: 10
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t {threads} {input.fasta} {input.fasta} | sort -k8n -k6 > {output.paf}"

rule filter_paf_for_long_indels:
    input:
        paf = "{sample}/{sample}.paf"
    output:
        paf = "{sample}/{sample}.filtered.paf"
    shell:
        "python3 workflow/scripts/filter_indels_in_paf.py {input.paf} --max_indel 50 > {output.paf}"

rule generate_graph_use_paf:
    input:
        fasta = get_pacbio_reads,
        paf = "{sample}/{sample}.filtered.paf"
    output:
        gfa = "{sample}/{sample}.reducted.gfa",
        contained = "{sample}/{sample}.contained.reads"
    threads: 3
    params:
        output = "{sample}/{sample}.gfa"
    #log:
    #    "regions/logs/r{region}.log.gz"
    shell:
        "python3 workflow/scripts/haplotype.py {input.fasta} -o {params.output} -t {threads} -p {input.paf} 2>/dev/null"

rule prune_graph:
    input:
        gfa = "{sample}/{sample}.reducted.gfa"
    output:        
        gfa = "{sample}/{sample}.reducted.pruned.gfa"
    params:
        tip_max_size = config["params"]["tip_max_size"],
        bubble_min_length = config["params"]["bubble_min_length"],
        degree_max_size = config["params"]["degree_max_size"]
    shell:
        """python3 workflow/scripts/prune_tips.py {input.gfa} remove --max_size {params.tip_max_size} | \
           python3 workflow/scripts/prune_ultrabubbles.py remove --min_length {params.bubble_min_length} | \
           python3 workflow/scripts/prune_tips.py remove --max_size {params.tip_max_size} | \
           python3 workflow/scripts/prune_ultrabubbles.py remove --min_length {params.bubble_min_length} | \
           python3 workflow/scripts/prune_tips.py remove --max_size {params.tip_max_size} | \
           python3 workflow/scripts/prune_degree3.py --max_size {params.degree_max_size} | \
           python3 workflow/scripts/prune_ultrabubbles.py remove --min_length {params.bubble_min_length} | \
           python3 workflow/scripts/prune_tips.py remove --max_size {params.tip_max_size} > {output.gfa}"""


############################
#Extract and polish contigs#
############################

rule convert_lemon:
    input:
        gfa = "{sample}/{sample}.reducted.pruned.gfa"
    output:
        lemon = "{sample}/{sample}.reducted.pruned.lemon",
        table = "{sample}/{sample}.reducted.pruned.tbl"
    shell:
        "python3 workflow/scripts/convert_to_lemon.py {input.gfa} {output.table} > {output.lemon}"

rule compute_path_cover:
    input:
        lemon = "{sample}/{sample}.reducted.pruned.lemon"
    output:
        cover = "{sample}/{sample}.reducted.pruned.cover"
    shell:
        "/project/heller_d-data/shilpa/scripts/MC-MPC/mc-mpc-solver/mc-mpc {input.lemon} {output.cover}"

rule ul_align_to_graph:
    input:
        graph = "{sample}/{sample}.reducted.pruned.gfa",
        nano = get_nanopore_reads
    output:
        aln = "{sample}/{sample}.reducted.pruned.json"
    log: "logs/{sample}/nanopore_alignment/{sample}.log"
    threads: 5
    conda:
        "../envs/graphaligner.yaml"
    shell:
        "GraphAligner -g {input.graph} -f {input.nano} -t {threads} --seeds-mum-count 30000 -a {output.aln} --seeds-mxm-length 10 -b 35 1>{log} 2>&1"

rule extract_contigs_from_cover:
    input:
        gfa = "{sample}/{sample}.reducted.pruned.gfa",
        json = "{sample}/{sample}.reducted.pruned.json",
        lemon = "{sample}/{sample}.reducted.pruned.lemon",
        table = "{sample}/{sample}.reducted.pruned.tbl",
        cover = "{sample}/{sample}.reducted.pruned.cover"
    output:
        contigs = "{sample}/{sample}.fa",
        reads = "{sample}/{sample}.reads.tsv"
    params:
        prefix = "{sample}"
    shell:
        "python3 workflow/scripts/contigs_from_nanopore_cover.py --prefix {params.prefix} {input.gfa} {input.lemon} {input.table} {input.cover} --json {input.json} --paths {output.reads} > {output.contigs}"

rule polish_contigs:
    input:
        contigs = "{sample}/{sample}.fa",
        reads = "{sample}/{sample}.reads.tsv",
        contained_reads = "{sample}/{sample}.contained.reads",
        sequences = get_pacbio_reads,
        sequence_index = get_pacbio_index
    output:
        contigs = "{sample}/{sample}.polished.fa"
    params:
        wd = "{sample}/polishing/"
    threads: 4
    conda:
        "../envs/polishing.yaml"
    script:
        "../scripts/polishing.py"

##################
#Match Haplotypes#
##################

rule self_align_contigs:
    input:
        contigs = "{sample}/{sample}.polished.fa"
    output:
        bam = "{sample}/{sample}.polished.selfaln.bam"
    threads: 4
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -x asm20 -Y -a --eqx -t {threads} {input.contigs} {input.contigs} | samtools view -F 4 -u - | samtools sort - > {output.bam}"

rule compute_identity_table:
    input:
        bam = "{sample}/{sample}.polished.selfaln.bam"
    output:
        tbl = "{sample}/{sample}.polished.selfaln.tbl"
    shell:
        "python3 workflow/scripts/samIdentity.py --header {input.bam} | awk '$1 != $5' > {output.tbl}" 

rule separate_paralogs:
    input:
        tbl = "{sample}/{sample}.polished.selfaln.tbl",
        contigs = "{sample}/{sample}.polished.fa",
        fai = "{sample}/{sample}.polished.fa.fai"
    output:
        contigs = "{sample}/{sample}.polished.grouped.fa"
    log:
        "logs/{sample}/merge_haplotypes/{sample}.log"
    shell:
        "python3 workflow/scripts/merge_haplotypes.py --max_similarity 99.95 {input.contigs} {input.tbl} {wildcards.sample} 2> {log} > {output.contigs}"

rule remove_duplicates_haploid:
    input:
        tbl = "{sample}/{sample}.polished.selfaln.tbl",
        contigs = "{sample}/{sample}.polished.fa",
        fai = "{sample}/{sample}.polished.fa.fai"
    output:
        contigs = "{sample}/{sample}.polished.haploid.fa"
    log:
        "logs/{sample}/merge_haplotypes/{sample}.log"
    shell:
        "python3 workflow/scripts/merge_haplotypes.py --haploid --max_similarity 99.95 {input.contigs} {input.tbl} {wildcards.region} 2> {log} > {output.contigs}"

rule split_to_haplotypes:
    input:
        fa = "{sample}/{sample}.polished.grouped.fa",
        fai = "{sample}/{sample}.polished.grouped.fa.fai"
    output:
        h0 = "{sample}/{sample}.polished.haplotype0.fa",
        h1 = "{sample}/{sample}.polished.haplotype1.fa",
        h2 = "{sample}/{sample}.polished.haplotype2.fa",
        dip = "{sample}/{sample}.polished.diploid.fa"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype0\") {input.fa} > {output.h0}
        samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype1\") {input.fa} > {output.h1}
        samtools faidx -r <(cat {input.fai} | cut -f 1 | grep \"_haplotype2\") {input.fa} > {output.h2}
        samtools faidx -r <(cat {input.fai} | cut -f 1 | grep -v \"_haplotype0\") {input.fa} > {output.dip}
        """


##############
#SV Detection#
##############

rule map_contigs_to_hg38:
    input:
        asm = "{sample}/{sample}.polished.{ploidy}.fa",
        ref = config["reference"]
    output:
        bam = "{sample}/{sample}.polished.{ploidy}.to.ref.bam"
    threads: 10
    conda:
        "../envs/minimap.yaml"
    shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
    -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
     {input.ref} {input.asm} | samtools view -F 260 -u - | samtools sort -@ {threads} - > {output.bam}
"""

rule index_bam:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule svim_call_diploid:
    input:
        bam1 =  "{sample}/{sample}.polished.haplotype1.to.ref.bam",
        bam2 =  "{sample}/{sample}.polished.haplotype2.to.ref.bam",
        bai1 =  "{sample}/{sample}.polished.haplotype1.to.ref.bam.bai",
        bai2 =  "{sample}/{sample}.polished.haplotype2.to.ref.bam.bai",
        ref = config["reference"]
    output:
        "{sample}/sv_calls_diploid/variants.vcf"
    params:
        working_dir = "{sample}/sv_calls_diploid/"
    shell:
        "/home/heller_d/.local/bin/svim-asm diploid {params.working_dir} {input.bam1} {input.bam2} {input.ref}"

rule svim_call_haploid:
    input:
        bam =  "{sample}/{sample}.polished.haploid.to.ref.bam",
        bai =  "{sample}/{sample}.polished.haploid.to.ref.bam.bai",
        ref = config["reference"]
    output:
        "{sample}/sv_calls_haploid/variants.vcf"
    params:
        working_dir = "{sample}/sv_calls_haploid/"
    shell:
        "/home/heller_d/.local/bin/svim-asm haploid --min_mapq 0 {params.working_dir} {input.bam} {input.ref}"
