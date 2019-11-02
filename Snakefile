configfile: "config.yaml"

big_regions = [17, 32, 96, 101, 102, 108, 162, 166, 216, 217, 218, 267, 268, 275, 278, 285, 296, 306, 307, 433, 467, 470, 471, 503, 504]
good_regions = [elem for elem in list(range(1, 505)) if not elem in big_regions]

rule all:
    input:
        expand("regions/pngs/r{i}.png", i=good_regions),
        expand("regions/pngs/r{i}_pruned_nanopore.min{min_cov}.png", i=good_regions, min_cov=[1, 2, 4, 6, 10]),
        expand("regions/pngs/r{i}_pruned_nanopore_pacbio.min{min_cov}.png", i=good_regions, min_cov=[1, 2, 4, 6, 10])
        #expand("regions/pngs/r{i}_pruned2.min1.max5.png", i=good_regions),
        #expand("regions/gfas/r{i}.reducted.gfa", i=good_regions)
        #"regions/contigs/pooled.sorted.bam"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]

wildcard_constraints:
    region="\d+"


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
    threads: 30
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
        primary_mapq_threshold = 30
    threads: 20
    shell:
        "chr=`awk '{{if (NR=={wildcards.region}) {{ print $1 }} }}' {input.regions}`; start=`awk '{{ if (NR=={wildcards.region}) {{ print $2 }} }}' {input.regions}`; end=`awk '{{ if (NR=={wildcards.region}) {{ print $3 }} }}' {input.regions}`; samtools view -b -@ {threads} -F 256 -q {params.primary_mapq_threshold} {input.bam} ${{chr}}:${{start}}-${{end}} > {output}"

rule recruit_secondary_reads:
    input:
        regions = config["regions"],
        bam = "alignment/pooled.bam",
        index = "alignment/pooled.bam.bai"
    output:
        temp("regions/bams/r{region,\d+}.secondary.bam")
    threads: 20
    shell:
        "chr=`awk '{{if (NR=={wildcards.region}) {{ print $1 }} }}' {input.regions}`; start=`awk '{{ if (NR=={wildcards.region}) {{ print $2 }} }}' {input.regions}`; end=`awk '{{ if (NR=={wildcards.region}) {{ print $3 }} }}' {input.regions}`; samtools view -b -@ {threads} -f 256 {input.bam} ${{chr}}:${{start}}-${{end}} > {output}"

rule merge_primary_and_secondary:
    input:
        primary = "regions/bams/r{region}.primary.bam",
        secondary = "regions/bams/r{region}.secondary.bam"
    output:
        "regions/bams/r{region,\d+}.bam"
    shell:
        "samtools merge {output} {input}"

rule get_read_names:
    input:
        "regions/bams/r{region}.bam"
    output:
        "regions/reads/r{region}.reads"
    shell:
        "samtools view {input} | cut -d$'\t' -f1 | sort | uniq > {output}"

rule get_fastq:
    input:
        allreads = config["reads"],
        names = "regions/reads/r{region}.reads"
    output:
        temp("regions/fastas/r{region}.fastq")
    shell:
        "LINES=$(wc -l < {input.names}) ; if [ $LINES -lt 20000 ]; then /project/pacbiosv/bin/seqtk/seqtk subseq {input.allreads} {input.names} > {output}; else echo '' > {output}; fi"

rule convert_fastq_to_fasta:
    input:
        "regions/fastas/r{region}.fastq"
    output:
        "regions/fastas/r{region}.fasta"
    shell:
        "/project/pacbiosv/bin/seqtk/seqtk seq -A {input} > {output}"

# rule generate_graph:
#     input:
#         fasta = "regions/fastas/r{region}.fasta"
#     output:
#         gfa = "regions/gfas/r{region}.reducted.gfa"
#     threads: 10
#     params:
#         output = "regions/gfas/r{region}.gfa"
#     log:
#         "regions/logs/r{region}.log.gz"
#     shell:
#         "python3 /project/pacbiosv/code/WHdenovo/paftest/haplotype.py {input.fasta} -o {params.output} -t {threads} 2>&1 | gzip > {log}"

rule all_vs_all_alignment:
    input:
        fasta = "regions/fastas/r{region}.fasta"
    output:
        paf = "regions/pafs/r{region}.fasta.paf"
    threads: 10
    shell:
        "minimap2 -c -x asm20 -DP --no-long-join --cs -n500 -t {threads} {input.fasta} {input.fasta} | sort -k8n -k6 > {output.paf}"

# rule generate_graph:
#     input:
#         fasta = "regions/fastas/r{region}.fasta"
#     output:
#         gfa = "regions/gfas/r{region}.reducted.gfa",
#         paf = "regions/pafs/r{region}.fasta.paf"
#     threads: 10
#     params:
#         output = "regions/gfas/r{region}.gfa"
#     log:
#         "regions/logs/r{region}.log.gz"
#     shell:
#         "cd regions/pafs && python3 /project/pacbiosv/code/WHdenovo/paftest/haplotype.py /project/heller_d-data/shilpa/reads_vs_wholegenome/{input.fasta} -o /project/heller_d-data/shilpa/reads_vs_wholegenome/{params.output} -t {threads} 2>&1 | gzip > /project/heller_d-data/shilpa/reads_vs_wholegenome/{log}"

rule generate_graph_use_paf:
    input:
        fasta = "regions/fastas/r{region}.fasta",
        paf = "regions/pafs/r{region}.fasta.paf"
    output:
        gfa = "regions/gfas/r{region}.reducted.gfa"
    threads: 3
    params:
        output = "regions/gfas/r{region}.gfa"
    #log:
    #    "regions/logs/r{region}.log.gz"
    shell:
        "python3 /project/pacbiosv/code/WHdenovo/paftest/haplotype.py {input.fasta} -o {params.output} -t {threads} -p {input.paf} 2>/dev/null"

# rule remove_edges_covering_entire_nodes:
#     input:
#         gfa = "regions/gfas/r{region}.reducted.gfa"
#     output:
#         gfa = "regions/gfas/r{region}.reducted.nofullmatches.gfa"
#     shell:
#         "awk 'OFS=\"\t\" {{ match($6, \"([0-9]+)M\", m); if($1 == \"S\") {{ print $0 }}; if(m[1] != $7 && m[1] != $8) {{ print $0 }} }}' {input.gfa} > {output.gfa}"

rule align_nanopore_reads:
    input:
        gfa = "regions/gfas/r{region}.reducted.gfa",
        fasta = "../nanopore_vs_hg38/regions/fastas/r{region}.fasta"
    output:
        json = "regions/json_nanopore/r{region}.json"
    shell:
        "/project/pacbiosv/bin/GraphAligner/bin/GraphAligner --seeds-minimizer-chunksize 500 \
         --seeds-minimizer-length 15 --seeds-minimizer-count 20 -g {input.gfa} -f {input.fasta} -a {output.json}"

rule align_pacbio_reads:
    input:
        gfa = "regions/gfas/r{region}.reducted.gfa",
        fasta = "regions/fastas/r{region}.fasta"
    output:
        json = "regions/json_pacbio/r{region}.json"
    shell:
        "/project/pacbiosv/bin/GraphAligner/bin/GraphAligner \
         --seeds-minimizer-count 20 -g {input.gfa} -f {input.fasta} -a {output.json}"


rule format_json:
    input:
        json = "regions/json_{technology}/r{region}.json"
    output:
        json = "regions/json_{technology}/r{region}.parsed.json"
    shell:
        "cat {input.json} | jq -s '[ [ [ .[] | {{ read: .name, identity: .identity, length: .sequence|length, nodes: [.path.mapping[].position.name]}} | \
         select(.nodes | length > 1)] | group_by(.read) | .[] | max_by(.length * .identity)] | .[] | select(.identity >= 0.8) ]' > {output.json}"

rule prune_with_nanopore:
    input:
        gfa = "regions/gfas/r{region}.reducted.gfa",
        json = "regions/json_nanopore/r{region}.parsed.json"
    output:
        gfa = "regions/gfas_pruned_nanopore/r{region}.min{min_coverage}.gfa"
    shell:
        "python3 ../scripts/use_read_alignments.py {input.gfa} {input.json} remove --min_cov {wildcards.min_coverage} > {output.gfa}"

rule prune_with_nanopore_and_pacbio:
    input:
        gfa = "regions/gfas/r{region}.reducted.gfa",
        nanopore = "regions/json_nanopore/r{region}.parsed.json",
        pacbio = "regions/json_pacbio/r{region}.parsed.json"
    output:
        gfa = "regions/gfas_pruned_nanopore_pacbio/r{region}.min{min_coverage}.gfa"
    shell:
        "python3 ../scripts/use_read_alignments.py {input.gfa} {input.nanopore} remove --additional_json {input.pacbio} --min_cov {wildcards.min_coverage} > {output.gfa}"


rule prune_small_bubbles:
    input:
        gfa = "regions/gfas_pruned_nanopore/r{region}.min{min_coverage}.gfa"
    output:
        gfa = "regions/gfas_pruned2/r{region}.min{min_coverage}.max{max_size}.gfa"
    shell:
        "python3 ../scripts/find_ultrabubbles.py {input.gfa} remove --max_size {wildcards.max_size} > {output.gfa}"

rule plot_bandage_raw:
    input:
        graph = "regions/gfas/r{region}.reducted.gfa"
    output:
        png = "regions/pngs/r{region}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then /project/pacbiosv/bin/bandage/Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"

rule plot_bandage_pruned:
    input:
        graph = "regions/gfas_{dir_suffix}/r{region}.{file_suffix}.gfa"
    output:
        png = "regions/pngs/r{region}_{dir_suffix,pruned_nanopore|pruned_nanopore_pacbio}.{file_suffix}.png"
    shell:
        "LINES=$(wc -l < {input.graph}) ; if [ $LINES -gt 0 ]; then /project/pacbiosv/bin/bandage/Bandage image {input.graph} {output.png} --height 4000; else echo '' > {output.png}; fi"


rule extract_contigs:
    input:
        graph = "regions/gfas/r{region}.reducted.gfa"
    output:
        "regions/contigs/r{region}.fa"
    shell:
        "python3 ../scripts/GFAextract_sequence.py {input.graph} > {output}"

rule pool_contigs:
    input:
        expand("regions/contigs/r{i}.fa", i=range(1, 505))
    output:
        "regions/contigs/pooled.fa"
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
