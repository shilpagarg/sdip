import networkx as nx
import pandas as pd
import numpy as np
from collections import defaultdict

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
        regions = config["regions"]["slop15k"],
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
        expand("regions/contigs/r{i}.t{{tip_max_size}}.b{{bubble_max_size}}.d{{degree_max_size}}.polished.fa", i=good_regions)
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa"
    shell:
        "cat {input} > {output}"

rule map_contigs:
    input:
        fa = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa",
        ref = config["genome"]
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bam"
    threads: 10
    shell:
        "minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
        {input.ref} {input.fa} | samtools view -F 260 -u - | samtools sort - > {output}"

rule index_bam:
    input:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bam"
    output:
        "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bam.bai"
    shell:
        "samtools index {input}"

##########
#Evaluate#
##########

rule bam_to_bed:
    input:
        bam = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bam"
    output:
        bed = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bed"
    shell:
        "bedtools bamtobed -i {input.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}"

rule sd_align_to_ref:
    input:
        bed = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.sorted.bed",
        regions = config["regions"]["original"],
        index = config["genome"] + ".fai"
    output:
        inter = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/inter.bed",
        rbed = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/SD.resolved.bed",
        ubed = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/SD.unresolved.bed",
        sd_bed = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/SD.status.tbl"
    shell:
        "bedtools intersect -a {input.bed} -b {input.regions} -wa -wb > {output.inter} && %s/PrintResolvedSegdups.py {output.inter} {input.regions} {input.index}\
        --extra 10000 --resolved {output.rbed} --unresolved {output.ubed} --allseg {output.sd_bed}" % (config["tools"]["paftest"])

rule align_bac_to_asm:
    input:
        asm = "regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.polished.fa",
        bacs = config["bacs"]["fasta"]
    output:
        bam = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/bacs/bacs.to.asm.bam"
    threads: 10
    shell:
        "minimap2 -I 8G -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
        {input.asm} {input.bacs} | samtools view -F 2308 -u - | samtools sort - > {output.bam}"

rule bac_to_asm_tbl:
    input:
        bam = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/bacs/bacs.to.asm.bam",
        names = config["bacs"]["names"]
    output:
        tbl = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/bacs/bacs.to.asm.tbl"
    shell:
        "python3 %s/samIdentity.py --header --mask {input.names} {input.bam} > {output.tbl}" % (config["tools"]["paftest"])

rule make_qv_sum:
    input:
        bac_tbl = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/bacs/bacs.to.asm.tbl"
    output:
        qv_sum = "regions/eval/t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}/bacs/qv_sum.txt"
    run:
        pd.options.mode.use_inf_as_na = True
        out = ""
        sys.stderr.write(input["bac_tbl"] + "\n")
        df = pd.read_csv(input["bac_tbl"], sep="\t")
        val = 1 - df["perID_by_all"]/100
        df["qv"] = -10 * np.log10( val )
        for mask in [True, False]:            
            tmp = df[df["mask"] == mask]
            perfect = tmp["qv"].isna()
            out += "{}\nPerfect\t{}\n{}\n\n".format(input["bac_tbl"], sum(perfect), tmp.qv.describe()   )
        open(output["qv_sum"], "w+").write(out)
