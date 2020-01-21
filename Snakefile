"""
This is a template Snakefile that should be copied to the respective working directory together with the config.yaml.
Please customize this Snakefile (particularly the regions and the target files in "rule all").
Please also customize the paths in the config.yaml
"""

#include configuration file from same directory
configfile: "config.yaml"

#Range of regions in regions file that is defined in config file
all_segdup_regions = range(1, 605)

#Range of regions in segdups.similar.tsv
filtered_segdup_regions = range(1, 505)
#Range of big regions: should be initially empty, later insert regions that take too long to process
big_regions = [17, 32, 96, 101, 102, 108, 162, 166, 216, 217, 218, 267, 268, 275, 278, 285, 296, 306, 307, 433, 467, 470, 471, 503, 504]
#Range of regions with cycles: should be initially empty, later insert regions that have an empty .lemon graph in regions/gfas/pruned
cycle_regions = [136, 187, 188, 189, 21, 227, 266, 269, 276, 333, 334, 34, 347, 361, 388, 423, 435, 440, 476, 488, 491, 500]
#Range of duplicate regions: should be initially empty, later insert regions that are contained by other regions (use find_redundant_segdup_regions.sh)
duplicate_regions = []
good_regions = [r for r in filtered_segdup_regions if (not r in big_regions) and (not r in cycle_regions) and (not r in duplicate_regions)]

#include rules file from local clone of WHdenovo
include: "/lustre/scratch115/realdata/mdt3/projects/graphs/yeast/erik/final_exp/segdup/david/WHdenovo/paftest/rules.smk"

rule all:
    input:
        "segdups.similar.tsv",
        expand("regions/pngs/r{i}.png", i=good_regions),
        expand("regions/pngs/r{i}.pruned.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.png", i=good_regions,
                                                                                                     tip_max_size=[5],
                                                                                                     bubble_max_size=[5],
                                                                                                     degree_max_size=[2]),
        "regions/stats/stats.t5.b5.d2.txt",
        "regions/eval/t5.b5.d2/bams/polished.to.bacs.bam",
        "regions/eval/t5.b5.d2/quast_to_assembly/report.html",
        "regions/eval/t5.b5.d2/quast_to_bacs/report.html",
        "regions/eval/t5.b5.d2/quast_to_hg38/report.html",
        "regions/eval/t5.b5.d2/resolved/inter.bed",
        "regions/eval/t5.b5.d2/bacs/qv_sum.txt",
        "regions/svim/t5.b5.d2/contigs.haploid/variants.vcf",
        "regions/svim/t5.b5.d2/bacs.haploid/variants.vcf",

        "regions/eval/sda/quast_to_assembly/report.html",
        "regions/eval/sda/quast_to_bacs/report.html",
        "regions/eval/sda/quast_to_hg38/report.html",
        "regions/eval/sda/resolved/inter.bed",
        "regions/eval/sda/bacs/qv_sum.txt",
        "regions/svim/sda/contigs.haploid/variants.vcf",
        "regions/svim/sda/bacs.haploid/variants.vcf"
