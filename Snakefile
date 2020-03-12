# The main entry point of this workflow.

################################################################
# 1. Before execution, first configure the config/config.yaml. #
################################################################
configfile: "config/config.yaml"

##########################################################################
# 2. Determine the number C of collapsed regions in your input bed file. #
#    Change the range below to (1, C+1)                                  #
##########################################################################

#Range of regions in regions file that is defined in config file
all_segdup_regions = range(1, 487)

#################################################
# 3. Then run first stage of the pipeline.      #
#    Determine the number S of region clusters. #
#    Change the range below to (1, S+1)         #
#################################################

#Range of regions in segdups.similar.tsv
filtered_segdup_regions = range(1, 467)

##########################################################################################
# 4. Finally execute the remainder of the pipeline.                                      #
#    Some regions might fail because of cycles in the graph or take too long to process. #
#    Add their numbers below to skip them.                                               #
##########################################################################################

#Range of big regions: should be initially empty, later insert regions that take too long to process
big_regions = []
#Range of regions with cycles: should be initially empty, later insert regions that have an empty .lemon graph in regions/gfas/pruned
cycle_regions = []
#Range of duplicate regions: should be initially empty, later insert regions that are contained by other regions (use find_redundant_segdup_regions.sh)
duplicate_regions = []
# Range of regions that produce empty contigs
empty_contigs = []
good_regions = [r for r in filtered_segdup_regions if (not r in big_regions) and (not r in cycle_regions) and (not r in duplicate_regions) and (not r in empty_contigs)]

#include collection of rules
include: "workflow/rules/rules.smk"

rule all:
    input:
        expand("regions/pngs/r{i}.png", i=good_regions),
        expand("regions/pngs/r{i}.pruned.t5.b5000.d2.png", i=good_regions),
        "regions/stats/stats.t5.b5000.d2.txt",
        expand("regions/eval/t5.b5000.d2/quast_to_bacs/{variant}/report.html", variant=["grouped", "haplotype1", "haplotype2"]),
        expand("regions/eval/t5.b5000.d2/quast_to_hg38/{variant}/report.html", variant=["grouped", "haplotype1", "haplotype2"]),
        "regions/eval/t5.b5000.d2/resolved.grouped/resolved.txt",
        "regions/eval/t5.b5000.d2/misassemblies.grouped/confirmed.txt",
        "regions/eval/t5.b5000.d2/tables/qv_sum.grouped.txt",
        "regions/eval/t5.b5000.d2/tables/polished.grouped.to.hg38.tbl",
        "regions/eval/t5.b5000.d2/tables/polished.grouped.to.bacs.tbl",
        "regions/svim/t5.b5000.d2/contigs.diploid/variants.norm.vcf.gz",
        "regions/svim/t5.b5000.d2/bacs.haploid/variants.norm.vcf.gz"
