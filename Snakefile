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
big_regions = []
#Range of regions with cycles: should be initially empty, later insert regions that have an empty .lemon graph in regions/gfas/pruned
cycle_regions = []
#Range of duplicate regions: should be initially empty, later insert regions that are contained by other regions (use find_redundant_segdup_regions.sh)
duplicate_regions = []
# Range of regions that produce empty contigs
empty_contigs = []
good_regions = [r for r in filtered_segdup_regions if (not r in big_regions) and (not r in cycle_regions) and (not r in duplicate_regions) and (not r in empty_contigs)]

#include rules.smk file from local clone of sdip repository
include: "/path/to/sdip/rules.smk"

rule all:
    input:
        expand("regions/pngs/r{i}.png", i=good_regions),
        expand("regions/pngs/r{i}.pruned.t5.b5000.d2.png", i=good_regions),
        "regions/stats/stats.t5.b5000.d2.txt",
        expand("regions/eval/t5.b5000.d2/quast_to_bacs/{variant}/report.html", variant=["grouped", "haplotype1", "haplotype2"]),
        expand("regions/eval/t5.b5000.d2/quast_to_hg38/{variant}/report.html", variant=["grouped", "haplotype1", "haplotype2"]),
        "regions/eval/t5.b5000.d2/resolved.grouped/inter.bed",
        "regions/eval/t5.b5000.d2/misassemblies.grouped/confirmed.txt",
        "regions/eval/t5.b5000.d2/tables/qv_sum.grouped.txt",
        "regions/eval/t5.b5000.d2/tables/polished.grouped.to.hg38.tbl",
        "regions/eval/t5.b5000.d2/tables/polished.grouped.to.bacs.tbl",
        "regions/svim/t5.b5000.d2/contigs.diploid/variants.norm.vcf.gz",
        "regions/svim/t5.b5000.d2/bacs.haploid/variants.norm.vcf.gz",

        "regions/eval/sda/quast_to_bacs/haploid/report.html",
        "regions/eval/sda/quast_to_hg38/haploid/report.html",
        "regions/eval/sda/resolved.haploid/inter.bed",
        "regions/eval/sda/misassemblies.haploid/confirmed.txt",
        "regions/eval/sda/tables/qv_sum.haploid.txt",
        "regions/eval/sda/tables/polished.haploid.to.hg38.tbl",
        "regions/eval/sda/tables/polished.haploid.to.bacs.tbl",
        "regions/svim/sda/contigs.haploid/variants.norm.vcf.gz",
        "regions/svim/sda/bacs.haploid/variants.norm.vcf.gz"
