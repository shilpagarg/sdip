"""
This is a template Snakefile that should be copied to the respective working directory together with the config.yaml.
Please customize this Snakefile (particularly the regions and the target files in "rule all").
Please also customize the paths in the config.yaml
"""

configfile: "config.yaml"

big_regions = [17, 32, 96, 101, 102, 108, 162, 166, 216, 217, 218, 267, 268, 275, 278, 285, 296, 306, 307, 433, 467, 470, 471, 503, 504]
cycle_regions = [136, 187, 188, 189, 21, 227, 266, 269, 276, 333, 334, 34, 347, 361, 388, 423, 435, 440, 476, 488, 491, 500]
regions = [elem for elem in list(range(1, 505)) if (not elem in big_regions) and (not elem in cycle_regions)]

rule all:
    input:
        expand("regions/pngs/r{i}.png", i=regions),
        expand("regions/pngs/r{i}.pruned.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.png", i=regions,
                                                                                                     tip_max_size=[5],
                                                                                                     bubble_max_size=[5],
                                                                                                     degree_max_size=[2]),
        expand("regions/contigs/pooled.t{tip_max_size}.b{bubble_max_size}.d{degree_max_size}.fa", tip_max_size=[5],
                                                                                                  bubble_max_size=[5],
                                                                                                  degree_max_size=[2]),
        "regions/contigs/pooled.sorted.bam.bai",
        "regions/stats/cover.t5.b5.d2.txt"
