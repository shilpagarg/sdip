# The main entry point of this workflow.
# After configuring the workflow in config/config.yaml, running snakemake --use-conda should successfully execute the workflow.

#include configuration file
configfile: "config/config.yaml"

#include collection of rules
include: "workflow/rules/rules.smk"

SAMPLES = list(config["reads"].keys())

rule all:
    input:
        expand("{sample}/{sample}.polished.diploid.fa", sample=SAMPLES),
        expand("{sample}/sv_calls_diploid/variants.vcf", sample=SAMPLES)
