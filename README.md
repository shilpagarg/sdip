# SDip: A novel graph-based approach to haplotype-aware structural variant calling in segmental duplication regions

Segmental duplications play an important role in human evolution and disease. The high similarity between allelic and paralogous sequences, however, has hindered their diploid assembly as well as the characterization of structural variants. We have developed a novel graph-based approach, SDip, that leverages single nucleotide differences in overlapping reads to distinguish allelic and duplicated sequences from accurate long-read PacBio HiFi sequencing. These differences enable the generation of an allele and paralogue-specific overlap graph to spell out a diploid assembly used for structural variant calling.

We have applied our method to three public human genome assemblies for CHM13, NA12878 and HG002. Our method fully resolved more than 80% of duplicated regions with a contig N50 of up to 82 kbps and produced >700 phased structural variant calls, outperforming state-of-the-art method SDA across all metrics. Furthermore, we demonstrate the importance of phased assemblies and variant calls for the characterization of biologically relevant duplicated genes such as SRGAP2C, SMN2, NPY4R and FAM72A. Our diploid assemblies and accurate variant calling specifically in duplicated regions will enable the study of the evolution and adaptation of various species.

See our preprint here: [https://doi.org/10.1101/2020.02.25.964445](https://doi.org/10.1101/2020.02.25.964445).

This branch (tool) of the repository contains only the SDip tool and can be used to generate diploid contigs and SV calls from a set of PacBio HiFi and Ultralong Nanopore reads. 

## Installation

SDip has been implemented as a [snakemake](https://snakemake.readthedocs.io) pipeline. All required dependencies are automatically deployed by snakemake in isolated [conda](https://docs.conda.io/en/latest/) environments.

1. To run SDip on a given base assembly you first need to install conda following the instructions [here](https://conda.io/en/latest/miniconda.html). Make sure to install the **Python 3** version of miniconda.
2. Then, you can [install](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) snakemake using conda: `conda install -c bioconda -c conda-forge snakemake`
3. Finally, you need to download the SDip pipeline from this github repository and switch to the `tool` branch:
```
# clone workflow into working directory
git clone https://github.com/shilpagarg/sdip.git path/to/workdir
cd path/to/workdir
# switch to branch "tool"
git checkout tool
```

## Configuration

Now that SDip is installed, you need to customize the configuration script under `config/config.yaml` for your base assembly. In particular, the following input files need to be specified:

* **reads**: FASTA files of PacBio HiFi and Ultralong Nanopore reads for each sample
* **reference**: the reference genome used for evaluation (e.g. from [here](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz))

## Execution

Now, the workflow can be executed:

```
# run the pipeline (-j specifies the number of threads)
snakemake --use-conda -j 30
```

While the pipeline is running, status message will be printed on stdout and saved to log files in `.snakemake/log/`. You can also keep an eye on the currently running jobs using `top` or `htop`.

Also have a look into the [documentation](https://snakemake.readthedocs.io/en/stable/) of snakemake to learn more about its different modes and options.

## Results

After the workflow completes, you will find the following output files in your working directory:

- Final polished contigs: `"<sample_name>/<sample_name>.polished.grouped.fa"`
- SV calls including genotypes:
`"<sample_name>/sv_calls_diploid/variants.vcf"`
