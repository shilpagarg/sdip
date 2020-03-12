# SDip: A novel graph-based approach to haplotype-aware structural variant calling in segmental duplication regions

Segmental duplications play an important role in human evolution and disease. The high similarity between allelic and paralogous sequences, however, has hindered their diploid assembly as well as the characterization of structural variants. We have developed a novel graph-based approach, SDip, that leverages single nucleotide differences in overlapping reads to distinguish allelic and duplicated sequences from accurate long-read PacBio HiFi sequencing. These differences enable the generation of an allele and paralogue-specific overlap graph to spell out a diploid assembly used for structural variant calling.

We have applied our method to three public human genome assemblies for CHM13, NA12878 and HG002. Our method fully resolved more than 80% of duplicated regions with a contig N50 of up to 82 kbps and produced >700 phased structural variant calls, outperforming state-of-the-art method SDA across all metrics. Furthermore, we demonstrate the importance of phased assemblies and variant calls for the characterization of biologically relevant duplicated genes such as SRGAP2C, SMN2, NPY4R and FAM72A. Our diploid assemblies and accurate variant calling specifically in duplicated regions will enable the study of the evolution and adaptation of various species.

See our preprint here: [https://doi.org/10.1101/2020.02.25.964445](https://doi.org/10.1101/2020.02.25.964445).

This repository contains the complete SDip pipeline and can be used to reproduce the results from our manuscript. 

## Installation

SDip has been implemented as a [snakemake](https://snakemake.readthedocs.io) pipeline. All required dependencies are automatically deployed by snakemake in isolated [conda](https://docs.conda.io/en/latest/) environments.

1. To run SDip on a given base assembly you first need to install conda following the instructions [here](https://conda.io/en/latest/miniconda.html). Make sure to install the **Python 3** version of miniconda.
2. Then, you can [install](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) snakemake using conda: `conda install -c bioconda -c conda-forge snakemake`
3. Finally, you need to download the SDip pipeline from this github repository:
```
# clone workflow into working directory
git clone https://github.com/shilpagarg/sdip.git path/to/workdir
cd path/to/workdir
```

## Configuration

Now that SDip is installed, you need to customize the configuration script under `config/config.yaml` for your base assembly. In particular, the following input files need to be specified:

* **base**: the base assembly whose collapsed regions shall be assembled 
* **hg38**: the reference genome used for evaluation (e.g. from [here](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz))
* **regions**: bed file of collapsed regions in base assembly (e.g. from SDA)
* **pacbio**: FASTA file of PacBio HiFi reads
* **nanopore**: FASTA file of Ultralong Nanopore reads
* **bacs**: FASTA file of BAC sequences for evaluation

You also need to determine the number of regions in the bed file of collapsed regions and paste it into the `Snakefile` as is indicated there (step 2).

## Execution

Now, the first stage of the workflow can be run. It aligns the collapsed regions to the base assembly to find similar regions:

```
# run pipeline until the file segdups.similar.tsv is generated (-j specifies the number of threads)
snakemake --use-conda segdups.similar.tsv -j 30
```

When this stage of the pipeline finished successfully, you need to determine the number of region clusters found. It is the number in the first column of the last line in the `segdups.similar.tsv` file. Paste the numnber into the `Snakefile` as is indicated there (step 3). Then execute the rest of the pipeline:

```
# run the remainder of the pipeline (-j specifies the number of threads)
snakemake --use-conda -j 30
```

While the pipeline is running, status message will be printed on stdout and saved to log files in `.snakemake/log/`. You can also keep an eye on the currently running jobs using `top` or `htop`.

The processing of some regions might fail e.g. because of cycles in the graph. Other regions might take too long to process. You can add the indices of these regions in the `Snakefile` (step 4) to skip them.

Also have a look into the [documentation](https://snakemake.readthedocs.io/en/stable/) of snakemake to learn more about its different modes and options.

## Results

After the workflow completes, you will find the final diploid contigs in `regions/contigs`:

- All contigs: `regions/contigs/pooled.t5.b5000.d2.polished.grouped.fa`
- Contigs from first chromosomal haplotype: `regions/contigs/pooled.t5.b5000.d2.polished.haplotype1.fa`
- Contigs from second chromosomal haplotype: `regions/contigs/pooled.t5.b5000.d2.polished.haplotype2.fa`
- Contigs from both chromosomal haplotypes: `regions/contigs/pooled.t5.b5000.d2.polished.diploid.fa`
- Contigs that could not be matched into haplotype pairs: `regions/contigs/pooled.t5.b5000.d2.polished.haplotype0.fa`

The diploid SV calls detected in alignments of both chromosomal haplotypes to the reference genome can be found here: `regions/svim/t5.b5000.d2/contigs.diploid/variants.norm.vcf`

Evaluation results for the contigs can be found in `regions/eval/t5.b5000.d2`:

- Contig statistics using quast: `regions/eval/t5.b5000.d2/quast_to_bacs/grouped/report.html`
- QV of contigs versus BACs: `regions/eval/t5.b5000.d2/tables/qv_sum.grouped.txt`
- Statistics on fully assembled regions (number of covering contigs: count): `regions/eval/t5.b5000.d2/resolved.grouped/stats.txt`
- List of misassembled contigs: `regions/eval/t5.b5000.d2/misassemblies.grouped/misassembled.txt`


## Contact

If you experience problems or have suggestions please create an issue or a pull request or contact heller_d@molgen.mpg.de.
