# SDip: A novel graph-based approach to haplotype-aware assembly based structural variant calling in segmental duplications 

Segmental duplications play an important role in human evolution and disease. The high similarity between allelic and paralogous sequences, however, has hindered their diploid assembly as well as the characterization of structural variants. Here we have developed a novel graph-based approach that leverages single nucleotide differences in overlapping reads to distinguish allelic and duplicated sequences from accurate long-read PacBio HiFi sequencing. These differences enable the generation of an allele and paralogue-specific overlap graph to spell out a diploid assembly used for structural variant calling. We have applied our method to three public human genome assemblies for CHM13, NA12878 and HG002. Our method fully resolved more than 80% of duplicated regions with a contig N50 of up to 82 kbps and produced >700 phased structural variant calls, outperforming state-of-the-art method SDA across all metrics. Furthermore, we demonstrate the importance of phased assemblies and variant calls for the characterization of biologically relevant duplicated genes such as SRGAP2C, SMN2, NPY4R and FAM72A. Our diploid assemblies and accurate variant calling specifically in duplicated regions will enable the study of the evolution and adaptation of various species.


## Installation

SDip has been implemented as a [snakemake](https://snakemake.readthedocs.io) pipeline.  To run it, you first need to clone this github repository:

```
# clone workflow into working directory
git clone https://github.com/shilpagarg/sdip.git path/to/workdir
# switch to branch "tool"
git checkout tool
cd path/to/workdir
```

Most dependencies of SDip are automatically deployed using conda. Two dependencies, however, need to be installed manually. To install svim-asm:

```
# clone repository for svim-asm
git clone https://github.com/eldariont/svim-asm.git path/to/install-dir
cd path/to/install-dir
# install using pip, make sure that svim-asm ends up to be in the PATH
pip install .
```

The second dependency to install is mc-mpc:

```
# clone repository for mc-mpc
git clone https://github.com/tobtobtob/MC-MPC.git path/to/install-dir2
cd path/to/install-dir2/mc-mpc-solver
# compile mc-mpc
make
# add to PATH
```

Now that SDip is installed, you need to customize the configuration script under `config/config.yaml`. In particular, the input files need to specified.

```
# edit config file: add paths to pacbio reads, nanopore reads and reference genome 
vim config/config.yaml

# execute workflow, deploy software dependencies via conda
snakemake --use-conda
```

Then, the workflow can be run:

```
# first, a dry run
snakemake -n --use-conda

# run for real
snakemake --use-conda
```

After the workflow completes, you will find the following output files in your working directory:

Final polished contigs: `"<sample_name>/<sample_name>.polished.grouped.fa"`
SV calls including genotypes:
`"<sample_name>/sv_calls_diploid/variants.vcf"`
