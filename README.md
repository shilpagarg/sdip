# SDip: A novel graph-based approach to haplotype-aware assembly based structural variant calling in segmental duplications 

Segmental duplications are important for understanding human diseases and evolution. The challenge to distinguish allelic and duplication sequences has hindered their phased assembly as well as characterization of structural variant calls. Here we have developed a novel graph-based approach that leverages single nucleotide differences in overlapping reads to distinguish allelic and duplication sequences information from long read accurate PacBio HiFi sequencing. These differences enable to generate allelic and duplication-specific overlaps in the graph to spell out phased assembly used for structural variant calling. We have applied our method to three public genomes: CHM13, NA12878 and HG002. Our method resolved 86% of collapsed duplicated regions with contig N50 up to 79 kb and 22 structural variant phased calls, outperforming state-of-the-part SDA method in terms of all metrics. Furthermore, we demonstrate the importance of phased assemblies and variant calls to the biologically-relevant duplicated genes such as SMN1, HYDIN, NPY4R and FAM72A. Thus, our phased assemblies and accurate variant calling enable the study of the impact of duplications on human evolution and diseases.


Installation using conda.
```
conda env create -n sdip -f config.yaml
```

Run whole analysis.
```
snakemake -j 30
```

Generate HiFi base graph. 
```
usage: haplotype.py [-h] [-p PAF] [-o GFA] [-t INT] [--het FLOAT]
                          [-c INT]
                          FASTA
positional arguments:
  FASTA                 Input all read records in one FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  -p PAF, --paf PAF     PAF file, should be consistent with the input FASTA
                        file; if not given, we will run minimap2 on the input
                        fasta and generate one.
  -o GFA, --output GFA  Output to FILE, default to stdout.
  -t INT, --threads INT
                        Maximum number of threads to use. [4]
  --het FLOAT           Heterozygosity rate, set as a filter for abnormal
                        alignment. [0.1]
  -c INT, --avg-coverage INT
                        average coverage, set as a filter for abnormal
                        alignment. [40]
```

An example run on our test data:
```
python haplotype.py -t 40 test.fasta > test.gfa
```
