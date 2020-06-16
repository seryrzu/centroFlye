# centroFlye

**Disclaimer:** The following description is aimed for reconstruction of assemblies presented in the Nature Biotechnology paper (see `Publications`).
The latest stable assemblies are produced by a more recent version of centroFlye that can be found in `master` branch.
[This](https://github.com/seryrzu/centroFlye_paper_scripts) github repository provides all supporting jupyter notebooks that are neccessary to replicate the results in the paper.
For convenience that repository includes current branch of the main centroFlye repository as a submodule.


### Version 0.8.3 cenX + 0.1.3 cen6
## Overview
centroFlye (Bzikadze et al., 2020) is an algorithm for centromere assembly using long error-prone reads.
Currently it supports assembly of a human centromere 6 (referred to as cen6) and X (referred to as cenX).
Here we show how to apply it for the cen6 and cenX of the CHM13hTERT human cell line.
The comparison of various cenX assemblies can be found in Bzikadze et al., 2020 and Mikheenko et al., 2020.

# Cloning
Please, clone the repository with submodules:
```
git clone --recurse-submodules --single-branch --branch cF_NatBiotech_paper_Xv0.8.3-6v0.1.3 git@github.com:seryrzu/centroFlye.git
```

## Dependencies

+ C++ 14
+ Python 3.6

### Python packages:
+ [Biopython](https://pypi.org/project/biopython/) (tested on version 1.70)
+ [Edlib](https://pypi.org/project/edlib/) (tested on version 1.2.3)
+ [Networkx](https://pypi.org/project/networkx/) (tested on version 2.2)
+ [Numpy](https://pypi.org/project/numpy/) (tested on version 1.16.1)
+ [Regex](https://pypi.org/project/regex/) (tested on version 2.4.104)

Required python packages can be installed through Conda with `conda install --file requirements.txt`.

### External software
+ [Flye](https://github.com/fenderglass/Flye) (*v2.5*, tested on commit: `315122d2ff58025aa4c38227239f431490b557ac`)
+ [Noise Cancelling Repeat Finder (NCRF)](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder) (Tested on commit : `758206f1689ad1338cf7a841482dbf12548c337a`)

Please note that all external software by default has to be in your `PATH`.

### Data
+ `<path_to_CHM13>` — path where the T2T ONT reads are located (rel2, Guppy flip-flop 2.3.1, used for cenX assembly and can be downloaded from [here](https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel2/rel2.fastq.gz); rel3, Guppy flip-flop 3.1.5, is used for cen6 assembly and can be downloaded from [here](https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/rel3.fastq.gz) ; also see [github](https://github.com/nanopore-wgs-consortium/CHM13)). The data is described in Miga, Koren et al., 2020.

## Availability
Final assembly and all intermediate results of the pipeline described below are published at [ZENODO](https://doi.org/10.5281/zenodo.3369553)

## Quick start guide for cenX (centroFlye)

If you want to run the whole centroFlye pipeline to generate cenX assembly with one command, please run
```
./run_all_cenX.sh <path_to_CHM13> results_cenX 50
```
You can customize the output directory and the number of threads:
```
./run_all_cenX.sh <path_to_CHM13> <output directory> <number of threads>
```
All intermediate and final results will then be placed in `<output directory>`.
If you want to start from scratch you can simply remove this directory.

**Required resources**:
+ Storage space: ~150GB (mostly from the first step "Recruitment of centromeric reads", see below)
+ Clock time: ~9 hours (mostly recruitment of unique k-mers)
+ RAM: peak usage up to 800GB


## Pipeline for cenX (centroFlye)
In this manual we go step-by-step demonstrating centroFlye algorithm.
The detailed information about the algorithm can be found in the paper.

Please, run all commands from the root of the repository.
Results of all steps will be stored at `results_cenX` directory.

### 1. Recruitment of cenX reads

We use 50x ultra-long Oxford Nanopore dataset generated by [Telomere2Telomere Consorsium](https://github.com/nanopore-wgs-consortium/CHM13), rel2.
This step is run directly on the reads at the [link](https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel2/rel2.fastq.gz).
The following bash script splits the input file in 50 files and runs DXZ1-based recruitment in 50 threads.
DXZ1 is supplied in the current repo at ``supplementary_data/DXZ1_rc.fasta``.
The result of this step is a fasta file with centromeric reads that is stored at `results_cenX/centromeric_reads/centromeric_reads.fasta`

From the root of the project run 
```
make -C scripts/read_recruitment
```
and start recruitment (`<path_to_CHM13>` is where the ONT reads are located, see section `Dependencies/Data`):
```
bash scripts/read_recruitment/run_read_recruitment.sh \
       <path_to_CHM13>/rel2.fastq.gz \
       results_cenX/centromeric_reads 50 11100000
```
**Required resources**:
+ Storage space: 150GB
+ Clock time: 1 hour
+ RAM: < 50MB

### 2. Applying centroFlye to cenX reads
At this step we are utilizing centromeric reads from step 1 and run centroFlye on them.
The result of this step is the final cenX assembly that is stored at `results_cenX/final_assembly.fasta`.
The following command uses 50 threads.

```
python centroFlye.py \
            --reads results_cenX/centromeric_reads/centromeric_reads.fasta \
            -t 50 \
            --outdir results_cenX \
            --unit supplementary_data/DXZ1_rc.fasta
```
**Required resources**:
+ Storage space: < 3GB
+ Clock time: ~9 hours (mostly recruitment of unique k-mers)
+ RAM: peak usage up to 800GB


## Quick start guide for cen6 (centroFlyeMono)

centroFlye has a special module centroFlyeMono that uses only alpha satellite structure to assemble centromeres.
The steps below describe how to apply it to cen6 assembly.

If you want to run the whole centroFlyeMono pipeline to generate cen6 assembly with one command, please run
```
./run_all_cen6.sh <path_to_CHM13> results_cenX 50
```
You can customize the output directory and the number of threads:
```
./run_all_cen6.sh <path_to_CHM13> <output directory> <number of threads>
```
All intermediate and final results will then be placed in `<output directory>`.
If you want to start from scratch you can simply remove this directory.


## Pipeline for cen6 (centroFlyeMono)
In this manual we go step-by-step demonstrating centroFlyeMono algorithm.
The detailed information about the algorithm can be found in the paper.

Please, run all commands from the root of the repository.
Results of all steps will be stored at `results_cen6` directory.


### 1. Recruitment of cen6 reads

We use 120x ultra-long Oxford Nanopore dataset generated by [Telomere2Telomere Consorsium](https://github.com/nanopore-wgs-consortium/CHM13), rel3.
This step is run directly on the reads at the [link](https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/rel3.fastq.gz).
The following bash script splits the input file in 50 files and runs D6Z1-based recruitment in 50 threads.
D6Z1 is supplied in the current repo at ``supplementary_data/D6Z1.fasta``.
The result of this step is a fasta file with centromeric reads that is stored at `results_cenX/centromeric_reads/centromeric_reads.fasta`

From the root of the project run 
```
make -C scripts/read_recruitment
```
and start recruitment (`<path_to_CHM13>` is where the ONT reads are located, see section `Dependencies/Data`):
```
bash scripts/read_recruitment/run_read_recruitment.sh \
       <path_to_CHM13>/rel3.fastq.gz \
       results_cen6/centromeric_reads 50 29000000 \
       supplementary_data/D6Z1.fasta
```
**Required resources**:
+ Storage space: 150GB
+ Clock time: 1 hour
+ RAM: < 50MB


### 2. Running String Decomposer on cen6 reads

We use String Decomposer (SD; Dvorkina et al., 2020) to partition cen6 reads into distinct monomers of D6Z1.
These monomers are supplied in the current repo at ``supplementary_data/D6Z1_monomers.fasta``.
The result of this step is the report of SD is stored at `results_cen6/string_decomposer_report`

The following commands take 50 threads
```
python scripts/ext/stringdecomposer/longreads_decomposer.py \
            -s results_cen6/centromeric_reads/centromeric_reads.fasta \
            -m supplementary_data/D6Z1_monomers.fasta \
            -t 50
mkdir results_cen6/string_decomposer_report
mv decomposition{.tsv,_alt.tsv} results_cen6/string_decomposer_report
```
**Required resources**:
+ Storage space: 6GB
+ Clock time: 9h
+ RAM: < 10GB


### 3. Running centroFlyeMono on cen6 reads using SD output

At this step we are utilizing string decomposer output on cen6 reads and run centromFlyeMono on them.
The final cen6 assembly is located in `results_cen6/centroFlyeMono_cen6/polishing/scaffold_0/scaffold_0.fasta`.
```
python scripts/centroFlyeMono.py \
       --sd-report results_cen6/string_decomposer_report/decomposition.tsv \
       --monomers supplementary_data/D6Z1_monomers.fasta \
       --centromeric-reads results_cen6/centromeric_reads/centromeric_reads.fasta \
       --outdir results_cen6/centroFlyeMono_cen6
```
**Required resources**:
+ Storage space: 120MB
+ Clock time: 30 mins
+ RAM: < 1GB


## Publications
- Bzikadze A.V., Pevzner P.A. centroFlye: Assembling Centromeres with Long Error-Prone Reads, *Nature Biotechnology, in press*, 2020
- Dvorkina T., Bzikadze A.V., Pevzner P.A. The String Decomposition Problem and its Applications to Centromere Assembly, *Bioinformatics, in press*, 2020
- Miga, K.H., Koren, S., Rhie, A., Vollger, M.R., Gershman, A., Bzikadze, A., Brooks, S., Howe, E., Porubsky, D., Logsdon, G.A., et al. Telomere-to-telomere assembly of a complete human X chromosome, *Nature, in press*, 2020
- Mikheenko, A., Bzikadze, A.V., Gurevich, A., Miga, K.H., and Pevzner, P.A. TandemMapper and TandemQUAST: mapping long reads and assessing/improving assembly quality in extra-long tandem repeats, *Bioinformatics, in press*, 2020


## Contacts
Please report any problems to the [issue tracker](https://github.com/seryrzu/centroFlye/issues).
Alternatively, you can write directly to [abzikadze@ucsd.edu](mailto:abzikadze@ucsd.edu).
