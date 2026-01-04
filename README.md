## Overview
The core idea of RAGA is to generate alternative long reads by leveraging the target assembly and PacBio HiFi sequencing reads, with the assistance of a reference genome. Users can follow the guidance provided in the "Usage" section to achieve this goal through our analysis pipeline.

![workflow](https://github.com/wzxie/RAGA/blob/main/workflow.jpg)

## Dependencies
1. minimap2 (2.22-r1101)
2. racon (v1.4.20)
3. ragtag.py (v2.1.0)
4. nucmer (4.0.0rc1)
5. awk (GNU Awk 4.0.2)
6. hifiasm (0.19.6-r595 or higher)
7. samtools (Version: 1.11)
8. bedtools (v2.26.0)
9. seqkit (Version: 2.5.1)
10. BWA (0.7.17-r1188)
11. pilon-1.24.jar (version 1.24)
12. spades.py (v4.2.0)

## Installation
Run the following command to install RAGA and its dependencies.
```
(1) Download RAGA from GitHub
1. git clone https://github.com/wzxie/RAGA.git
2. chmod 755 /path/to/RAGA/bin/*
3. export PATH=/path/to/RAGA/bin/:$PATH

(2) Download dependencies
a. install with mamda
mamba create -y -n RAGA
mamba activate RAGA
mamba install bioconda::minimap2 racon mummer4 hifiasm samtools bedtools seqkit pysam bwa spades -y
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
git clone https://github.com/malonge/RagTag.git
cd RagTag
chmod +x *
export PATH=/path/to/RagTag:$PATH
# Add RagTag to the PATH.

b. install with source
# minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# racon
git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
# ragtag
git clone https://github.com/malonge/RagTag.git
# mummer
git clone https://github.com/mummer4/mummer.git
cd mummer
./configure --prefix=/path/to/installation
make
make install
# hifiasm
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
# samtools
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
make install
# bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
# seqkit
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_arm64.tar.gz
tar -zxvf *.tar.gz
# pilon-1.24.jar
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
# Add all the software to the PATH.
```
## Example
### 1. Download example data
```
wget https://download.cncb.ac.cn/gsa2/CRA008584/CRR591673/CRR591673.fastq.gz
wget https://cbi.gxu.edu.cn/downloadFiles/rpzhao/RAGA-test/Col-CEN_v1.2.fasta.gz
gzip -d CRR591673.fastq.gz
gzip -d Col-CEN_v1.2.fasta.gz
```
### 2. Run RAGA
```
RAGA.sh -r Col-CEN_v1.2.fasta -c CRR591673.fastq -t 8 &> RAGA.log
```
## Usage
### Quick start
```
Usage: RAGA.sh [-r reference genome] [-c target PacBio HiFi reads] [options]

Options:
    Input/Output:
    -r          reference genome
    -c          target PacBio HiFi reads
    -hic1       target Hi-C_r1 reads
    -hic2       target Hi-C_r2 reads

    Assembly type:
    -homo       assemble inbred/homozygous genomes
    -hetero     assemble heterozygous genomes
    -haplotype  generates a pair of haplotype-resolved assemblies with paired-end Hi-C reads

    Polish:
    -n INT      number of Polishing Rounds [>=3], default 3

    Filter:
    -i FLOAT    set the minimum alignment identity [0, 100], default 99
    -l INT      set the minimum alignment length, default 20,000
    -p FLOAT    extract the target PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
    -P FLOAT    extract the target longAlt read which aligns length is >= *% of its own length [0-1), default 0.5

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```
### Proceed in a stepwise fashion
Assembly may involve multiple types of reads, multiple reference genomes, closely related species as the reference genome, or non-haploid genomes.

(i) Initial_assembly.
```
# Select the appropriate de novo assembler for assembly.
```
(ii) Generate alternative long reads.

A. Homologous species as a reference.
```
Usage: RAGA-same.sh [-r reference genome] [-q target assembly] [-c target PacBio HiFi reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -q          target assembly
    -c          target PacBio HiFi reads
    -o          output directory

    Polish:
    -n INT      number of Polishing Rounds [>=3], default 3

    Filter:
    -i FLOAT    set the minimum alignment identity [0, 100], default 99
    -l INT      set the minimum alignment length, default 20,000
    -p FLOAT    extract the target PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
    -P FLOAT    extract the target longAlt read which aligns length is >= *% of its own length [0-1), default 0.5

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```
B. Closely related species as a reference.
```
Usage: RAGA-diff.sh [-r reference genome] [-c target PacBio HiFi reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -c          target PacBio HiFi reads
    -o          output directory

    Polish:
    -n INT      number of Polishing Rounds [>=10], default 10

    Filter:
    -l INT      set the minimum alignment length, default 10,000

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```
(iii) Optimized_assembly.
```
# Incorporate alternative long reads into the original reads and reassemble.
```
## Inputs and Outputs
### Input files
The reference genome sequence file
```
>Chr01
ATCGATCGATCGATCGATCGATCG...
```
The target assembly sequence file
```
>Contig1
ATCGATCGATCGATCGATCGATCG...
>Contig2
ATCGATCGATCGATCGATCGATCG...
```
The target PacBio HiFi sequence file
```
@*/ccs1
ATCGATCGATCGATCGATCGATCG...
+
~{m~pI~~e~P~c~~o~[o~q_~~...
@*/ccs2
ATCGATCGATCGATCGATCGATCG...
+
=<G~~n~~~~V~~~~~~~}~~~~~...
```
### Output files
```
* Initial_assembly/initial.fa        # de novo assembly
* Alternative_reads/longAlt_tgt.fa   # Alternative long reads by RAGA
* Optimized_assembly/optimized.fa    # Optimized assembly
* The 'gapA_area.svg'                # Visualization of GAP area
* The 'longAlt_tgt_lenDis.svg'       # statistics of the 'longAlt_tgt.fa'
```

## NGS testing
### 1. Download example data
```
wget https://cbi.gxu.edu.cn/downloadFiles/rpzhao/RAGA-test/NGS-test.zip
```

### 2. Run RAGA-ngs.sh
```
RAGA-ngs.sh -r GCA_949124345.1.fa -read1 read1.fq -read2 read2.fq -t 8 &> RAGA-ngs.log
```

### Output files
```
* longAlt_ngs/longAlt_ngs.fa         # Alternative long reads by RAGA-ngs
```

## Note
* RAGA aligns the reference with the target contigs to identify alignment blocks near the gaps when the input is the reference of the homologous species. Hence, RAGA's output will be more reliable if the input target assembly is of higher quality.
* In the absence of a high-quality reference for the target species, users can input PacBio HiFi reads from the target species and the reference FASTA file of closely related species directly, eliminating the need for target assembly.

## Citing RAGA
Ru-Peng Zhao#, Yu-Hong Luo#, Wen-Zhao Xie#, Zu-Wen Zhou, Yong-Qing Qian, Si-Long Yuan, Dong-Ao Li, Jiana Li, Kun Lu, Xingtan Zhang, Jia-Ming Song*, Ling-Ling Chen*, RAGA: a reference-assisted genome assembly tool for efficient population-scale assembly, Horticulture Research, 2025;, uhaf207, https://doi.org/10.1093/hr/uhaf207

## Contact
We hope this tools could be helpful for the groups which focused on plants genome assembly, you can use the GitHub page to report issues or email us with any suggestions.
* Zhao rupeng:    2247290650@qq.com
* Luo yuhong:     luoyuhong0720@163.com
* Xie wenzhao:    wzxie@hebtu.edu.cn
* Song jiaming:   jmsong@swu.edu.cn
* Chen lingling:  llchen@gxu.edu.cn
