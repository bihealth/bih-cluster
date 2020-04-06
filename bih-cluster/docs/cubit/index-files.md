# Precomputed Index Files

Index files for

* BWA version 0.7.12 and 0.7.15,
* bowtie2 version 2.2.5 and
* STAR version 2.4.1d

have been precomputed. The index corresponding to each genome is stored in the
following directory structure with the above mentioned reference genomes as
subfolders (listed here only for `Bowtie/1.1.2`, same subfolders for the
remaining programs):

```
static_data/precomputed
├── Bowtie
│   └── 1.1.2
│       ├── danRer10
│       ├── dm6
│       ├── ecoli
│       ├── GRCh37
│       ├── GRCh38
│       ├── GRCm38
│       ├── hg18
│       ├── hg19
│       ├── hg38
│       ├── mm10
│       ├── mm9
│       ├── NCBIM37
│       ├── phix
│       ├── sacCer3
│       ├── UniVec
│       └── UniVec_Core
├── Bowtie2
│   └── 2.2.5
│       └── [see Bowtie/1.1.2]
├── BWA
│   ├── 0.7.12
│   │   └── [see Bowtie/1.1.2]
│   └── 0.7.15
│       └── [see Bowtie/1.1.2]
└── STAR
    └── 2.4.1d
        └── default
            └── [see Bowtie/1.1.2]
```
