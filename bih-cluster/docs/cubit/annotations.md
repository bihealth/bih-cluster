# Annotation Data

The following Ensembl and GENCODE versions corresponding to the indicated
reference genomes will be made available on the cluster.

| Database | Version | Reference Genome |
|---        |---      |---      |
|Ensembl | 65 | NCBIM37 (Ensembl release corresponding to GENCODE M1) |
|Ensembl | 67 | NCBIM37 (Ensembl release for sanger mouse genome assembly) |
|Ensembl | 68 | GRCm38 (Ensembl release for sanger mouse genome assembly) |
|Ensembl | 74 | GRCh37 (Ensembl release for GENCODE 19) |
|Ensembl | 75 | GRCh37 (Latest release for GRCh37) |
|Ensembl | 79 | GRCh38 (Ensembl release for GENCODE 22) |
|Ensembl | 80 | GRCh38 (Ensembl release corresponding to GENCODE 22) |
|Ensembl | 80 | GRCm38 (Ensembl release corresponding to GENCODE M1) |
|GENCODE | M1 | NCBIM37 (No gff3 file) |
|GENCODE | M5 | GRCm38 |
|GENCODE | 19 | current for GRCh37 |
|GENCODE | 22 | current for GRCh38 |

The annotation files associated with the indicated genomes can be accessed in
the following directories:

```
static_data/annotation
├── ENSEMBL
│   ├── 65
│   │   └── NCBIM37
│   ├── 67
│   │   └── NCBIM37
│   ├── 68
│   │   └── GRCm38
│   ├── 74
│   │   └── GRCh37
│   ├── 75
│   │   └── GRCh37
│   ├── 79
│   │   └── GRCh38
│   └── 80
│       ├── GRCh38
│       └── GRCm38
└── GENCODE
    ├── 19
    │   └── GRCh37
    ├── 22
    │   └── GRCh38
    ├── M1
    │   └── NCBIM37
    └── M5
        └── GRCm38
```
