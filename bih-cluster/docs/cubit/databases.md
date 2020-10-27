# Databases

The file formats in the `static_data/db` folder are mostly `.vcf` or `.bed`
files. We provide the following databases:

| Database | Version | Reference genome |
|----------|---------|------------------|
| COSMIC   | v72     | GRCh37           |
| dbNSFP   | 2.9     | GRCh37/hg19      |
| dbSNP    | b128    | mm9              |
| dbSNP    | b128    | NCBIM37          |
| dbSNP    | b142    | GRCh37           |
| dbSNP    | b144    | GRCh38           |
| dbSNP    | b147    | GRCh37           |
| dbSNP    | b147    | GRCh38           |
| dbSNP    | b150    | GRCh37           |
| dbSNP    | b150    | GRCh38           |
| DGV      | 2015-07-23 | GRCh37        |
| ExAC     | release0.3 | GRCh37/hg19   |
| ExAC     | release0.3.1 | GRCh37/hg19 |
| giab     | NA12878_HG001/NISTv2.19 | GRCh37 |
| goldenpath | variable | GRCh37        |
| goldenpath | variable | hg19          |
| goldenpath | variable | mm9           |
| goldenpath | variable | NCBIM37       |
| SangerMousegenomesProject | REL-1211-SNPs_Indels | mm9 |
| SangerMousegenomesProject | REL-1211-SNPs_Indels | NCBIM37 |
| UK10K cohort | REL-2012-06-02 | GRCh37 |

The directory structure is as follows:

```
static_data/db
├── COSMIC
│   └── v72
│       └── GRCh37
├── dbNSFP
│   └── 2.9
├── dbSNP
│   ├── b128
│   │   ├── mm9
│   │   └── NCBIM37
│   ├── b142
│   │   └── GRCh37
│   ├── b144
│   │   └── GRCh38
│   └── b147
│       ├── GRCh37
│       └── GRCh38
├── DGV
│   └── 2015-07-23
│       └── GRCh37
├── ExAC
│   ├── release0.3
│   └── release0.3.1
├── giab
│   └── NA12878_HG001
│       └── NISTv2.19
├── goldenpath
│   └── variable
│       ├── GRCh37
│       ├── hg19
│       ├── mm9
│       └── NCBIM37
├── SangerMouseGenomesProject
│   └── REL-1211-SNPs_Indels
│       ├── mm9
│       └── NCBIM37
└── UK10K_cohort
    └── REL-2012-06-02
```
