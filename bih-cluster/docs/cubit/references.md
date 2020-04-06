# Reference Sequences

## NCBI mouse reference genome assemblies

We provide the NCBI mouse reference assembly used by the [Sanger Mouse Genomics
group](http://www.sanger.ac.uk/resources/mouse/genomes) for NCBIM37 and GRCm38.
This is a reliable source where the appropriate contigs have already been
selected by experts. NCBIM37 is annotated with Ensembl release 67 and GRCm38
with Ensembl release 68.

## UCSC mouse reference genome assemblies

The assembly sequence is in one file per chromosome and is available for mm9
and mm10. We concatenated all the chromosome files to one final fasta file for
each genome assembly.

## NCBI human reference genome assemblies

* **GRCh37:** We provide the version used by the 1000genomes project as it is
widely used and recommended. The chromosomes and contigs are already
concatenated.
    - **g1k_phase1/hs37:** This reference sequence contains the autosomal and both
    sex chromosomes, an updated mitochondrial chromosome as well as
    "non-chromosomal supercontigs". The
    [README](http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/README.human_g1k_v37.fasta.txt)
    explains the method of construction.
    - **g1k_phase2/hs37d5:** In addition to these sequences the phase 2 reference
    sequence contains the herpes virus genome and decoy sequences for improving
    SNP calling.
* **GRCh38:** The GRCh38 assembly offers an "analysis set" that was created to
accommodate next generation sequencing read alignment pipelines. We provide the
three analysis sets from the NCBI.
    - **hs38/no_alt_analysis_set:** The chromosomes, mitochondrial genome,
    unlocalized scaffolds, unplaced scaffolds and the Epstein-Barr virus
    sequence which has been added as a decoy to attract contamination in
    samples.
    - **hs38a/full_analysis_set:** the alternate locus scaffolds in addition to all
    the sequences present in the no_alt_analysis_set.
    - **hs38DH/full_plus_hs38d1_analysis_set:** contains the human decoy sequences
    from hs38d1 in addition to all the sequences present in the full_analysis
    set.  More detailed information is available in the
    [README](https://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/README).

## UCSC human reference genome assemblies

The assembly sequence is in one file per chromosome is available for hg18, hg19
and hg38. We concatenated all the chromosome files to one  final fasta file for
each genome assembly. Additionally, in the subfolder `chromosomes` we keep the
chromosome fasta files separately for hg18 and hg19.

## Other reference genomes

* **danRer10:** UCSC/GRC zebrafish build 10
* **dm6:** UCSC/GRC Drosophila melanogaster build 6
* **ecoli:**
    - **GCA_000005845.2_ASM584v2:** Genbank Escherichia coli K-12 subst. MG1655 genome
* **genomemedley:**
    - **1:** Concatenated genome of hg19, dm6, mm10; Chromosomes are tagged with corresponding organism
* **PhiX:** Control genome that is used by Illumina for sequencing runs
* **sacCer3:** UCSC's Saccharomyces cerevisiae genome build 3
* **UniVec:**
    - **9:** NCBI's non redundant reference of vector sequences, adapters, linkers and primers commonly used in the process of cloning cDNA or genomic DNA (build 9)
* **UniVec_Core**
    - **9:** A subset of UniVec build 9

The following directory structure indicates the available genomes. Where there
isn't a name for the data set, either the source (e.g. sanger - from the Sanger
Mouse Genomes project) or the download date is used to name the sub-directory.

```
static_data/reference
├── danRer10
│   └── ucsc
├── dm6
│   └── ucsc
├── ecoli
│   └── GCA_000005845.2_ASM584v2
├── genomemedley
│   └── 1
├── GRCh37
│   ├── g1k_phase1
│   ├── g1k_phase2
│   ├── hs37
│   └── hs37d5
├── GRCh38
│   ├── hs38
│   ├── hs38a
│   └── hs38DH
├── GRCm38
│   └── sanger
├── hg18
│   └── ucsc
├── hg19
│   └── ucsc
├── hg38
│   └── ucsc
├── mm10
│   └── ucsc
├── mm9
│   └── ucsc
├── NCBIM37
│   └── sanger
├── phix
│   └── illumina
├── sacCer3
│   └── ucsc
├── UniVec
│   └── 9
└── UniVec_Core
    └── 9
```
