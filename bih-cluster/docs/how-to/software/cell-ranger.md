# How-To: Run CellRanger

# what is Cell Ranger?
from the official [website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger):
"Cell Ranger is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis"

# installation

requires registration before download from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

to unpack Cell Ranger, its dependencies and the `cellranger` script:

```
cd /data/cephfs-1/home/users/$USER/work
mv /path/to/cellranger-3.0.2.tar.gz .
tar -xzvf cellranger-3.0.2.tar.gz
```

# reference data

will be provided in `/data/cephfs-1/work/projects/cubit/current/static_data/app_support/cellranger`

# cluster support SLURM

add a file `slurm.template` to `/data/cephfs-1/home/users/$USER/work/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/sge.template` with the following contents:

```
#!/usr/bin/env bash
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# =============================================================================
# Setup Instructions
# =============================================================================
#
# 1. Add any other necessary Slurm arguments such as partition (-p) or account
#    (-A). If your system requires a walltime (-t), 24 hours (24:00:00) is
#    sufficient.  We recommend you do not remove any arguments below or Martian
#    may not run properly.
#
# 2. Change filename of slurm.template.example to slurm.template.
#
# =============================================================================
# Template
# =============================================================================
#
#SBATCH -J __MRO_JOB_NAME__
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=__MRO_THREADS__
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --partition=medium
#SBATCH --time=24:00:00
### Alternatively: --ntasks=1 --cpus-per-task=__MRO_THREADS__
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=__MRO_MEM_GB__G
#SBATCH -o __MRO_STDOUT__
#SBATCH -e __MRO_STDERR__

__MRO_CMD__
```

**note**: on newer cellranger version, `slurm.template` needs to go to `/data/cephfs-1/home/users/$USER/work/cellranger-XX/external/martian/jobmanagers/`

# demultiplexing

if that hasn't been done yet, you can use `cellranger mkfastq` (details to be added)

# run the pipeline (`count`)

create a script `run_cellranger.sh` with these contents (consult the [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) for help:

```
#!/bin/bash

/data/cephfs-1/home/users/$USER/work/cellranger-3.0.2/cellranger count \
  --id=sample_id \
  --transcriptome=/data/cephfs-1/work/projects/cubit/current/static_data/app_support/cellranger/refdata-cellranger-${species}-3.0.0\
  --fastqs=/path/to/fastqs \
  --sample=sample_name \
  --expect-cells=n_cells \
  --jobmode=slurm \
  --maxjobs=100 \
  --jobinterval=1000
```

and then submit the job via

```
sbatch --ntasks=1 --mem-per-cpu=4G --time=8:00:00 -p medium -o cellranger.log run_cellranger.sh
```

# cluster support SGE (outdated)

add a file `sge.template` to `/data/cephfs-1/home/users/$USER/work/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/sge.template` with the following contents:

```
# =============================================================================
# Template
# =============================================================================
#
#$ -N __MRO_JOB_NAME__
#$ -V
#$ -pe smp __MRO_THREADS__
#$ -cwd
#$ -P medium
#$ -o __MRO_STDOUT__
#$ -e __MRO_STDERR__
#$ -l h_vmem=__MRO_MEM_GB_PER_THREAD__G
#$ -l h_rt=08:00:00

#$ -m a
#$ -M user@email.com

__MRO_CMD__
```

and submit the job via

```
 qsub -cwd -V -pe smp 1 -l h_vmem=8G -l h_rt=24:00:00 -P medium -m a -j y run_cellranger.sh
```
