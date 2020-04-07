# How-To: Run CellRanger

!!! todo "Not yet updated to Slurm"

    TODO: This still needs to be updated to Slurm.

# what is Cell Ranger?
from the official [website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger):
"Cell Ranger is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis"

# installation

requires registration before download from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

to unpack Cell Ranger, its dependencies and the `cellranger` script:

```
cd /fast/users/$USER/scratch
mv /path/to/cellranger-3.0.2.tar.gz .
tar -xzvf cellranger-3.0.2.tar.gz
```

# reference data

will be provided in `/fast/projects/cubit/current/static_data/app_support/cellranger`

# cluster support

add a file `sge.template` to `/fast/users/$USER/scratch/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/sge.template` with the following contents:

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

# demultiplexing

if that hasn't been done yet, you can use `cellranger mkfastq` (details to be added)

# run the pipeline (`count`)

create a script `run_cellranger.sh` with these contents (consult the [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) for help:

```
#!/bin/bash

/fast/users/$USER/scratch/cellranger-3.0.2/cellranger count \
  --id=sample_id \
  --transcriptome=/fast/projects/cubit/current/static_data/app_support/cellranger/refdata-cellranger-${species}-3.0.0\
  --fastqs=/path/to/fastqs \
  --sample=sample_name \
  --expect-cells=n_cells \
  --jobmode=sge \
  --maxjobs=100 \
  --jobinterval=1000
```

and then submit the job via

```
qsub -cwd -V -pe smp 1 -l h_vmem=4G -l h_rt=8:00:00 -P medium -m a -o cellranger.log run_cellranger.sh
```