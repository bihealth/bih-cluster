# First Steps: Episode 1

|Episode|Topic|
|---|---|
| 0 | [How can I install the tools?](episode-0.md) |
| **1** | **How can I use the static data?** |
| 2 | [How can I distribute my jobs on the cluster (Slurm)?](episode-2.md) |
| 3 | [How can I organize my jobs with Snakemake?](episode-3.md) |
| 4 | [How can I combine Snakemake and Slurm?](episode-4.md) |

This is part one of the "First Steps" BIH Cluster Tutorial.
Here we will build a small pipeline with alignment and variant calling.
The premise is that you have the tools installed as described in [Episode 0](episode-0.md). For this episode, please make sure that you
are on a compute node. As a reminder, the command to access a compute node with the required resources is

```
$ srun --time 7-00 --mem=8G --cpus-per-task=8 --pty bash -i
```

## Tutorial Input Files

We will provide you with some example FASTQ files, but you can use your own if you like.
You can find the data here:

- `/data/cephfs-1/work/projects/cubit/tutorial/input/test_R1.fq.gz`
- `/data/cephfs-1/work/projects/cubit/tutorial/input/test_R2.fq.gz`

## Creating a Project Directory

First, you should create a folder where the output of this tutorial will go.
It would be good to have it in your `work` directory in `/data/cephfs-1/home/users/$USER`, because it is faster and there is more space available.

```terminal
(first-steps) $ mkdir -p /data/cephfs-1/home/users/$USER/work/tutorial/episode1
(first-steps) $ pushd /data/cephfs-1/home/users/$USER/work/tutorial/episode1
```

!!! important "Quotas / File System limits"

    - Note well that you have a quota of 1 GB in your home directory at `/data/cephfs-1/home/users/$USER`.
      The reason for this is that nightly snapshots and backups are created for this directory which are precious resources.
    - This limit does not apply to your work directory at `/data/cephfs-1/home/users/$USER/work`.
      The limits are much higher here but no snapshots or backups are available.
    - There is no limit on your scratch directory at `/data/cephfs-1/home/users/$USER/scratch`.
      However, **files placed here are automatically removed after 2 weeks.**
      This is only appropriate for files during download or temporary files.

## Creating a Directory for Temporary Files

In general it is advisable to have a proper temporary directory available.
You can create one in your `~/scratch` folder and make it available to the system.

```terminal
(first-steps) $ export TMPDIR=/data/cephfs-1/home/users/$USER/scratch/tmp
(first-steps) $ mkdir -p $TMPDIR
```

## Using the Cubit Static Data

The static data is located in `/data/cephfs-1/work/projects/cubit/current/static_data`.
For our small example, the required reference genome and index can be found at:

- `/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta`
- `/data/cephfs-1/work/projects/cubit/current/static_data/precomputed/BWA/0.7.17/GRCh37/g1k_phase1/human_g1k_v37.fasta`

## Aligning the Reads

Let's align our data:

```terminal
(first-steps) $ bwa mem -t 8 \
    -R "@RG\tID:FLOWCELL.LANE\tPL:ILLUMINA\tLB:test\tSM:PA01" \
    /data/cephfs-1/work/projects/cubit/current/static_data/precomputed/BWA/0.7.17/GRCh37/g1k_phase1/human_g1k_v37.fasta \
    /data/cephfs-1/work/projects/cubit/tutorial/input/test_R1.fq.gz \
    /data/cephfs-1/work/projects/cubit/tutorial/input/test_R2.fq.gz \
| samtools view -b \
| samtools sort -O BAM -T $TMPDIR -o aln.bam

(first-steps) $ samtools index aln.bam
```

## Perform Structural Variant Calling

And do the structural variant calling:

```terminal
(first-steps) $ delly call \
    -g /data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta \
    aln.bam
```

Note that delly will not find any variants.

## Small Variant Calling (SNV, indel)

And now for the SNP calling (this step will take ~ 20 minutes):

```terminal
(first-steps) $ gatk HaplotypeCaller \
    -R /data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta \
    -I aln.bam \
    -ploidy 2 \
    -O test.GATK.vcf
```

## Outlook: More Programs and Static Data

So this is it!
We used the tools that we installed previously, accessed the reference data and ran a simple alignment and variant calling pipeline.
You can access a list of all static data through this wiki, follow this link to the [Static Data](../cubit/index.md).
You can also have a peek via:

```terminal
(first-steps) $ tree -L 3 /data/cephfs-1/work/projects/cubit/current/static_data | less
```
