# First Steps: Episode 0

|Episode|Topic|
|---|---|
| **0** | **How can I install the tools?** |
| 1 | [How can I use the static data?](episode-1.md) |
| 2 | [How can I distribute my jobs on the cluster (Slurm)?](episode-2.md) |
| 3 | [How can I organize my jobs with Snakemake?](episode-3.md) |
| 4 | [How can I combine Snakemake and Slurm?](episode-4.md) |

## Prerequisites

This tutorial assumes familiarity with Linux/Unix operating systems.
It also assumes that you have already connected to the cluster.
We have collected some links to [tutorials and manuals on the internet](../misc/external-resources.md).

## Legend

Before we start with our first steps tutorial, we would like to
introduce the following convention that we use throughout the series:

```terminal
$ Commands are prefixed with a little dollar sign
```

While file paths are highlighted like this: `/data/cephfs-1/work/projects/cubit/current`.

## Instant Gratification

After connecting to the cluster, you are located on a login node.
To get to your first compute node, type `srun --time 7-00 --mem=8G --cpus-per-task=8 --pty bash -i` which will launch an interactive Bash session on a free remote node running up to 7 days, enabling you to use 8 cores and 8 Gb memory. Typing `exit` will you bring back to the login node.

```terminal
hpc-login-1$ srun -p long --time 7-00 --mem=8G --cpus-per-task=8 --pty bash -i
hpc-cpu-1$ exit
$
```

See?
That was easy!

## Preparation

In preparation for our first steps tutorial series, we would like you to install the software for this tutorial.
In general the users on the cluster will manage their own software with the help of conda.
If you haven't done so so far, please [follow the instructions in installing conda](../best-practice/software-installation-with-conda.md) first.
The only premise is that you are able to [log into the cluster](../connecting/advanced-ssh/unix.md).
Make also sure that you are logged in to a computation node using `srun -p medium --time 1-00 --mem=4G --cpus-per-task=1 --pty bash -i`.

Now we will create a new environment, so as to not interfere
with your current or planned software stack, and install into it all the
software that we need during the tutorial. Run the following commands:

```terminal
$ conda create -n first-steps python=3 snakemake bwa delly samtools gatk4
$ conda activate first-steps
(first-steps) $
```
