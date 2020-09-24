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

While file paths are highlighted like this: `/fast/projects/cubit/current`.

## Instant Gratification

After connecting to the cluster, you are located on a login node.
To get to your first compute node, type `srun --time 7-00 --pty bash -i` which will launch an interactive Bash session on a free remote node running up to 7 days.
Typing `exit` will you bring back to the login node.

```terminal
$ srun -p long --time 7-00 --pty bash -i
med0107 $ exit
$
```

See?
That was easy!

## Preparation

In preparation for our first steps tutorial series, we would like you to install the software for this tutorial.
In general the users on the cluster will manage their own software with the help of conda.
If you haven't done so so far, please [follow the instructions in installing conda](../best-practice/software-installation-with-conda.md) first.
The only premise is that you are able to [log into the cluster](../connecting/configure-ssh/linux.md).
Make also sure that you are logged in to a computation node using ``.

Now we will create a new environment, so as to not interfere
with your current or planned software stack, and install into it all the
software that we need during the tutorial. Run the following commands:

```terminal
$ conda create -n first-steps python=3 snakemake drmaa bwa delly samtools gatk4
$ conda activate first-steps
(first-steps) $
```

As you can see, we also installed a tool called DRMAA (which is in fact a software library). DRMAA is an API that provides more stable job distribution on the cluster than the native SGE implementation. We will use this during the tutorial and we recommend it for your every day usage. If you are interested, we have a [wiki page](../slurm/snakemake.md#snakemake-and-slurm) about this library and how to use it. For now there is no need for you to set up anything more.
