# Snakemake with Slurm

This page describes how to use Snakemake with Slurm.

## Prerequisites

- This assumes that you have Miniconda properly setup with Bioconda.
- Also it assumes that you have already activated the Miniconda base environment with `source miniconda/bin/activate`.

## Environment Setup

We first create a new environment `snakemake-slurm` and activate it.
We need the `snakemake` package and the `drmaa` Python library for interfacing with the SLURM scheduler via the DRMAA interface.

```bash
host:~$ conda create -y -n snakemake-slurm snakemake drmaa
[...]
#
# To activate this environment, use
#
#     $ conda activate snakemake-slurm
#
# To deactivate an active environment, use
#
#     $ conda deactivate
host:~$ conda activate snakemake-slurm
(snakemake-slurm) host:~$
```

## Snakemake Workflow Setup

We create a workflow and ensure that it works properly with multi-threaded Snakemake (no cluster submission here!)

```bash
host:~$ mkdir -p snake-slurm
host:~$ cd snake-slurm
host:snake-slurm$ cat <<"EOF"
rule default:
    input: "the-result.txt"

rule mkresult:
    output: "the-result.txt"
    shell: r"sleep 1m; touch the-result.txt"
EOF
host:snake-slurm$ snakemake --cores=1
[...]
host:snake-slurm$ ls
Snakefile  the-result.txt
host:snake-slurm$ rm the-result.txt
```

## Snakemake and :tada: Slurm

It's really simple:

0. unset the `DRMAA_LIBRARY_PATH` variable (that might be set and point to the SGE DRMAA library); you could set it to `/usr/lib64/libdrmaa.so` but that would not be necessary, consider all mentions of it from `~/.bashrc`.
1. Use `snakemake --drmaa " [params]"` and use SLURM syntax here.

```bash
host:snake-slurm$ unset DRMAA_LIBRARY_PATH
host:snake-slurm$snakemake --drmaa " -t 05:00" --jobs 2
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 2
Job counts:
        count   jobs
        1       default
        1       mkresult
        2

[Wed Mar 25 23:06:01 2020]
rule mkresult:
    output: the-result.txt
    jobid: 1

Submitted DRMAA job 1 with external jobid 325.
[Wed Mar 25 23:07:11 2020]
Finished job 1.
1 of 2 steps (50%) done

[Wed Mar 25 23:07:11 2020]
localrule default:
    input: the-result.txt
    jobid: 0

[Wed Mar 25 23:07:11 2020]
Finished job 0.
2 of 2 steps (100%) done
Complete log: /fast/home/users/holtgrem_c/snake-slurm/.snakemake/log/2020-03-25T230601.353735.snakemake.log
```

Note that we sneaked in a `sleep 1m`? In a second terminal session, we can see that the job has been submitted to SLURM indeed.

```bash
host:~$ squeue  -u holtgrem_c
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
               325     debug snakejob holtgrem  R       0:47      1 med0127
```

## Limitations

The DRMAA interface to Slurm has a few limitations:

- memory has to be given as an integer in the unit Megabytes (1024 * 1024 bytes) and as a plain number without unit,
- running time has to be given as `hh:mm`.

A full list of supported parameters can be found [in the officical documentation](http://apps.man.poznan.pl/trac/slurm-drmaa#Nativespecification).

... that's all, folks!