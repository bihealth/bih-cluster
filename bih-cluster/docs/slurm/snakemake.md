# Snakemake with Slurm

This page describes how to use Snakemake with Slurm.

## Prerequisites

- This assumes that you have Miniconda properly setup with Bioconda.
- Also it assumes that you have already activated the Miniconda base environment with `source miniconda/bin/activate`.

## Environment Setup

We first create a new environment `snakemake-slurm` and activate it.
We need the `snakemake` package for this.

```bash
host:~$ conda create -y -n snakemake-slurm snakemake
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
host:snake-slurm$ cat >Snakefile <<"EOF"
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

You have two options:

1. Simply use `snakemake --profile=cubi-v1` and the Snakemake resource configuration as shown below. **STRONGLY PREFERRED**
2. Use the `snakemake --cluster='sbatch ...'` command.

Note that we sneaked in a `sleep 1m`? In a second terminal session, we can see that the job has been submitted to SLURM indeed.

```bash
host:~$ squeue  -u holtgrem_c
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
               325     debug snakejob holtgrem  R       0:47      1 med0127
```

## Threads & Resources

The `cubi-v1` profile (stored in `/etc/xdg/snakemake/cubi-v1` on all cluster nodes) supports the following specification in your Snakemake rule:

* `threads`: the number of threads to execute the job on
* memory in a syntax understood by Slurm, EITHER
    * `resources.mem`/`resources.mem_mb`: the memory to allocate **for the whole job**, OR 
    * `resources.mem_per_thread`: the memory to allocate **for each thread**.
* `resources.time`: the running time of the rule, in a syntax supported by Slurm, e.g. `HH:MM:SS` or `D-HH:MM:SS`
* `resources.partition`: the partition to submit your job into (Slurm will pick a fitting partition for you by default)
* `resources.nodes`: the number of nodes to schedule your job on (defaults to `1` and you will want to keep that value unless you want to use MPI)

You will need Snakemake >=7.0.2 for this.

So for example:

```python
rule myrule:
    threads: 1
    resources:
        mem='8G',
        time='04:00:00',
    input: # ...
    output: # ...
    shell: # ...
```

You can combine this with Snakemake [resource callables](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=resources#resources), of course:

```python
def myrule_mem(wildcards, attempt):
    mem = 2 * attempt
    return '%dG' % mem

rule snps:
    threads: 1
    resources:
        mem=myrule_mem,
        time='04:00:00',
    input: # ...
    output: # ...
    shell: # ...
```

Here is how to call Snakemake:

```bash
# snakemake --profile=cubi-v1 -j1
```

Note that right now you will need an unreleased version of Snakemake v7:

```bash
# pip install git+https://github.com/snakemake/snakemake.git@6d19f7e5585d9e9ee93c66222b7883ef09c1bc1d
```
