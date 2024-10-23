# Snakemake with Slurm

This page describes how to use Snakemake with Slurm.

## Prerequisites

- This assumes that you have Miniforge properly setup with Bioconda.
- Also it assumes that you have already activated the Miniforge base environment with `source miniforge/bin/activate`.

## Environment Setup

We first create a new environment `snakemake-slurm` and activate it.
We need the `snakemake` package for this.
For snakemake 8, we additionally need [`snakemake-executor-plugin-slurm`](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).

```bash
host:~$ conda create -y -n snakemake-slurm 'snakemake>=8.24.1' snakemake-executor-plugin-slurm
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
    shell: r"sleep 2m; touch the-result.txt"
EOF
host:snake-slurm$ snakemake --cores=1
[...]
host:snake-slurm$ ls
Snakefile  the-result.txt
host:snake-slurm$ rm the-result.txt
```

## Snakemake and :tada: Slurm

Simply use `snakemake --profile=cubi-v1` and the Snakemake resource configuration as shown below.
This works for both snakemake version 7 *and* 8.

(For version 7, you can also use `snakemake --cluster='sbatch ...'` command instead, but this is discouraged.)

Note that we sneaked in a `sleep 2m`? In a second terminal session, we can see that the job has been submitted to SLURM indeed.

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

For more information on resources, see [snakemake resources](https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#standard-resources) and [slurm resources](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html#advanced-resource-specifications).

Here is how to call Snakemake:

```bash
# snakemake --profile=cubi-v1 --jobs 1
```
When using conda, additionally use `--sdm conda` (snakemake 8) or `--with-conda` (snakemake 7).

To set rule-specific resources:

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

## Custom logging directory

By default, slurm will write log files into the working directory of snakemake, which will look like `slurm-$jobid.out`.

To change this behaviour, the environment variable `SBATCH_DEFAULTS` can be set to re-route the `--output` parameter. If you want to write your files into `slurm_logs` with a filename pattern of `$name-$jobid` for instance, consider the following snippet for your submission script:

```bash
#!/bin/bash
#
#SBATCH --job-name=snakemake_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:10:00
#SBATCH --mem-per-cpu=300M
#SBATCH --output=slurm_logs/%x-%j.log

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date
srun snakemake --use-conda -j1 --profile=cubi-v1
date

```

The name of the snakemake slurm job will be `snakemake_main_job`, the name of the jobs spawned from it will be called after the rule name in the Snakefile.

