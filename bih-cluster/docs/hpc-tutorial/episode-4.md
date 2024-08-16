# First Steps: Episode 4

|Episode|Topic|
|---|---|
| 0 | [How can I install the tools?](episode-0.md) |
| 1 | [How can I use the static data?](episode-1.md) |
| 2 | [How can I distribute my jobs on the cluster (Slurm)?](episode-2.md) |
| 3 | [How can I organize my jobs with Snakemake?](episode-3.md) |
| **4** | **How can I combine Snakemake and Slurm?** |

In the last episodes we learned about distributing a job among the cluster nodes using `sbatch` and
how to automate and parallelize our pipeline with Snakemake. We are lucky that those two
powerful commands can be combined. What is the result? You will have an automated pipeline
with Snakemake that uses `sbatch` to distribute jobs among the cluster nodes instead of
running only the same node.

The best thing is that we can reuse our `Snakefile` as it is and just write a wrapper script
to call Snakemake. We run the script and the magic will start.

First, create a new folder for this episode:

```terminal
(first-steps) $ mkdir -p /data/cephfs-1/home/users/${USER}/work/tutorial/episode4/logs
(first-steps) $ pushd /data/cephfs-1/home/users/${USER}/work/tutorial/episode4
```

And copy the wrapper script to this folder as well as the Snakefile (you can also reuse the one with the adjustments from the previous [episode](episode-3.md)):

```terminal
(first-steps) $ cp /data/cephfs-1/work/projects/cubit/tutorial/skeletons/submit_snakejob.sh .
(first-steps) $ cp /data/cephfs-1/work/projects/cubit/tutorial/skeletons/Snakefile .
(first-steps) $ chmod u+w submit_snakejob.sh Snakefile
```

The `Snakefile` is already known to you but let me explain the wrapper script `submit_snakejob.sh`:

```bash
#!/bin/bash

# Set a name for the job (-J or --job-name).
#SBATCH --job-name=tutorial

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=2

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the total memory. Units can be given in T|G|M|K.
#SBATCH --mem=1G

# Optionally, set the partition to be used (-p or --partition).
#SBATCH --partition=medium

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=30:00


export TMPDIR=/data/cephfs-1/home/users/${USER}/scratch/tmp
export LOGDIR=logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
mkdir -p $LOGDIR

eval "$($(which conda) shell.bash hook)"
conda activate first-steps

set -x

snakemake --profile=cubi-v1 -j 2 -k -p --restart-times=2
```

In the beginning you see the `#SBATCH` that introduces the parameters when you provide this script to `sbatch` as described in the [second episode](episode-2.md).
Please make sure that the `logs` folder exists before starting the run!
We then set and export the `TMPDIR` and `LOGDIR` variables.
Note that `LOGDIR` has a subfolder named `$SLURM_JOB_NAME-$SLURM_JOB_ID` that will be created for you.
Snakemake will store its logfiles for this very Snakemake run in this folder.
The next new thing is `set -x`. This simply prints to the terminal every command that is executed within the script.
This is useful for debugging.

Finally, the Snakemake call takes place.
With the `--profile` option we define that Snakemake uses the [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) at `/etc/xdg/snakemake/cubi-v1`.
The profile will take create appropriate calls to `sbatch` and interpret the following settings from your Snakemake rule:

* `threads`: the number of threads to execute the job on
* memory in megabytes or with a suffix of `k`, `M`, `G`, or `T`. You can specify EITHER
    * `resources.mem`/`resources.mem_mb`: the memory to allocate **for the whole job**, OR
    * `resources.mem_per_thread`: the memory to allocate **for each thread**.
* `resources.time`: the running time of the rule, in a syntax supported by Slurm, e.g. `HH:MM:SS` or `D-HH:MM:SS`
* `resources.partition`: the partition to submit your job into (Slurm will pick a fitting partition for you by default)
* `resources.nodes`: the number of nodes to schedule your job on (defaults to `1` and you will want to keep that value unless you want to use MPI)

The other options to `snakemake` have the meaning:

* `-j 2`: run at most two jobs at the same time
* `-k`: keep going even if a rule execution fails
* `-p`: print the executed shell commands
* `--restart-times=2`: restart failing jobs up to two times

It is now time to update your `Snakefile` such that it actually specifies the resources mentioned above:

```python
rule all:
    input:
        'snps/test.vcf',
        'structural_variants/test.vcf'

rule alignment:
    input:
        '/data/cephfs-1/work/projects/cubit/tutorial/input/{id}_R1.fq.gz',
        '/data/cephfs-1/work/projects/cubit/tutorial/input/{id}_R2.fq.gz',
    output:
        bam='alignment/{id}.bam',
        bai='alignment/{id}.bam.bai',
    threads: 8
    resources:
        mem='8G',
        time='12:00:00',
    shell:
        r"""
        export TMPDIR=/data/cephfs-1/home/users/${{USER}}/scratch/tmp
        mkdir -p ${{TMPDIR}}

        BWAREF=/data/cephfs-1/work/projects/cubit/current/static_data/precomputed/BWA/0.7.17/GRCh37/g1k_phase1/human_g1k_v37.fasta

        bwa mem -t 8 \
            -R "@RG\tID:FLOWCELL.LANE\tPL:ILLUMINA\tLB:test\tSM:PA01" \
            ${{BWAREF}} \
            {input} \
        | samtools view -b \
        | samtools sort -O BAM -T ${{TMPDIR}} -o {output.bam}

        samtools index {output.bam}
        """

rule structural_variants:
    input:
        'alignment/{id}.bam'
    output:
        'structural_variants/{id}.vcf'
    threads: 1
    resources:
        mem='4G',
        time='2-00:00:00',
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        delly call -o {output} -g ${{REF}} {input}
        """

def snps_mem(wildcards, attempt):
    mem = 2 * attempt
    return '%dG' % mem

rule snps:
    input:
        'alignment/{id}.bam'
    output:
        'snps/{id}.vcf'
    threads: 1
    resources:
        mem=snps_mem,
        time='04:00:00',
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        gatk HaplotypeCaller \
            -R ${{REF}} \
            -I {input} \
            -ploidy 2 \
            -O {output}
        """
```

We thus configure the resource consumption of the rules as follows:

- `alignment` with 8 threads and up to 8GB of memory in total with a running time of up to 12 hours,
- `structural_variants` with one thread and up to 4GB of memory in with a running time of up to 2 days,
- `snps` with one thread and running up to four hours.
  Instead of passing a static amount of memory, we pass a [resource callable](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=resources#resources).
  The `attempt` parameter will be passed a value of `1` on the initial invocation.
  If variant calling with the GATK HaplotypeCaller fails then it will retry and `attempt` will have an incremented value on each invocation (`2` on the first retry and so on).
  Thus, we try to do small variant calling with 2, 4, 6, and 8 GB.

Finally, run the script:

```terminal
(first-steps) $ sbatch submit_snakejob.sh
```

If you watch `squeue --me` now, you will see that the jobs are distributed to the system:

```terminal
(first-steps) $ squeue --me
```

Please refer to the Snakemake documentation for more details on using Snakemake, in particular [how to use the cluster configuration](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration) on how to specify the resource requirements on a per-rule base.
