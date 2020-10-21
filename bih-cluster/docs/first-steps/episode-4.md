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
(first-steps) $ mkdir -p /fast/users/${USER}/work/tutorial/episode4/logs
(first-steps) $ pushd /fast/users/${USER}/work/tutorial/episode4
```

And copy the wrapper script to this folder as well as the Snakefile (you can also reuse the one with the adjustments from the previous [episode](episode-3.md)):

```terminal
(first-steps) $ cp /fast/projects/cubit/work/tutorial/skeletons/submit_snakejob.sh .
(first-steps) $ cp /fast/projects/cubit/work/tutorial/skeletons/Snakefile .
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


export TMPDIR=/fast/users/${USER}/scratch/tmp
export LOGDIR=logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
mkdir -p $LOGDIR

unset DRMAA_LIBRARY_PATH
eval "$($(which conda) shell.bash hook)"
conda activate first-steps

set -x

# Note that Slurm DRMAA differs slightly from original Slurm syntax
# --mem-per-cpu doesn't accept units and the default unit here is MB
# -t only accepts HH:MM
snakemake \
    --drmaa " \
        -p medium \
        -t 01:00 \
        --nodes=1 \
        --mem=8192 \
        -n 8 \
        -o $LOGDIR/%x-%j.log" \
    -j 2 \
    -k \
    -p
```

In the beginning you see the `#SBATCH` that introduces the parameters when you provide this script to `sbatch`
as described in the [second episode](episode-2.md). Please make sure that the `logs` folder exists before starting the run!
We then set and export the `TMPDIR` and `LOGDIR` variables.
Note that `LOGDIR` has a subfolder named `$SLURM_JOB_NAME-$SLURM_JOB_ID` that will be created
for you. Snakemake will store its logfiles for this very Snakemake run in this folder.
The next new thing is `set -x`. This simply prints to the terminal every command that is executed within
the script. This is useful for debugging.

Finally, the Snakemake call takes place. With the `--drmaa` option we define how the jobs inside Snakemake should be started on the nodes. Please note that the `sbatch` command does not appear here, but the argument string provided to the `--drmaa` option is the same as the parameters for `sbatch`, except for some minor differences as descrived above or [here](../slurm/snakemake.md#limitations). The `drmaa` library provides more stable job distribution than using plain `sbatch` in this scenario. Note that we don't write the parameter string to a bash file (which we don't have in this case). There are three new parameters that are defining the hardware requirements for our run and the project:

* `--mem=X`: How much memory ONE job (=Snakemake rule) should get in Megabyte
* `-t HH:MM`: How long a job (=Snakemake rule) is allowed to run
* `-n X`: How many cores of a node a job (=Snakemake rule) should get
* `--nodes=1`: Force allocation of all cores on a single node.
* `-p medium`: The partition name (this system was introduced in [Episode 2](episode-2.md#job-queues) and is described [here](../overview/job-scheduler.md))

Note that the memory is not shared among the cores. This means the final memory on a node is defined by

```
<mem-per-cpu> * <n>
```

So in our example the total memory requested on one node would be 8GB.

The other parameters provided to Snakemake are:

* `-j`: Use 2 cores for this Snakemake run. In this scenario the parameter determines how many jobs will be submitted to SLURM at a time.
* `-k`: Keep going if a job fails
* `-p`: Print out shell commands

Finally, run the script:

```terminal
(first-steps) $ sbatch submit_snakejob.sh
```

If you watch `squeue --me` now, you will see that the jobs are distributed to the system:

```terminal
(first-steps) $ watch -n 60 squeue --me
```

Please refer to the Snakemake documentation for more details on using Snakemake, in particular [how to use the cluster configuration](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration) on how to specify the resource requirements on a per-rule base.
