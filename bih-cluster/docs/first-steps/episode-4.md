# First Steps: Episode 4

|Episode|Topic|
|---|---|
| 0 | [How can I install the tools?](episode-0.md) |
| 1 | [How can I use the static data?](episode-1.md) |
| 2 | [How can I distribute my jobs on the cluster (`qsub`)?](episode-2.md) |
| 3 | [How can I organize my jobs with Snakemake?](episode-3.md) |
| **4** | **How can I combine Snakemake and `qsub`?** |

In the last episodes we learned about distributing a job among the cluster nodes using `qsub` and
how to automate and parallelize our pipeline with Snakemake. We are lucky that those two
powerful commands can be combined. What is the result? You will have an automated pipeline
with Snakemake that uses `qsub` to distribute jobs among the cluster nodes instead of
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

#$ -V
#$ -j y
#$ -o logs
#$ -r yes
#$ -cwd
#$ -S /bin/bash
#$ -P control

export TMPDIR=/fast/users/${USER}/scratch/tmp
export LOGDIR=logs/${JOB_ID}
mkdir -p $LOGDIR

set -x

snakemake \
    --drmaa " \
        -V \
        -cwd \
        -P medium \
        -l h_vmem=5g \
        -l h_rt=72:00:00 \
        -pe smp 8 \
        -j yes \
        -o $LOGDIR/" \
    -j 2 \
    -k \
    -p
```

In the beginning you see the `#$` that introduces the parameters when you provide this script to `qsub`
as described in the [second episode](episode-2.md). Please make sure that the log folder exists before starting the run!
We then set and export the `TMPDIR` and `LOGDIR` variables.
Note that `LOGDIR` has a subfolder named `$JOB_ID` that will be created
for you. Snakemake will store its logfiles for this very Snakemake run in this folder.
The next new thing is `set -x`. This simply prints to the terminal every command that is executed within
the script. This is useful for debugging.

Finally, the Snakemake call takes place. With the `--drmaa` option we define how the jobs inside
Snakemake should be started on the nodes. Please note that the `qsub` command does not appear here, but the argument string provided to the `--drmaa` option is the same as the parameters for `qsub`. The `drmaa` library provides more stable job distribution than using plain `qsub` in this scenario. Note that we don't write the parameter string to a bash file (which we don't have in this case). There are three new parameters that are defining the hardware requirements for
our run and the project:

* `-l h_vmem`: How much memory ONE job (=Snakemake rule) should get [Xg]
* `-l h_rt`: How long a job (=Snakemake rule) is allowed to run [HH:MM:SS]
* `-pe smp 8`: How many cores of a node a job (=Snakemake rule) should get [N]
* `P medium`: The project name (this system was introduced in [Episode 2](episode-2.md#job-queues) and is described [here](../overview/job-scheduler.md))

Note that the memory is not shared among the cores. This means the final memory on a node is defined by

```
h_vmem=Xg * smp N
```

So in our example the total memory requested on one node would be 40GB.

The other parameters provided to Snakemake are:

* `-j`: Use 2 cores for this Snakemake run
* `-k`: Keep going if a job fails
* `-p`: Print out shell commands

Finally, run the script:

```terminal
(first-steps) $ qsub submit_snakejob.sh
```

If you watch `qstat` now, you will see that the jobs are distributed to the system:

```terminal
(first-steps) $ watch -n 60 qstat
```

Please refer to the Snakemake documentation for more details on using Snakemake, in particular [how to use the cluster configuration](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration) on how to specify the resource requirements on a per-rule base.
