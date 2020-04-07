# First Steps: Episode 2

|Episode|Topic|
|---|---|
| 0 | [How can I install the tools?](episode-0.md) |
| 1 | [How can I use the static data?](episode-1.md) |
| **2** | **How can I distribute my jobs on the cluster (`qsub`)?** |
| 3 | [How can I organize my jobs with Snakemake?](episode-3.md) |
| 4 | [How can I combine Snakemake and `qsub`?](episode-4.md) |

Welcome to the second episode of our tutorial series!

Once you are logged in to the cluster, you have the possibility to distribute your jobs to all the nodes that are available.
But how can you do this easily?
The key command to this magic is `qsub`.
This tutorial will show you how you can use this efficiently.

!!! hint
    Before using `qsub`, make sure that you logged in using `qrsh` once you entered the login node via ssh, as described in [Connecting via SSH](../connecting/configure-ssh/linux.md).

## The `qsub` Command

So what is `qsub` doing for you?

You use the `qsub` command in front of the script you actually want to run.
`qsub` then puts your job into the job queue.
The job scheduler looks at the current status of the whole system and will assign the first job in the queue to a node that is free in terms of computational load.
If all machines are busy, yours will wait.
But your job will sooner or later get assigned to a free node.

We strongly recommend using this process for starting your computationally intensive tasks because you will get the best performance for your job and the
whole system won't be disturbed by jobs that are locally blocking nodes.
Thus, everybody using the cluster benefits.

You may have noticed that you run `qsub` with a script, not with regular commands.
The reason is that `qsub` only accepts bash scripts.
If you give `qsub` a normal shell command or binary, it won't work.
This means that we have to put the command(s) we want to use in a bash script.
A skeleton script can be found at `/fast/projects/cubit/work/tutorial/skeletons/submit_job.sh`

The content of the file:

```bash
#!/bin/bash

# Define the expected running time of your job.
#$ -l h_rt=24:00:00

# Define the project you want to enqueue your job to,
# according to the expected running time.
#$ -P medium

# Combine stderr and stdout log files into the stdout log file.
#$ -j y

# Keep current environment variables.
#$ -V

# Use current working directory as working directory of the job.
#$ -cwd

# Set the log directory.
#$ -o logs

export TMPDIR=/fast/users/${USER}/scratch/tmp
mkdir -p ${TMPDIR}
```

The lines starting with `#$` are actually setting parameters for a `qsub` command, so `#$ -j y` is equal to `qsub -j y`.
If the `logs` value in the `-o` parameter refers to an existing directory, `qsub` will create a logfile inside this directory.
If it does not exist, `qsub` will create a file named `logs` and put the log content there.
We recommend having the directory available before starting the job.
We also emphasize the importance of the log files!
They are the first place to look if anything goes wrong.

To start now with our tutorial, create a new tutorial directory with a log directory, e.g.,

```terminal
(first-steps) $ mkdir -p /fast/users/$USER/work/tutorial/episode2/logs
```

and copy the wrapper script to this directory:

```terminal
(first-steps) $ pushd /fast/users/$USER/work/tutorial/episode2
(first-steps) $ cp /fast/projects/cubit/work/tutorial/skeletons/submit_job.sh .
(first-steps) $ chmod u+w submit_job.sh
```

Now open this file and copy the same commands we executed in the last tutorial to this file.

To keep it simple, we will put everything into one script.
This is perfectly fine because the alignment and indexing are sequential.
But there are two steps that could be run in parallel, namely the variant calling, because they don't depend on each other.
We will learn how to do that in a later tutorial.
Your file should look something like this:

```bash
#!/bin/bash

# Define the expected running time of your job.
#$ -l h_rt=24:00:00

# Define the project you want to enqueue your job to,
# according to the expected running time.
#$ -P medium

# Combine stderr and stdout log files into stdout log file.
#$ -j y

# Keep current environment variables.
#$ -V

# Use current working directory as working directory of the job.
#$ -cwd

# Set the log directory.
#$ -o logs

export TMPDIR=/fast/users/${USER}/scratch/tmp
mkdir -p ${TMPDIR}

BWAREF=/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.15/GRCh37/g1k_phase1/human_g1k_v37.fasta
REF=/fast/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

bwa mem -t 8 \
    -R "@RG\tID:FLOWCELL.LANE\tPL:ILLUMINA\tLB:test\tSM:PA01" \
    $BWAREF \
    /fast/projects/cubit/work/tutorial/input/test_R1.fq.gz \
    /fast/projects/cubit/work/tutorial/input/test_R2.fq.gz \
| samtools view -b - \
| samtools sort -O BAM -T $TMPDIR -o aln.bam

samtools index aln.bam

delly call -g \
    $REF \
    aln.bam

gatk-launch HaplotypeCaller \
    -R $REF \
    -I aln.bam \
    -ploidy 2 \
    -O test.GATK.vcf
```

Let's run it (make sure that you are in the `tutorial/episode2` directory!):

```terminal
(first-steps) $ qsub submit_job.sh
```

And wait for the response which will tell you that your job was submitted and which job id number it was assigned. Note that `qsub` only tells you that the job has started, but nothing about finishing. You won't get any response at the terminal when the job finishes. It will take approximately 20 minutes to finish the job.

## Monitoring Jobs

You'll probably want to see how your job is doing. You can get a list of your jobs using:

```terminal
(first-steps) $ qstat
```

Note that logins are also considered as jobs, albeit "interactive" and not "batch" jobs.

Identify your job by the `<JOBID>` (1st column) or the name of the script (3rd column).
The most likely states you will see (4th column of the table):

* `qw` waiting to be submitted
* `r` running
* `Eqw` in error queue
* disappeared, either because of an error or because it finished

In the 8th column (queue) you can see that your job is very likely running on a different machine than the one you are on!

Get more information about your jobs by either passing the job id when it is in the `r` or `Eqw` state:

```terminal
(first-steps) $ qstat -j <JOBID>
```

Once a job is complete, you can use `qacct -j <JOBID>` to see information on your job.

And of course, watch what the logs are telling you:

```terminal
(first-steps) $ tail -f logs/submit_job.sh.o<JOBID>
```

There will be no notification when your job is done, so it is best to watch the `qstat` command.
To watch the `qsub` command there is a linux command `watch` that you give a command to execute every few seconds.
This is useful for looking for changes in the output of a command. The seconds between two executions can be set with the `-n` option.
:warning: It is best to use `-n 60` to minimize unnecessary load on the file system:

```terminal
(first-steps) $ watch -n 60 qstat
```
If for some reason your job is hanging (eg. file_not_found situations), you can delete your job using `qdel` with your job-ID:
```terminal
(first-steps) $ qdel <job-ID>
```

## Job Queues

You might have noticed that the job you submitted when you monitor it with `qstat`
is running in the `medium.q` queue (in the 7th column named `queue` you will see something like `medium.q@med0...`).

The cluster has a special way of organizing itself and by telling the cluster how long and with which priority you want your jobs to run, you can help it in this.
There is a system set up on the cluster where you can enqueue your jobs to so-called *projects*. *projects* have different prioritites and are allowed for different running times.
Those *projects* have queues assigned carrying the same name as the project, e.g. submitting a job to project `long` might end up in queue `long.q`, but not necessarily.
To get to know why this is, what projects are available, and how to use them properly, we highly encourage you to read the [cluster queues wiki page](../overview/job-scheduler.md).

