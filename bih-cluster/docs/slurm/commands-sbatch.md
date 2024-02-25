# Slurm Command: `sbatch`

The `sbatch` command allows you to put a job into the scheduler's queue to be executed at a later time.

!!! info "Representative Example"

    ```bash
    # Execute job.sh in partition medium with 4 threads and 4GB of RAM total for a
    # running time of up to one day.
    med-login:~$ sbatch --partition=medium --mem=4G --ntasks 4 --time=1-00 job.sh
    Submitted batch job JOB_ID
    ```

The command will create a batch job and add it to the queue to be executed at a later point in time.

!!! info "Slurm Documentation: sbatch"

    Please also see the official [Slurm documentation on sbatch](https://slurm.schedmd.com/sbatch.html).

## Important Arguments

- `--array`
    -- Submit jobs as array jobs.
    Also see the section [#array-jobs] below.
- `--nodes`
    -- The number of nodes to allocate.
    This is only given here as an important argument as the maximum number of nodes allocatable to any partition but `mpi` is set to one (1).
    This is done as there are few users on the BIH HPC that actually use multi-node paralleilsm.
    Rather, most users will use multi-core parallelism and might forget to limit the number of nodes which causes inefficient allocation of resources.
- `--ntasks`
    -- This corresponds to the number of threads allocated to each node.
- `--mem`
    -- The memory to allocate for the job.
    As you can define minimal and maximal number of tasks/CPUs/cores, you could also specify `--mem-per-cpu` and get more flexible scheduling of your job.
- `--gres`
    -- Generic resource allocation.
    On the BIH HPC, this is only used for allocating GPUS, e.g., with `--gres=gpu:tesla:2`, a user could allocate two NVIDIA Tesla GPUs on the same host (use `a40` instead of `tesla` for the A40 GPUs).
- `--licenses`
    -- On the BIH HPC, this is used for the allocation of MATLAB 2016b licenses only.
- `--partition`
    -- The partition to run in.
    Also see the [Job Scheduler](../overview/job-scheduler.md) section.
- `--time`
    -- Specify the running time, see `man sbatch` or the official [Slurm documentation on srun](https://slurm.schedmd.com/srun.html) for supported formats.
    **Please note that the DRMA API only accepts the `hours:minutes` format.
- `--dependency`
    -- Specify dependencies on other jobs, e.g., using `--dependency afterok:JOBID` to only execute if the job with ID `JOBID` finished successfully or `--dependency after:JOBID` to wait for a job to finish regardless of its termination status.
- `--constraint`
    -- Require one or more features from your node.
    On the BIH HPC, the processor generation is defined as a feature on the nodes, e.g., `haswell`, or special networking such as `infiniband`.
    You can have a look at `/etc/slurm/slurm.conf` on all configured features.
- `--output`
    -- The path to the output log file (by default joining stdout and stderr, see the man page on `--error` on how to redirect stderr separately).
    A various number of placeholders is available, see the "filename pattern" section of `man sbatch` or the official [Slurm documentation on srun](https://slurm.schedmd.com/srun.html).
- `--mail-type=<type>`
  -- Send out notifications by email when an event occurs.
  Use `FAIL` to get emails when your job fails.
  Also see the documentation of [sbatch in the Slurm manual](https://slurm.schedmd.com/sbatch.html).
- `--mail-user=<email>`
  -- The email address to send to.
  Must end in `@charite.de`, `@mdc-berlin.de`, or `@bih-charite.de`.

!!! important "Ensure your `--output` directory exists!"

    In the case that the path to the log/output file does not exist, the job will just fail.
    `scontrol show job ID` will report `JobState=FAILED Reason=NonZeroExitCode`.
    Regrettably, no further information is displayed to you as the user.
    Always check that the path to the directories in `StdErr` and `StdOut` exists when checking `scontrol show job ID`.

## Other Arguments

- `--job-name`

## Job Scripts

Also see the section [Slurm Job Scripts](job-scripts.md) on how to embed the `sbatch` parameters in `#SBATCH` lines.

## Array Jobs

If you have many (say, more than 10) similar jobs (e.g., when performing a grid search), you can also use array jobs.
However, you should also consider whether it would make sense to increase the time of your jobs, e.g, to be at least ~10min.

You can submit array jobs by specifying `-a EXPR` or `--array EXPR` where `EXPR` is a range or a list (of course, you can also add this as an `#SBATCH` header in your job script).
For example:

```bash
hpc-login-1 ~# sbatch -a 1-3 grid_search.sh
hpc-login-1 ~# sbatch -a 1,2,5-10 grid_search.sh
```

This will submit `grid_search.sh` with certain variables set:

- `SLURM_ARRAY_JOB_ID` -- the ID of the first job
- `SLURM_ARRAY_TASK_ID` -- the index of the job in the array
- `SLURM_ARRAY_TASK_COUNT` -- number of submitted jobs in array
- `SLURM_ARRAY_TASK_MAX` -- higehst job array index value
- `SLURM_ARRAY_TASK_MIN` -- lowest job array index value

Using array jobs has several advantages:

- It greatly reduces the load on the Slurm scheduler.
- You do not need to submit in a loop, but rather
- You can use a single command line.

Also see [Slurm documentation on job arrays](https://slurm.schedmd.com/job_array.html).

For example, if you submit `sbatch --array=1-3 grid_search.sh` and slurm responsds with `Submitted batch job 36` then the script will be run three times with the following prameters set:

```

SLURM_JOB_ID=36
SLURM_ARRAY_JOB_ID=36
SLURM_ARRAY_TASK_ID=1
SLURM_ARRAY_TASK_COUNT=3
SLURM_ARRAY_TASK_MAX=3
SLURM_ARRAY_TASK_MIN=1

SLURM_JOB_ID=37
SLURM_ARRAY_JOB_ID=36
SLURM_ARRAY_TASK_ID=2
SLURM_ARRAY_TASK_COUNT=3
SLURM_ARRAY_TASK_MAX=3
SLURM_ARRAY_TASK_MIN=1

SLURM_JOB_ID=38
SLURM_ARRAY_JOB_ID=36
SLURM_ARRAY_TASK_ID=3
SLURM_ARRAY_TASK_COUNT=3
SLURM_ARRAY_TASK_MAX=3
SLURM_ARRAY_TASK_MIN=1
```

## Notes

- This is the primary entry point for creating batch jobs to be executed at a later point in time.
- As with all jobs allocated by Slurm, interactive sessions executed with `sbatch` are governed by resource allocations, in particular:
    - `sbatch` jobs have a maximal running time set,
    - `sbatch` jobs have a maximal memory and number of cores set, and
    -  also see [`scontrol show job JOBID`](commands-scontrol.md).
