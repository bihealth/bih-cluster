# Slurm Command: `sbatch`

The `srun` command allows you to put a job into the scheduler's queue to be executed at a later time.

!!! info "Representative Example"

    ```bash
    # Execute job.sh in partition medium with 4 threads and 4GB of RAM total for a
    # running time of up to one day.
    med-login:~$ sbatch --partition=medium --mem=4G --ntasks 4 --time=1-00 job.sh
    Submitted batch job JOB_ID
    ```

The command will create a batch job and add it to the queue to be executed at a later point in time.

!!! info "Slurm Documentation: srun"

    Please also see the official [Slurm documentation on srun](https://slurm.schedmd.com/srun.html).

## Important Arguments

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
    On the BIH HPC, this is only used for allocating GPUS, e.g., with `--gres=gpu:tesla:2`, a user could allocate two NVIDIA Tesla GPUs on the same hsot.
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

!!! important "Ensure your `--output` directory exists!"

    In the case that the path to the log/output file does not exist, the job will just fail.
    `scontrol show job ID` will report `JobState=FAILED Reason=NonZeroExitCode`.
    Regrettably, no further information is displayed to you as the user.
    Always check that the path to the directories in `StdErr` and `StdOut` exists when checking `scontrol show job ID`.

## Other Arguments

- `--job-name`

## Job Scripts

Also see the section [Slurm Job Scripts](job-scripts.md) on how to embed the `sbatch` parameters in `#SBATCH` lines.

## Notes

- This is the primary entry point for creating batch jobs to be executed at a later point in time.
- As with all jobs allocated by Slurm, interactive sessions executed with `sbatch` are governed by resource allocations, in particular:
    - `sbatch` jobs have a maximal running time set,
    - `sbatch` jobs have a maximal memory and number of cores set, and
    -  also see [`scontrol show job JOBID`](commands-scontrol.md).
