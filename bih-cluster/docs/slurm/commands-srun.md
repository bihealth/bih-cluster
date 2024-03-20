# Slurm Command: `srun`

The `srun` command allows you to run a command **now**.

!!! info "Representative Example"

    ```bash
    hpc-login-1:~$ srun --pty bash -i
    med0201:~$
    ```

The command will perform a resource allocation with the scheduler (and wait until it has allocated the requested resources) first.
Most importantly, you can specify the `--pty` argument which will connect the current terminal's standard output, error, and input to your current one.
This allows you to run interactive jobs such as shells with `srun --pty bash -i`.

!!! info "Slurm Documentation: srun"

    Please also see the official [Slurm documentation on srun](https://slurm.schedmd.com/srun.html).

## Important Arguments

Also see all important arguments of the [`sbatch`](commands-sbatch.md) command.

- `--pty`
  -- Connect current terminal to the job's stdoud/stderr/stdin.
- `--x11`
  -- Setup X11 forwarding.
- `--immediate`
  -- Immediately terminate if the resources to run the job are not available, do not wait.
- `--test-only`
  -- Don't run anything, but only estimate when the job would be scheduled.

## Notes

- This is the primary entry point for creating interactive shell sessions on the cluster.
- As with all jobs allocated by Slurm, interactive sessions executed with `srun` are governed by resource allocations, in particular:
    - `srun` jobs have a maximal running time set,
    - `srun` jobs have a maximal memory and number of cores set, and
    -  also see [`scontrol show job JOBID`](commands-scontrol.md).
