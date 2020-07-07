# Slurm Command: `squeue`

The `squeue` command allows you to view currently running and pending jobs.

!!! info "Representative Example"

    ```bash
    med-login:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           1583165   highmem 20200702 usr      PD       0:00      1 (DependencyNeverSatisfied)
           1605901  critical variant_ holtgrem PD       0:00      1 (DependencyNeverSatisfied)
           1605902  critical variant_ holtgrem PD       0:00      1 (Dependency)
           1605905  critical variant_ holtgrem PD       0:00      1 (DependencyNeverSatisfied)
           1605916  critical wgs_sv_c holtgrem PD       0:00      1 (Dependency)
           1607103    medium wgs_sv_a holtgrem PD       0:00      1 (DependencyNeverSatisfied)
    [...]
    ```

!!! info "Slurm Documentation: squeue"

    Please also see the official [Slurm documentation on squeue](https://slurm.schedmd.com/squeue.html).

## Important Arguments

- `--nodelist`
    -- Only display jobs running on certain nodes (e.g., GPU nodes).
- `--format`
    -- Define the format to print, see `man squeue` for details.
    See below for a format string that includes the jobid, partition, job name, user name, job status, running time, number of nodes, number of CPU cores, and allocated GPUs.

## Notes

The following aliases in `~/.bashrc` will allow you to print a long and informative `squeue` output with `sq`, pipe it into less with `sql`, get only your jobs (adjust the `alias` to your account) using `sqme` and pipe that into less with `sqmel`.

```bash
alias sq='squeue -o "%.10i %9P %60j %10u %.2t %.10M %.6D %.4C %10R %b" "$@"'
alias sql='sq "$@" | less -S'
alias sqme='sq -u YOURUSER_c_or_m "$@"'
alias sqmel='sqme "$@" | less -S'
```
