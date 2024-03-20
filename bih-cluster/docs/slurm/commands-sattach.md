# Slurm Command: `sattach`

The `sattach` command allows you to connect the standard input, output, and error streams to your current terminals ession.

!!! info "Representative Example"

    ```bash
    hpc-login-1:~$ sattach 12345.0
    [...output of your job...]
    med0211:~$ [Ctrl-C]
    hpc-login-1:~$
    ```

Press `Ctrl-C` to detach from the current session.
Please note that you will have to give the job ID as well as step step ID.
For most cases, simply append `".0"` to your job ID.

!!! info "Slurm Documentation: sattach"

    Please also see the official [Slurm documentation on srun](https://slurm.schedmd.com/sbatch.html).

## Important Arguments

- `--pty`
    -- Execute task zero in pseudo terminal.
- `--verbose`
    -- Increase verbosity of `sattach`.
