# Slurm Command: `sinfo`

The `sinfo` command allows you to query the current cluster status.

!!! info "Representative Example"

    ```bash
    med-login1:~$ sinfo
    PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
    [...]
    medium       up 7-00:00:00     10 drain* med[0101-0103,0125-0126,0128-0132]
    medium       up 7-00:00:00      1  down* med0243
    medium       up 7-00:00:00     31    mix med[0104,0106-0122,0124,0133,0232-0233,0237-0238,0241-0242,0244,0263-0264,0503,0506]
    medium       up 7-00:00:00      5  alloc med[0105,0123,0127,0239-0240]
    medium       up 7-00:00:00    193   idle med[0134-0164,0201-0231,0234-0236,0245-0262,0501-0502,0504-0505,0507-0516,0601-0632,0701-0764]
    [...]
    med-login1:$  sinfo --summarize
    PARTITION AVAIL  TIMELIMIT   NODES(A/I/O/T) NODELIST
    debug*       up    8:00:00    38/191/11/240 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
    medium       up 7-00:00:00    38/191/11/240 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
    long         up 28-00:00:0    38/191/11/240 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
    critical     up 7-00:00:00    25/141/10/176 med[0101-0164,0501-0516,0601-0632,0701-0764]
    highmem      up 14-00:00:0          1/2/1/4 med[0401-0404]
    gpu          up 14-00:00:0          3/0/1/4 med[0301-0304]
    mpi          up 14-00:00:0    38/191/11/240 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
    ```

This command will summaries the state of nodes by different criteria (e.g., by partition or globally).

!!! info "Slurm Documentation: sinfo"

    Please also see the official [Slurm documentation on srun](https://slurm.schedmd.com/sinfo.html).

## Important Arguments

Also see all important arguments of the [`sinfo`](commands-sinfo.md) command.

- `--summarize`
    -- Summarize the node state by partition.
- `--nodes`
    -- Select the nodes to show the status for, e.g., display the status of all GPU nodes with `sinfo -n med030[1-4]`.

## Node States

The most important node states are:

- `down` -- node is marked as offline
- `draining` -- node will not accept any more jobs but has jobs running on it
- `drained` -- node will not accept any more jobs and has no jobs running on it, but is not offline yet
- `idle` -- node is ready to run jobs
- `allocated` -- node is fully allocated (e.g., CPU, RAM, or GPU limit has been reached)
- `mixed` -- node is running jobs but there is space for more

## Notes

- Also see the [Slurm Format Strings](format-strings.md) section.
