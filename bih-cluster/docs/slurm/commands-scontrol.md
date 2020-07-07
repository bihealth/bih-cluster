# Slurm Command: `scontrol`

The `scontrol` allows to query detailed information from the scheduler and perform manipulation.
Object manipulation is less important for normal users.

!!! info "Representative Example"

    ```bash
    med-login:~$ scontrol show job 1607103
    JobId=1607103 JobName=wgs_sv_annotation
        UserId=holtgrem_c(100131) GroupId=hpc-ag-cubi(5272) MCS_label=N/A
        Priority=748 Nice=0 Account=(null) QOS=normal
        [...]
    med-login:~$ scontrol show node med02[01-32]
    NodeName=med0201 Arch=x86_64 CoresPerSocket=8
        CPUAlloc=0 CPUTot=32 CPULoad=0.01
        AvailableFeatures=ivybridge,infiniband
        ActiveFeatures=ivybridge,infiniband
        [...]
    med-login:~$ scontrol show partition medium
    PartitionName=medium
        AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
        AllocNodes=ALL Default=NO QoS=medium
        DefaultTime=NONE DisableRootJobs=NO ExclusiveUser=NO GraceTime=0 Hidden=NO
        [...]
    ```

This command allows to query all information for an object from Slurm, e.g., jobs, nodes, or partitions.
The command also accepts ranges of jobs and hosts.
It is most useful to get the information of one or a few objects from the scheduler.

!!! info "Slurm Documentation: scontrol"

    Please also see the official [Slurm documentation on scontrol](https://slurm.schedmd.com/scontrol.html).

## Important Sub commands

- `scontrol show job`
    -- Show details on jobs.
- `scontrol show partition`
    -- Show details on partitions.
- `scontrol show node`
    -- Show details on nodes.
- `scontrol help`
    -- Show help.
- `scontrol`
    -- Start an interactive scontrol shell / REPL (read-eval-print loop).

## Notes

- `scontrol` can only work on jobs that are pending (in the queue), running, or in "completing' state.
- For jobs that have finished, you have to use Slurm's accounting features, e.g., with the [`sacct`](commands-sacct.md) command.
