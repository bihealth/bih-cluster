# Slurm and Temporary Files

This section describes how Slurm handles temporary files on the local disk.

!!! info "Temporary Files Best Practices"

    See [Best Practices: Temporary Files](../best-practices/temp-files/) for information how to use temporary files effectively.

## Slurm Behaviour

Our Slurm configuration has the following behaviour.

### Environment Variable TMPDIR

Slurm itself will by default **not** change the `TMPDIR` environment variable but **retain** the variable's value from the `srun` or `sbatch` call.

### Private Local `/tmp` Directories

The only place where users can write data to on local storage of the compute nodes is `/tmp`.

Storage is a consumable shared resource as the storage used by one job cannot use another job.
It is thus critical that Slurm cleans up after each job such that all space on the local node is available to the next job.
This is done using the [job_container/tmpfs](https://slurm.schedmd.com/job_container.conf.html) Slurm plugin.

This plugin creates a so-called Linux namespace for each job and creates a [bind mount](https://unix.stackexchange.com/questions/198590/what-is-a-bind-mount) of `/tmp` to a location on the local storage.
This mount is only visible to the currently running job and each job, even of the same user, get their own `/tmp`.
After a job terminates, Slurm will remove the directory and all of its content.

There is a notable exception.
If you use `ssh` to connect to a node rather than using `srun` or `sbatch`, you will see the system `/tmp` directory and can also write to it.
This usage of storage is not tracked and consequently you can circumvent the Slurm quota management.
Using `/tmp` in this fashion (i.e., outside of Slurm-controlled jobs) is prohibited.
If it cannot be helped (e.g., if you need to run some debugging application that needs to create FIFO or socket files) then keep usage of `/tmp` outside of Slurm job below 100MB.

### Tracking Local Storage `localtmp`

!!! warning "Enforcing `localtmp` Gres"

    From January 31, we will enforce the allocated storage in `/tmp` on the local disk with quotas.
    Jobs writing to `/tmp` beyond the quota in the job allocation will not function properly and probably crash with "out of disk quota" messages.

Slurm tracks the available local storage above 100MB on nodes in the `localtmp` generic resource (aka Gres).
The resource is counted in steps of 1MB, such that a node with 350GB of local storage would look as follows in `scontrol show node`:

```
hpc-login-1 # scontrol show node hpc-cpu-1
NodeName=hpc-cpu-1 Arch=x86_64 CoresPerSocket=24
   [...]
   Gres=localtmp:350K
   [...]
   CfgTRES=cpu=96,mem=360000M,billing=96,gres/localtmp=358400
   [...]
```

Each job is automaticaly granted 100MB of storage on the local disk which is sufficient for most standard programs.
If your job needs more temporary storage then you should either

- use the `$HOME/scratch` volume (see [Best Practices: Temporary Files](../best-practices/temp-files/))
- specify a `localtmp` generic resource (described here)

You can allocate the resource with `--gres=localtmp:SIZE` where `SIZE` is given in MB.

```
hpc-login-1 # srun --gres=localtmp:100k --pty bash -i
hpc-cpu-1 # scontrol show node hpc-cpu-1
NodeName=hpc-cpu-1 Arch=x86_64 CoresPerSocket=24
   [...]
   Gres=localtmp:250K
   [...]
   CfgTRES=cpu=96,mem=360000M,billing=96,gres/localtmp=358400
   [...]
   AllocTRES=cpu=92,mem=351G,gres/localtmp=102400
   [...]
```

The first output tells us about the resource configured to be available to user jobs and the last line show us that `100k=102400` MB of local storage are allocated.

You can also see the used resources in the details of your job:

```
scontrol show job 14848
JobId=14848 JobName=example.sh
   [...]
   TresPerNode=gres:localtmp:100k
```