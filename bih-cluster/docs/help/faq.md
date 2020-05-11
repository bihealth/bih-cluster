# Frequently Asked Questions

## What is this Website?

This is the BIH cluster documentation that was created and is maintained by BIH Core Unit Bioinformatics (CUBI) and BIH HPC IT with contributions by BIH HPC Users.
The aim is to gather the information for using the cluster efficiently and helping common issues andproblems.

## Where can I get help?

- First, contact [bih-cluster@charite.de](mailto:bih-cluster@charite.de) with your question.
  Here, administrators, CUBI, and other users alike are subscribed and can answer your question.
- For problems while connecting and logging in, please contact [helpdesk@mdc-berlin.de](mailto:helpdesk@mdc-berlin.de) or [helpdesk@charite.de](mailto:helpdesk@charite.de).
- For problems with BIH HPC please contact [hpc-helpdesk@bihealth.de].

## I cannot connect to the cluster. What's wrong?

Please see the section [Connection Problems](/connecting/configure-ssh/connection-problems).

## I'd like to learn more about Slurm

- Some documentation is available on this website, e.g., start at [Slurm Quickstart](../slurm/quickstart).

## What is the difference between MAX and BIH cluster? What is their relation?

**Administrativa**

- The BIH cluster is the cluster of the Berlin Institute of Health (BIH) and is located in Buch and operated by MDC IT.
  MDC IT performs the administration for the BIH and the BIH cluster uses the authentication (user/password/groups) infrastructure of the MDC.
  Thus you have to get an MDC account for using the BIH cluster.
  Systems support can be requested from helpdesk@mdc-berlin.de as for the MAX cluster.
- The MAX cluster is the cluster of the Max Delbrueck Center (MDC) in Buch.
  This cluster is used by the researchers at MDC and integrates a lot of infrastructure of the MDC.

Request for both systems are handled separately, depending on the user's affiliation with research/service groups.

**Hardware and Systems**

- Both clusters consist of similar hardware for the compute nodes and both feature a DDN system at different number of nodes and different storage volume.
- Both clusters run CentOS but at potentially different version.
- BIH HPC uses the Slurm workload manager whereas MAX uses Univa Grid Engine.
- The BIH cluster has a significantly faster internal network (40GB/s optical).

**Bioinformatics Software**

- On the BIH cluster, users can install their own bioinformatics software in their user directory.
- On the MAX cluster, users can also install their own software or use [software provided by Altuna Akalin's group at MDC](http://bioinformatics.mdc-berlin.de/resources.html#other-it-services).

## My SSH sessions break with "`packet_write_wait: Connection to XXX : Broken pipe`". How can I fix this?

Try to put the following line at the top of your `~/.ssh/config`.

```
ServerAliveInterval 30
```

This will make `ssh` send an empty network package to the server.
This will prevent network hardware from thinking your connection is unused/broken and terminating it.

If the problem persists, please report it to hpc-helpdesk@bihealth.de.

## My job terminated before being done. What happened?

First of all, look into your job logs.
In the case that the job was terminated by Slurm (e.g., because it ran too long), you will find a message at the bottom.

Otherwise, you can use `sacct -j JOBID` to read the information that the job accounting system has recorded for your job.
Use the `--long` flag to see all fields (and probably pipe it into `less` as: `sacct -j JOBID --long | less -S`).
Things to look out for:

- What is the exit code?
- Is the highest recorded memory usage too high/higher than expected (field `MaxRSS`)?
- Is the running time too long/longer than expected (field `Elapsed`)?

## How can I create a new project?

You can create a project if you are either a group leader of an AG or a delegate of an AG.

In this case, please send an email to hpc-admin@bihealth.de and request a new project.

## I have a problem with my SSH key, help!

Please contact [MDC Helpdesk](mailto:helpdesk@mdc-berlin.de) for help.
CUBI is also only a user on the BIH cluster and has no access to the user, password, or SSH keys registries.

## I have a problem with logging into the cluster, help!

See ["I have a problem with my SSH key, help!"](#i-have-a-problem-with-my-ssh-key-help)

## I cannot create PNGs in R

For using the `png` method, you need to have an X11 session running.
This might be the case if you logged into a cluster node using `srun --x11` if configured correctly but is not the case if you submitted a bash job.
The solution is to use `xvfb-run` (xvfb = X11 virtual frame-buffer).

Here is the content of an example script:

```terminal
$ cat img.R 
#!/usr/bin/env Rscript

png('cars.png')
cars <- c(1, 3, 6, 4, 9)
plot(cars)
dev.off()
```

Here, it fails without X11:

```terminal
$ ./img.R 
Error in .External2(C_X11, paste("png::", filename, sep = ""), g$width,  : 
  unable to start device PNG
Calls: png
In addition: Warning message:
In png("cars.png") : unable to open connection to X11 display ''
Execution halted
```

Here, it works with  `xvfb-run`:

```terminal
$ xvfb-run ./img.R 
null device 
          1 
$ ls
cars.png  foo.png  img.R  Rplots.pdf
```

## My jobs don't get scheduled

You can use `scontrol show job JOBID` to get the details displayed about your jobs.
In the example below, we can see that the job is in the `PENDING` state.
The `Reason` field tells us that the job did not scheduled because the specified dependency was neverfulfilled.
You can find a list of all job reason codes in the [Slurm `squeue` documentation](https://slurm.schedmd.com/squeue.html#lbAF).

```bash hl_lines="4"
JobId=863089 JobName=pipeline_job.sh
   UserId=holtgrem_c(100131) GroupId=hpc-ag-cubi(5272) MCS_label=N/A
   Priority=1 Nice=0 Account=(null) QOS=normal
   JobState=PENDING Reason=DependencyNeverSatisfied Dependency=afterok:863087(failed)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=08:00:00 TimeMin=N/A
   SubmitTime=2020-05-03T18:57:34 EligibleTime=Unknown
   AccrueTime=Unknown
   StartTime=Unknown EndTime=Unknown Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2020-05-03T18:57:34
   Partition=debug AllocNode:Sid=med-login1:28797
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=(null)
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,node=1,billing=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/wgs_cnv_export/pipeline_job.sh
   WorkDir=/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/wgs_cnv_export
   StdErr=/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/wgs_cnv_export/slurm-863089.out
   StdIn=/dev/null
   StdOut=/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/wgs_cnv_export/slurm-863089.out
   Power=
   MailUser=(null) MailType=NONE

```

## My jobs don't run in the partition I expect

You can see the partition that your job runs in with `squeue -j JOBID`:

```bash
med-login1:~$ squeue -j 877092
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            877092    medium snakejob holtgrem  R       0:05      1 med0626
```

See [Job Scheduler](/overview/job-scheduler/) for information about the partition's properties.
To get your job to run in the `medium` partition, for example use the `--partition=medium` or `-p medium` arguments to your `srun` or `sbatch` commands.

## My jobs get killed after four hours

This is probably answered by the answer to [My jobs don't run in the partition I expect](#my-jobs-dont-run-in-the-partition-i-expect).

## How can I mount a network volume from elsewhere on the cluster?

You cannot.
Also see the [For the Impatient](../overview/for-the-impatient) section of the manual.

## Why can I not mount a network volume from elsewhere on the cluster?

For performance and stability reasons.
Network volumes are notorious for degrading performance, depending on the used protocol, even stability.

## How can I then access the files from my workstation/server?

You can transfer files to the cluster through Rsync over SSH or through SFTP to the `med-transfer1` or `med-transfer2` node.

**Do not transfer files through the login nodes.**
Large file transfers through the login nodes can cause performance degradation for the users with interactive SSH connections.

## How can I circumvent "invalid instruction" (signal 4) errors?

Make sure that software is compiled with "sandy bridge" optimizations and no later one.
E.g., use the `-march=sandybridge` argument to the GCC/LLVM compiler executables.

If you absolutely need it, there are some boxes with more recent processors in the cluster (e.g., Haswell architecture).
Look at the `/proc/cpuinfo` files for details.

## Where should my (Mini)conda install go?

As conda installations are big and contain many files, they should go into your `work` directory.
**E.g., `/fast/users/$USER/work/miniconda` is appropriate.**

## I have problems connecting to the GPU node! What's wrong?

Please check whether there might be other jobs waiting in front of you!
The following `squeue` call will show the allocated GPUs of jobs in the `gpu` queue.
This is done by specifying a format string and using the `%b` field.

```bash
squeue -o "%.10i %9P %20j %10u %.2t %.10M %.6D %10R %b" -p gpu
     JOBID PARTITION NAME                 USER       ST       TIME  NODES NODELIST(R TRES_PER_NODE
    872571 gpu       bash                 user1       R   15:53:25      1 med0303    gpu:tesla:1
    862261 gpu       bash                 user2       R 2-16:26:59      1 med0304    gpu:tesla:4
    860771 gpu       kidney.job           user3       R 2-16:27:12      1 med0302    gpu:tesla:1
    860772 gpu       kidney.job           user3       R 2-16:27:12      1 med0302    gpu:tesla:1
    860773 gpu       kidney.job           user3       R 2-16:27:12      1 med0302    gpu:tesla:1
    860770 gpu       kidney.job           user3       R 4-03:23:08      1 med0301    gpu:tesla:1
    860766 gpu       kidney.job           user3       R 4-03:23:11      1 med0303    gpu:tesla:1
    860767 gpu       kidney.job           user3       R 4-03:23:11      1 med0301    gpu:tesla:1
    860768 gpu       kidney.job           user3       R 4-03:23:11      1 med0301    gpu:tesla:1
```

In the example above, user1 has one job with one GPU running on med0303, user2 has one job running with 4 GPUs on med0304 and user3 has 7 jobs in total running of different machines with one GPU each.

## How can I access graphical user interfaces (such as for Matlab) on the cluster?

1. First of all, you will need an X(11) server on your local machine (see [Wikipedia: X Window System](https://en.wikipedia.org/wiki/X_Window_System).
  This server offers a "graphical surface" that the programs on the cluster can then paint on.
2. You need to make sure that the programs running on the cluster can access this graphical surface.
    - Generally, you need to connect to the login nodes with X forwarding.
      Refer to the manual of your SSH client on how to do this (`-X` for Linux/Mac `ssh`
    - As you should not run compute-intensive programs on the login node, connect to a cluster node with X forwarding.
      With Slurm, this is done using `srun --pty --x11 bash -i` (instead of `srun --pty --x11 bash -i`).

Also see:

- [Running graphical(X11) applications on Windows](../connecting/configure-ssh/windows.md#x11)
- [Running graphical(X11) applications on Linux](../connecting/configure-ssh/linux.md#x11)

## How can I log into a node outside of the scheduler?

This is sometimes useful, e.g., for monitoring the CPU/GPU usage of your job interactively.

!!! warning "No Computation Outside of Slurm"

    Do **not** perform any computation outside of the scheduler as (1) this breaks the purpose of the scheduling system and (2) administration is not aware and might kill you jobs.

The answer is simple, just SSH into this node.

```bash
med-login1:~$ ssh med0XXX
```

## Snakemake DRMAA doesn't accept my Slurm parameters!?

Yes. Unfortunately, [Slurm DRMAA differs slightly](../slurm/snakemake.md#limitations) from the original Slurm syntax.

## Why am I getting multiple nodes to my job?

Classically, jobs on HPC systems are written in a way that they can run on multiple nodes at once, using the network to communicate.
Slurm comes from this world and when allocating more than one CPU/core, it might allocate them on different nodes.
Please use `--nodes=1` to force Slurm to allocate them on a single node.

## How can I select a certain CPU architecture?

You can select the CPU architecture by using the `-C`/`--constraint` flag to `sbatch` and `srun`.
The following are available (as detected by the Linux kernel):

- `ivybridge` (96 nodes, plus 4 high-memory nodes)
- `haswell` (16 nodes)
- `broadwell` (112 nodes)
- `skylake` (16 nodes, plus 4 GPU nodes)

You can specify contraints with OR such as `--constraints=haswell|broadwell|skylake`.
You can see the assignment of architectures to nodes using the `sinfo -o "%8P %.5a %.10l %.6D %.6t %10f %N"` command.
This will also display node partition, availability etc.
