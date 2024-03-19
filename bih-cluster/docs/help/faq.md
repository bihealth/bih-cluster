# Frequently Asked Questions

## What is this Website?

This is the BIH cluster documentation that was created and is maintained by BIH Core Unit Bioinformatics (CUBI) and BIH HPC IT with contributions by BIH HPC Users.
The aim is to gather the information for using the cluster efficiently and helping common issues andproblems.

## Where can I get help?

- Talk to your colleagues!
- Have a look at our forums at [HPC-talk](https://hpc-talk.cubi.bihealth.org/) to see if someone already solved the same problem.
  If not, create a new topic. Administrators, CUBI, and other users can see and answer your question.
- For problems while connecting and logging in, please contact [helpdesk@mdc-berlin.de](mailto:helpdesk@mdc-berlin.de) or [helpdesk@charite.de](mailto:helpdesk@charite.de).
- For problems with BIH HPC please contact [hpc-helpdesk@bih-charite.de].

## I cannot connect to the cluster. What's wrong?

Please see the section [Connection Problems](../connecting/configure-ssh/connection-problems.md).

## I'd like to learn more about Slurm

- Some documentation is available on this website, e.g., start at [Slurm Quickstart](../slurm/quickstart.md).

## What is the difference between MAX and BIH cluster? What is their relation?

**Administrativa**

- The BIH HPC 4 Research cluster of the Berlin Institute of Health (BIH) is located in Buch and operated by BIH HPC IT.
  The cluster is open for users of both BIH/Charite and MDC.
- The MAX cluster is the cluster of the Max Delbrueck Center (MDC) in Buch.
  This cluster is used by the researchers at MDC and integrates with a lot of infrastructure of the MDC.

Request for both systems are handled separately, depending on the user's affiliation with research/service groups.

**Hardware and Systems**

- Both clusters consist of similar hardware for the compute nodes and both feature a DDN system at different number of nodes and different storage volume.
- Both clusters run CentOS/rocky but at potentially different version.
- BIH HPC uses the Slurm workload manager whereas MAX uses Univa Grid Engine.
- The BIH cluster has a significantly faster internal network (40GB/s optical).

**Bioinformatics Software**

- On the BIH cluster, users can install their own (bioinformatics) software in their user directory.
- On the MAX cluster, users can also install their own software or use [software provided by Altuna Akalin's group at MDC](http://bioinformatics.mdc-berlin.de/resources.html#other-it-services).

## My SSH sessions break with "`packet_write_wait: Connection to XXX : Broken pipe`". How can I fix this?

Try to put the following line at the top of your `~/.ssh/config`.

```
ServerAliveInterval 30
```

This will make `ssh` send an empty network package to the server.
This will prevent network hardware from thinking your connection is unused/broken and terminating it.

If the problem persists, please report it to hpc-helpdesk@bih-charite.de.

## My job terminated before being done. What happened?

First of all, look into your job logs.
In the case that the job was terminated by Slurm (e.g., because it ran too long), you will find a message like this at the bottom.
Please look at the end of the last line in your log file.

```
slurmstepd: error: *** JOB <your job id> ON med0xxx CANCELLED AT 2020-09-02T21:01:12 DUE TO TIME LIMIT ***
```

This indicates that you need to need to adjust the `--time` limit to your `sbatch` command.

```
slurmstepd: error: Detected 2 oom-kill event(s) in step <your job id>.batch cgroup.
Some of your processes may have been killed by the cgroup out-of-memory handler
```

This indicates that your job tries to use more memory than has been allocated to it.
Also see [Slurm Scheduler: Memory Allocation](../slurm/memory-allocation.md)

Otherwise, you can use `sacct -j JOBID` to read the information that the job accounting system has recorded for your job.
A job that was canceled (indicated by `CANCELED`) by the Slurm job scheduler looks like this (ignore the `COMPLETED` step that is just some post-job step added by Slurm automatically).

```
# sacct -j _JOBID_
       JobID    JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
_JOBID_      snakejob.+     medium hpc-ag-xx+          4    TIMEOUT      0:0
_JOBID_.bat+      batch            hpc-ag-xx+          4  CANCELLED     0:15
_JOBID_.ext+     extern            hpc-ag-xx+          4  COMPLETED      0:0
```

Use the `--long` flag to see all fields (and probably pipe it into `less` as: `sacct -j JOBID --long | less -S`).
Things to look out for:

- What is the exit code?
- Is the highest recorded memory usage too high/higher than expected (field `MaxRSS`)?
- Is the running time too long/longer than expected (field `Elapsed`)?

Note that `--long` does not show all fields.
For example, the following tells us that the given job was above its elapsed time which caused it to be killed.

```
# sacct -j _JOBID_ --format Timelimit,Elapsed
 Timelimit    Elapsed
---------- ----------
  01:00:00   01:00:12
             01:00:13
             01:00:12
```

Use `man sacct`, `sacct --helpformat`, or see the [Slurm Documentation](https://slurm.schedmd.com/sacct.html) for options for the `--format` field of `sacct`.

## I'm getting a "Bus error (core dumped)"

This is most probably caused by your job being allocated insufficient memory.
Please see the memory part of the answer to [My job terminated before being done. What happened?](faq.md#my-job-terminated-before-being-done-what-happened)

## How can I create a new project?

You can create a project if you are either a group leader of an AG or a delegate of an AG.
If this is the case, please follow [these instructions](../../admin/getting-access/#projects).

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

If you see a `Reason=ReqNodeNotAvail,_Reserved_for_maintenance` then also see [Reservations / Maintenances](../../slurm/reservations/).

For GPU jobs also see "My GPU jobs don't get scheduled".

## My GPU jobs don't get scheduled

There are only four GPU machines in the cluster (with four GPUs each, hpc-gpu-1 to hpc-gpu-4).
Please inspect first the number of running jobs with GPU resource requests:

```bash
res-login-1:~$ squeue -o "%.10i %20j %.2t %.5D %.4C %.10m %.16R %.13b" "$@" | grep hpc-gpu- | sort -k7,7
   1902163 ONT-basecalling       R     1    2         8G          hpc-gpu-1   gpu:tesla:2
   1902167 ONT-basecalling       R     1    2         8G          hpc-gpu-1   gpu:tesla:2
   1902164 ONT-basecalling       R     1    2         8G          hpc-gpu-2   gpu:tesla:2
   1902166 ONT-basecalling       R     1    2         8G          hpc-gpu-2   gpu:tesla:2
   1902162 ONT-basecalling       R     1    2         8G          hpc-gpu-3   gpu:tesla:2
   1902165 ONT-basecalling       R     1    2         8G          hpc-gpu-3   gpu:tesla:2
   1785264 bash                  R     1    1         1G          hpc-gpu-4   gpu:tesla:2
```

This indicates that there are two free GPUs on hpc-gpu-4.

Second, inspect the node states:

```bash
res-login-1:~$ sinfo -n hpc-gpu-[1-4]
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
debug*       up    8:00:00      0    n/a
medium       up 7-00:00:00      0    n/a
long         up 28-00:00:0      0    n/a
critical     up 7-00:00:00      0    n/a
highmem      up 14-00:00:0      0    n/a
gpu          up 14-00:00:0      1   drng hpc-gpu-4
gpu          up 14-00:00:0      3    mix med[0301-0303]
mpi          up 14-00:00:0      0    n/a
```

This tells you that hpc-gpu-1 to hpc-gpu-3 have jobs running ("mix" indicates that there are free resources, but these are only CPU cores not GPUs).
hpc-gpu-4 is shown to be in "draining state".
Let's look what's going on there.

```bash hl_lines="10 18"
res-login-1:~$ scontrol show node hpc-gpu-4
NodeName=hpc-gpu-4 Arch=x86_64 CoresPerSocket=16
   CPUAlloc=2 CPUTot=64 CPULoad=1.44
   AvailableFeatures=skylake
   ActiveFeatures=skylake
   Gres=gpu:tesla:4(S:0-1)
   NodeAddr=hpc-gpu-4 NodeHostName=hpc-gpu-4 Version=20.02.0
   OS=Linux 3.10.0-1127.13.1.el7.x86_64 #1 SMP Tue Jun 23 15:46:38 UTC 2020
   RealMemory=385215 AllocMem=1024 FreeMem=347881 Sockets=2 Boards=1
   State=MIXED+DRAIN ThreadsPerCore=2 TmpDisk=0 Weight=1 Owner=N/A MCS_label=N/A
   Partitions=gpu
   BootTime=2020-06-30T20:33:36 SlurmdStartTime=2020-07-01T09:31:51
   CfgTRES=cpu=64,mem=385215M,billing=64
   AllocTRES=cpu=2,mem=1G
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s
   Reason=deep power-off required for PSU [root@2020-07-17T13:21:02]
```

The "State" attribute indicates the node has jobs running but is currenlty being "drained" (accepts no new jobs).
The "Reason" gives that it has been scheduled for power-off for maintenance of the power supply unit.

## When will my job be scheduled?

You can use the `scontrol show job JOBID` command to inspect the scheduling information for your job.
For example, the following job is scheduled to start at `2022-09-19T07:53:29` (`StartTime`) and will be terminated if it does not stop before `2022-09-19T15:53:29` (`EndTime`)
For further information, it has been submitted at `2022-09-15T12:24:57` (`SubmitTime`) and has been last considered by the scheduler at `2022-09-19T07:53:15` (`LastSchedEval`).

```hl_lines="10"
# scontrol show job 4225062
JobId=4225062 JobName=C2371_2
   UserId=user_c(133196) GroupId=hpc-ag-group(1030014) MCS_label=N/A
   Priority=805 Nice=0 Account=hpc-ag-group QOS=normal
   JobState=PENDING Reason=QOSMaxCpuPerUserLimit Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=08:00:00 TimeMin=N/A
   SubmitTime=2022-09-15T12:24:57 EligibleTime=2022-09-15T12:24:57
   AccrueTime=2022-09-15T12:24:57
   StartTime=2022-09-19T07:53:29 EndTime=2022-09-19T15:53:29 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2022-09-19T07:53:15 Scheduler=Main
   Partition=medium AllocNode:Sid=hpc-login-1:557796
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=(null)
   NumNodes=1-1 NumCPUs=25 NumTasks=25 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=25,mem=150G,node=1,billing=25
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=150G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=YES Contiguous=0 Licenses=(null) Network=(null)
   Command=/fast/work/users/user_c/SCZ_replic/JR_sims/GS_wrapy/wrap_y0_VP_2371_GS_chunk2_C02.sh
   WorkDir=/fast/work/users/user_c/SCZ_replic/JR_sims
   StdErr=/fast/work/users/user_c/SCZ_replic/JR_sims/E2371_2.txt
   StdIn=/dev/null
   StdOut=/fast/work/users/user_c/SCZ_replic/JR_sims/slurm-4225062.out
   Power=
```

## My jobs don't run in the partition I expect

You can see the partition that your job runs in with `squeue -j JOBID`:

```bash
res-login-1:~$ squeue -j 877092
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            877092    medium snakejob holtgrem  R       0:05      1 med0626
```

See [Job Scheduler](../overview/job-scheduler.md) for information about the partition's properties and how jbos are routed to partitions.
You can force jobs to run in a particular partition by specifying the `--partition` parameter, e.g., by adding `--partition=medium` or `-p medium` to your `srun` and `sbatch` calls.

## My jobs get killed after four hours

This is probably answered by the answer to [My jobs don't run in the partition I expect](#my-jobs-dont-run-in-the-partition-i-expect).

## How can I mount a network volume from elsewhere on the cluster?

You cannot.
Also see the [For the Impatient](../overview/for-the-impatient.md) section of the manual.

## Why can I not mount a network volume from elsewhere on the cluster?

For performance and stability reasons.
Network volumes are notorious for degrading performance, depending on the used protocol, even stability.

## How can I then access the files from my workstation/server?

You can transfer files to the cluster through Rsync over SSH or through SFTP to the `transfer-1` or `transfer-2` node.

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
    872571 gpu       bash                 user1       R   15:53:25      1 hpc-gpu-3    gpu:tesla:1
    862261 gpu       bash                 user2       R 2-16:26:59      1 hpc-gpu-4    gpu:tesla:4
    860771 gpu       kidney.job           user3       R 2-16:27:12      1 hpc-gpu-2    gpu:tesla:1
    860772 gpu       kidney.job           user3       R 2-16:27:12      1 hpc-gpu-2    gpu:tesla:1
    860773 gpu       kidney.job           user3       R 2-16:27:12      1 hpc-gpu-2    gpu:tesla:1
    860770 gpu       kidney.job           user3       R 4-03:23:08      1 hpc-gpu-1    gpu:tesla:1
    860766 gpu       kidney.job           user3       R 4-03:23:11      1 hpc-gpu-3    gpu:tesla:1
    860767 gpu       kidney.job           user3       R 4-03:23:11      1 hpc-gpu-1    gpu:tesla:1
    860768 gpu       kidney.job           user3       R 4-03:23:11      1 hpc-gpu-1    gpu:tesla:1
```

In the example above, user1 has one job with one GPU running on hpc-gpu-3, user2 has one job running with 4 GPUs on hpc-gpu-4 and user3 has 7 jobs in total running of different machines with one GPU each.

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
res-login-1:~$ ssh med0XXX
```

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

You can specify contraints with OR such as `--constraint=haswell|broadwell|skylake`.
You can see the assignment of architectures to nodes using the `sinfo -o "%8P %.5a %.10l %.6D %.6t %10f %N"` command.
This will also display node partition, availability etc.

## Help, I'm getting a Quota Warning Email!

No worries!

As documented in the [Storage Locations](../storage/storage-locations.md) section, each user/project/group has three storage volumes:
A small `home`, a larger `work` and a large (but temporary) `scratch`.
There are limits on the size of these volumes.
You get a nightly warning email in case you are over the soft limit and you will not be able to write any more data if you get above the hard limit.
When you login to the login nodes, the quotas and current usage is displayed to you.

Please note that not all files will be displayed when using `ls`.
You have to add the `-a` parameter to also show files and directory starting with a dot.
Often, users are confused if these dot directories take up all of their `home` quota.

Use the following command to list **all** files and directories in your home:

```bash
res-login-1:~$ ls -la ~/
```

You can use the following command to see how space each item takes up, including the hidden directories

```bash
res-login-1:~$ du -shc ~/.* ~/* --exclude=.. --exclude=.
```

In the case that, e.g., the `.cpan` directory is large, you can move it to `work` and create a symlink in its original place.

```bash
res-login-1:~$ mv ~/.cpan ~/work/.cpan
res-login-1:~$ ln -sr ~/work/.cpan ~/.cpan
```

## I'm getting a "Disk quota exceeded" error.

Most probably you are running into the same problem as described and solved in the entry [Help, I'm getting a Quota Warning Email!](#help-im-getting-a-quota-warning-email)

## Environment modules don't work and I get "module: command not found"

First of all, ensure that you are on a compute node and not on one of the login nodes.
One common reason is that the system-wide Bash configuration has not been loaded, try to execute `source /etc/bashrc` and then re-try using `module`.
In the case that the problem persists, please contact hpc-helpdesk@bih-charite.de.

## What should my ~/.bashrc look like?

All users get their home directory setup using a skelleton files.
These file names start with a dot `.` and are hidden when you type `ls`, you have to type `ls -a` to see them.
You can find the current skelleton in `/etc/skel.bih` and inspect the content of the Bash related files as follows:

```bash
res-login-1:~$ head /etc/skel.bih/.bash*
==> /etc/skel.bih/.bash_logout <==
# ~/.bash_logout

==> /etc/skel.bih/.bash_profile <==
# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin

==> /etc/skel.bih/.bashrc <==
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=
```

There actually are a couple of more files by default.
The original copy in `/etc/skel.bih` might slightly change over time during improvements but we will not touch your home directory in an unsolicited way at any time!

```bash
res-login-1:~$ tree -a /etc/skel.bih/
/etc/skel.bih/
├── .bash_logout
├── .bash_profile
├── .bashrc
├── .screenrc
└── .vimrc
```

## My program crashes! What should I do?

Have a look at our [How-To: Debug Software](../how-to/misc/debug-software.md) and [How-To: Debug Software on HPC Systems](../how-to/misc/debug-at-hpc.md) guides!

**But it works on my workstation!**

Yes, please also refer to these guides on possible approaches to find the problem.

## My R sessions is always crashing!

We are aware of issues for macOS users that are related to locales.
You might see something like the following in your R session or logs.

Note that this is not macOS specific by itself but by macOS using a locale that R does not understand.

Try using `export LC_ALL=C` before calling R.
If the problem goes away, you can add the following lines to your `~/.bashrc` file to make this permanent.
Note that this is actually a good idea so your locale settings are fixed.

```bash
# Fix all locale settings to "C" for robustness.
export LC_ALL=C
```

!!! note "ChatGPT on Unix Locales"

    Let's se what ChatGPT 3.5 tells us about Unix locale.

    **Me: Please explain the concept of Unix locales to me.**

    > Certainly! Unix locales are a set of environment variables and settings that define how various aspects of software and user interactions should be handled to accommodate different languages, regions, and cultural conventions. They play a crucial role in enabling software to adapt to local preferences for displaying text, formatting numbers, handling dates, and more.
    >
    > In a Unix-like operating system (such as Linux or macOS), locales provide a way to customize how programs behave based on the user's language, country, and cultural norms. Here are some of the key components of Unix locales:
    >
    > 1. **Language:** The primary language used for messages, menus, prompts, and system output. For instance, if you set your locale to "en_US", messages and interfaces will be in English (United States).
    >
    > 2. **Region:** The country or region that defines conventions for aspects like date and time formatting, number formatting (decimal point or comma, thousands separator, etc.), and currency symbols. For example, "en_GB" represents English (Great Britain).
    >
    > 3. **Character Encoding:** Specifies how characters are encoded, which is vital for correctly displaying and handling text in various scripts and languages. Common encodings include UTF-8 for Unicode text.
    >
    > 4. **Collation:** Defines how strings should be sorted. Different languages have different sorting orders for characters, and the locale determines the appropriate collation rules.
    >
    > 5. **Numeric and Date Formats:** Specifies how numbers, dates, and times should be formatted. Different regions have different conventions for date formats (dd/mm/yyyy vs. mm/dd/yyyy) and numerical representations (commas vs. periods as decimal points).
    >
    > 6. **Currency Symbols:** Defines the symbol or abbreviation used for representing currency. This is important for currency formatting in software.
    >
    > Unix locales are typically identified by a combination of language and region codes, like "en_US" for American English or "fr_FR" for French (France). These locales are stored in files on the system, and they contain the necessary rules and data for proper localization.
    >
    > Users and administrators can set the system's default locale, and individual users can override this setting in their environment. Here are some common commands used to manage locales:
    >
    > - `locale`: Displays the current locale settings.
    > - `locale -a`: Lists all available locales on the system.
    > - `export LC_ALL=locale_code`: Temporarily sets the locale for the current session.
    > - `update-locale locale_code`: Changes the system-wide default locale.
    >
    > Programs that are designed to be locale-aware use these settings to display messages, format numbers and dates, and handle other language-dependent features. This allows software to be more user-friendly and relevant to different linguistic and cultural contexts.
    >
    > In summary, Unix locales provide a framework for adapting software behavior to match the language, region, and cultural preferences of users, ensuring a more inclusive and globally accessible computing experience.

## Which CUDA version is installed?

For this, connect to the node you want to query (via SSH but do not perform any computation via SSH!)

```bash
res-login-1:~$ ssh hpc-gpu-1
hpc-gpu-1:~$ yum list installed 2>/dev/null | grep cuda.x86_64
cuda.x86_64                               10.2.89-1                  @local-cuda
nvidia-driver-latest-dkms-cuda.x86_64     3:440.64.00-1.el7          @local-cuda
```

## Can I use Docker on the Cluster?

No, as Docker essentially gives you access as the root user.

However, you can use Apptainer (former Singularity) to run containers (and even many Docker contains if they are "properly built").
Also see [Using Apptainer (with Docker Images)](../how-to/software/apptainer.md).

## How can I copy data between the MAX Cluster (MDC Network) and BIH HPC?

The MAX cluster is the HPC system of the MDC.
It is located in the MDC network.
The BIH HPC is located in the BIH network.

In general, connections can only be initiated from the MDC network to the BIH network.
The reverse does not work.
In other words, you have to log into the MAX cluster and then initiate your file copies to or from the BIH HPC from there.
E.g., use `rsync -avP some/path user_m@hpc-transfer-1.cubi.bihealth.org:/another/path` to copy files from the MAX cluster to BIH HPC and `rsync -avP user_m@hpc-transfer-1.cubi.bihealth.org:/another/path some/path` to copy data from the BIH HPC to the MAX cluster.

## How can I copy data between the Charite Network and BIH HPC?

In general, connections can only be initiated from the Charite network to the BIH network.
The reverse does not work.
In other words, you have to be on a machine inside the Charite network and then initiate your file copies to or from the BIH HPC from there.
E.g., use `rsync -avP some/path user_c@hpc-transfer-1.cubi.bihealth.org:/another/path` to copy files from the MAX cluster to BIH HPC and `rsync -avP user_c@hpc-transfer-1.cubi.bihealth.org:/another/path some/path` to copy data from the BIH HPC to the MAX cluster.

## My jobs are slow/die on the login/transfer node!

As of December 3, 2020 we have established a policy to limit you to 512 files and 128MB of RAM.
Further, you are limited to using the equivalent of one core.
This limit is enforced for all processes originating from an SSH session and the limit is enforced on all jobs.
This was done to prevent users from thrashing the head nodes or using SSH based sessions for computation.

## Slurm complains about `execve` / "No such file or directory"

This means that the program that you want to execute does not exist.
Consider the following example:

```
[user@res-login-1 ~]$ srun --time 2-0 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=12 --mem 96G --partition staging --immediate 5 \
  --pty bash -i
slurmstepd: error: execve(): 5: No such file or directory
srun: error: hpc-cpu-2: task 0: Exited with exit code 2
```

Can you spot the problem?
In this case, the problem is that for long arguments such as `--mem` you **must use the equal sign for `--arg=value`** with Slurm.
This means that instead of writing `--mem 96G --partition staging --immediate 5`, you must use ``--mem=96G --partition=staging --immediate=5`.

In this respect, Slurm deviates from the [GNU argument syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) where the equal sign is optional for long arguments.

## `slurmstepd` says that `hwloc_get_obj_below_by_type` fails

You can ignore the following problem:

```
slurmstepd: error: hwloc_get_obj_below_by_type() failing, task/affinity plugin may be required to address bug fixed in HWLOC version 1.11.5
slurmstepd: error: task[0] unable to set taskset '0x0'
```

This is a minor failure related to Slurm and cgroups.
Your job **should** run through successfully despite this error (that is more of a warning for end-users).

## My login stalls / something weird is happening

You can try to run `ssh -l USER hpc-login-1.cubi.bihealth.org bash -iv`.
This will run `bash -iv` instead of the normal login shell.
The parameter `-i` is creating an interactive shell (which is what you want) and `-v` to see every command that is executed.
This way you will see **every command** that is executed.
You will also be able to identify at which point there is any stalling (e.g., activating conda via `source .../conda` when the fiel system is slow).

## How can I share files/collaborate with users from another work group?

Please use [projects as documented here](../admin/getting-access.md#projects).
Projects were created for this particular purpose.

## What's the relation of Charite, MDC, and cluster accounts?

For HPC 4 Research either an active and working Charite or MDC account is required (that is, you can login e.g., into email.charite.de or mail.mdc-berlin.de).
The system has a separate meta directory that is used for the authorization of users (in other words, whether the user is active, has access to the system, and which groups the user belongs to).
Charite and MDC accounts map to accounts `<Charite user name>_c` and `<MDC user name>_m` accounts in this meta directory.
In the case that a user has both Charite and MDC accounts these are completely separate entities in the meta directory.
For authentication (veryfing that a user has acccess to an account), the Charite and MDC account systems (MS Active Directory) are used.
Authentication currently only uses the SSH keys deposited into Charite (via zugang.charite.de) and MDC (via MDC persdb).
Users have to obtain a suitable Charite/MDC account via Charite and MDC central IT departments and upload their SSH keys through the host organization systems on their own.
The hpc-helpdesk process is then used for getting their accounts setup on the HPC 4 Research system (the home/work/scratch shares being setup), becoming part of the special `hpc-users` group that controls access to the system and organizing users into work groups and projects.

The process of submitting keys to Charite and MDC is documented in the "Connecting" section.

## How do Charite/MDC/Cluster accounts interplay with VPN and the MDC jail node?

Charite users have to obtain a VPN account with the appropriate VPN access permissions, i.e., [Zusatzantrag B as documented here](/bih-cluster/connecting/from-external/#for-charite-users).
For Charite VPN, as for all Charite IT systems, users must use their Charite user name (e.g., `jdoe` and not `jdoe_c`).

MDC users either have to use MDC VPN or the MDC jail node, as [documented here](/bih-cluster/connecting/from-external/#for-mdc-users).
For MDC VPN and jail node, as for all MDC IT systems, users must use their MDC user name (e.g., `jdoe` and not `jdoe_m`).

For help with VPN or jail node, please contact the central Charite or MDC helpdesks as appropriate.

Only when connecting **from** the host organizations' VPN or **from** the host organizations' jail node, the users use the HPC 4 Research user name that is `jdoe_c` or `jdoe_m` and not `jdoe`!

## How can I exchange data with external collaborators?

BIH HPC IT does not have the resources to offer such a service to normal users.

In particular, for **privacy sensitive** data this comes with a large number of strings attached to fulfill all regulatory requirements.
If you need to exchange such data then you need to contact the central IT departments of your home organisation:

- Charite GB IT: heldpesk@charite.de
- MDC: helpdesk@mdc-berlin.de

If your data is not **privacy sensitive** or you can **guarantee strong encryption of the data** then the Gigamove service of RWTH Aachen might come in handy:

- https://gigamove.rwth-aachen.de/en
- https://help.itc.rwth-aachen.de/en/service/1jeqhtat4k0o3/faq/

You can login via Charite/MDC credentials (or most German academic institutions) and store up to 1TB of data at a time in the account with each file having up to 100GB.

As a note, Charite GB IT has a (German) manual on how to use 7-Zip with AES256 and strong passwords for **encrypting data** such that it is fit for transfer over unencrypted channels.
You can find it here (Charite Intranet only) at point 2.12.

- https://intranet.charite.de/it/helpdesk/anleitungen/

The key point is using a strong password (e.g. with the `pwgen` utility), creating an encrypted file with AES256 encryption, using distinct password for each recipient, and exchanging the password over a second channel (SMS or voice phone).
Note that the central manual remains the ground truth of information and this FAQ entry may not reflect the current process recommended by GB IT if it changes without us noticing.
