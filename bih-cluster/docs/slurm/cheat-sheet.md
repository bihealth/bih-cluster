# Slurm Cheat Sheet

This page contains assorted Slurm commands and Bash snippets that should be helpful.

**`man` pages!**

```bash
$ man sinfo
$ man scontrol
$ man squeue
# etc...
```

interactive sessions

```bash
med-login1:~$ srun --pty bash
med0740:~$ echo "Hello World"
med0740:~$ exit
```

batch submission

```bash
med-login1:~$ sbatch script.sh
Submitted batch job 2
med-login1:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                27     debug script.s holtgrem  R       0:06      1 med0703
```

**listing nodes**

```bash
$ sinfo -N
NODELIST   NODES PARTITION STATE
med0740        1    debug* idle
med0741        1    debug* down*
med0742        1    debug* down*

$ scontrol show nodes
NodeName=med0740 Arch=x86_64 CoresPerSocket=8
   CPUAlloc=0 CPUTot=32 CPULoad=0.06
   AvailableFeatures=(null)
[...]

$ scontrol show nodes med0740
NodeName=med0740 Arch=x86_64 CoresPerSocket=8
   CPUAlloc=0 CPUTot=32 CPULoad=0.06
   AvailableFeatures=(null)
   ActiveFeatures=(null)
   Gres=(null)
   NodeAddr=med0740 NodeHostName=med0740 Version=20.02.0
   OS=Linux 3.10.0-1062.12.1.el7.x86_64 #1 SMP Tue Feb 4 23:02:59 UTC 2020
   RealMemory=1 AllocMem=0 FreeMem=174388 Sockets=2 Boards=1
   State=IDLE ThreadsPerCore=2 TmpDisk=0 Weight=1 Owner=N/A MCS_label=N/A
   Partitions=debug
   BootTime=2020-03-05T00:54:15 SlurmdStartTime=2020-03-05T16:23:25
   CfgTRES=cpu=32,mem=1M,billing=32
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s
```

**queue states**

```bash
$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
$ squeue -u holtgrem_c
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

**node resources**

```bash
$ sinfo -o "%20N  %10c  %10m  %25f  %10G "
```

**additional resources such as GPUs**

```bash
$ sinfo -o "%N %G"
```

**listing job details
```
$ # scontrol show job 225
JobId=225 JobName=bash
   UserId=XXX(135001) GroupId=XXX(30069) MCS_label=N/A
   Priority=4294901580 Nice=0 Account=(null) QOS=normal
   JobState=FAILED Reason=NonZeroExitCode Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=0 Reboot=0 ExitCode=130:0
   RunTime=00:16:27 TimeLimit=14-00:00:00 TimeMin=N/A
   SubmitTime=2020-03-23T11:34:26 EligibleTime=2020-03-23T11:34:26
   AccrueTime=Unknown
   StartTime=2020-03-23T11:34:26 EndTime=2020-03-23T11:50:53 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2020-03-23T11:34:26
   Partition=gpu AllocNode:Sid=med-login1:1918
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=med0301
   BatchHost=med0301
   NumNodes=1 NumCPUs=2 NumTasks=0 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,node=1,billing=2
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=bash
   WorkDir=XXX
   Power=
   TresPerNode=gpu:tesla:4
   MailUser=(null) MailType=NONE
```
