# Slurm Cheat Sheet

This page contains assorted Slurm commands and Bash snippets that should be helpful.

**`man` pages!**

```bash
$ man sinfo
$ man scontrol
$ man squeue
# etc...
```

**interactive sessions**

```bash
hpc-login-1:~$ srun --pty bash
med0740:~$ echo "Hello World"
med0740:~$ exit
```

**batch submission**

```bash
hpc-login-1:~$ sbatch script.sh
Submitted batch job 2
hpc-login-1:~$ squeue
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

**listing job details**

```
$ scontrol show job 225
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
   Partition=gpu AllocNode:Sid=hpc-login-1:1918
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


```bash
host:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
              1177    medium     bash jweiner_  R 4-21:52:24      1 med0127 
              1192    medium     bash jweiner_  R 4-07:08:40      1 med0127 
              1209   highmem     bash mkuhrin_  R 2-01:07:17      1 med0402 
              1210       gpu     bash hilberta  R 1-10:30:34      1 med0304 
              1213      long     bash schubacm  R 1-09:42:27      1 med0127 
              2401       gpu     bash ramkem_c  R 1-05:14:53      1 med0303 
              2431    medium ngs_mapp holtgrem  R 1-05:01:41      1 med0127 
              2437  critical snakejob holtgrem  R 1-05:01:34      1 med0135 
              2733     debug     bash schubacm  R    7:36:42      1 med0127 
              3029  critical ngs_mapp holtgrem  R    5:59:07      1 med0127 
              3030  critical snakejob holtgrem  R    5:56:23      1 med0134 
              3031  critical snakejob holtgrem  R    5:56:23      1 med0137 
              3032  critical snakejob holtgrem  R    5:56:23      1 med0137 
              3033  critical snakejob holtgrem  R    5:56:23      1 med0138 
              3034  critical snakejob holtgrem  R    5:56:23      1 med0138 
              3035  critical snakejob holtgrem  R    5:56:20      1 med0139 
              3036  critical snakejob holtgrem  R    5:56:20      1 med0139 
              3037  critical snakejob holtgrem  R    5:56:20      1 med0140 
              3038  critical snakejob holtgrem  R    5:56:20      1 med0140 
              3039  critical snakejob holtgrem  R    5:56:20      1 med0141 
              3040  critical snakejob holtgrem  R    5:56:20      1 med0141 
              3041  critical snakejob holtgrem  R    5:56:20      1 med0142 
              3042  critical snakejob holtgrem  R    5:56:20      1 med0142 
              3043  critical snakejob holtgrem  R    5:56:20      1 med0143 
              3044  critical snakejob holtgrem  R    5:56:20      1 med0143 
              3063      long     bash schubacm  R    4:12:37      1 med0127 
              3066      long     bash schubacm  R    4:11:47      1 med0127 
              3113    medium ngs_mapp holtgrem  R    1:52:33      1 med0708 
              3118    medium snakejob holtgrem  R    1:50:38      1 med0133 
              3119    medium snakejob holtgrem  R    1:50:38      1 med0703 
              3126    medium snakejob holtgrem  R    1:50:38      1 med0706 
              3127    medium snakejob holtgrem  R    1:50:38      1 med0144 
              3128    medium snakejob holtgrem  R    1:50:38      1 med0144 
              3133    medium snakejob holtgrem  R    1:50:35      1 med0147 
              3134    medium snakejob holtgrem  R    1:50:35      1 med0147 
              3135    medium snakejob holtgrem  R    1:50:35      1 med0148 
              3136    medium snakejob holtgrem  R    1:50:35      1 med0148 
              3138    medium snakejob holtgrem  R    1:50:35      1 med0104 
```

```bash
host:~$ squeue -o "%.10i %9P %20j %10u %.2t %.10M %.6D %10R %b"
     JOBID PARTITION NAME                 USER       ST       TIME  NODES NODELIST(R TRES_PER_NODE
      1177 medium    bash                 jweiner_m   R 4-21:52:22      1 med0127    N/A
      1192 medium    bash                 jweiner_m   R 4-07:08:38      1 med0127    N/A
      1209 highmem   bash                 mkuhrin_m   R 2-01:07:15      1 med0402    N/A
      1210 gpu       bash                 hilberta_c  R 1-10:30:32      1 med0304    gpu:tesla:4
      1213 long      bash                 schubacm_c  R 1-09:42:25      1 med0127    N/A
      2401 gpu       bash                 ramkem_c    R 1-05:14:51      1 med0303    gpu:tesla:1
      2431 medium    ngs_mapping          holtgrem_c  R 1-05:01:39      1 med0127    N/A
      2437 critical  snakejob.ngs_mapping holtgrem_c  R 1-05:01:32      1 med0135    N/A
      2733 debug     bash                 schubacm_c  R    7:36:40      1 med0127    N/A
      3029 critical  ngs_mapping          holtgrem_c  R    5:59:05      1 med0127    N/A
      3030 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:21      1 med0134    N/A
      3031 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:21      1 med0137    N/A
      3032 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:21      1 med0137    N/A
      3033 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:21      1 med0138    N/A
      3034 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:21      1 med0138    N/A
      3035 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0139    N/A
      3036 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0139    N/A
      3037 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0140    N/A
      3038 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0140    N/A
      3039 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0141    N/A
      3040 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0141    N/A
      3041 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0142    N/A
      3042 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0142    N/A
      3043 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0143    N/A
      3044 critical  snakejob.ngs_mapping holtgrem_c  R    5:56:18      1 med0143    N/A
      3063 long      bash                 schubacm_c  R    4:12:35      1 med0127    N/A
      3066 long      bash                 schubacm_c  R    4:11:45      1 med0127    N/A
      3113 medium    ngs_mapping          holtgrem_c  R    1:52:31      1 med0708    N/A
      3118 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:36      1 med0133    N/A
      3119 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:36      1 med0703    N/A
      3126 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:36      1 med0706    N/A
      3127 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:36      1 med0144    N/A
      3128 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:36      1 med0144    N/A
      3133 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:33      1 med0147    N/A
      3134 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:33      1 med0147    N/A
      3135 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:33      1 med0148    N/A
      3136 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:33      1 med0148    N/A
      3138 medium    snakejob.ngs_mapping holtgrem_c  R    1:50:33      1 med0104    N/A
```


```bash
host:~$ sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
debug*       up    8:00:00     11  drain med[0707,0709-0710,0740-0742,0744-0745,0749,0752,0755] 
debug*       up    8:00:00      8    mix med[0104,0127,0133-0135,0703,0706,0708] 
debug*       up    8:00:00     10  alloc med[0137-0144,0147-0148] 
debug*       up    8:00:00    103   idle med[0105-0124,0136,0145-0146,0151-0164,0201-0264,0704-0705] 
medium       up 7-00:00:00     11  drain med[0707,0709-0710,0740-0742,0744-0745,0749,0752,0755] 
medium       up 7-00:00:00      8    mix med[0104,0127,0133-0135,0703,0706,0708] 
medium       up 7-00:00:00     10  alloc med[0137-0144,0147-0148] 
medium       up 7-00:00:00    103   idle med[0105-0124,0136,0145-0146,0151-0164,0201-0264,0704-0705] 
long         up 28-00:00:0     11  drain med[0707,0709-0710,0740-0742,0744-0745,0749,0752,0755] 
long         up 28-00:00:0      8    mix med[0104,0127,0133-0135,0703,0706,0708] 
long         up 28-00:00:0     10  alloc med[0137-0144,0147-0148] 
long         up 28-00:00:0    103   idle med[0105-0124,0136,0145-0146,0151-0164,0201-0264,0704-0705] 
critical     up 7-00:00:00     11  drain med[0707,0709-0710,0740-0742,0744-0745,0749,0752,0755] 
critical     up 7-00:00:00      8    mix med[0104,0127,0133-0135,0703,0706,0708] 
critical     up 7-00:00:00     10  alloc med[0137-0144,0147-0148] 
critical     up 7-00:00:00    103   idle med[0105-0124,0136,0145-0146,0151-0164,0201-0264,0704-0705] 
highmem      up 14-00:00:0      1    mix med0402 
highmem      up 14-00:00:0      3   idle med[0401,0403-0404] 
gpu          up 14-00:00:0      2    mix med[0303-0304] 
gpu          up 14-00:00:0      2   idle med[0301-0302] 
```