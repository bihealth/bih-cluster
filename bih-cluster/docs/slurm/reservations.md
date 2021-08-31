# Reservations / Maintenances


!!! hint
    Read this in particular if you want to know why your job does not get scheduled and you see `Reason=ReqNodeNotAvail,_Reserved_for_maintenance` in `scontrol show job` :scream:.

Administration registers maintenances with the Slurm scheduler as so-called **reservations**.
You can see the current reservations with `scontrol show reservation`.
The following is a scheduled reservation affecting ALL nodes of the cluster.

```
# scontrol show reservation
ReservationName=root_13 StartTime=2021-09-07T00:00:00 EndTime=2021-09-09T00:00:00 Duration=2-00:00:00
   Nodes=hpc-cpu-[1-36],med[0101-0116,0201-0264,0301-0304,0401-0404,0501-0516,0601-0632,0701-0764]
   NodeCnt=236 CoreCnt=5344 Features=(null) PartitionName=(null)
   Flags=MAINT,IGNORE_JOBS,SPEC_NODES,ALL_NODES TRES=cpu=10176
   Users=root Groups=(null) Accounts=(null) Licenses=(null) State=INACTIVE BurstBuffer=(null) Watts=n/a
   MaxStartDelay=(null)
```

You will also be notified when logging into the login nodes, e.g.,

```
--
             ***NOTE: 1 scheduled maintenance(s)***

 1: 2021-09-07 00:00:00 to 2021-09-09 00:00:00 ALL nodes

You jobs do not start because of "Reserved_for_maintenance"?
Slurm jobs will only start if they do not overlap with scheduled reservations.
More information:

  - https://bihealth.github.io/bih-cluster/slurm/reservations/
  - https://bihealth.github.io/bih-cluster/admin/maintenance/
--
```

## What is the Effect of a Reservation?

Maintenance reservations will block the affected nodes (or even the whole cluster) for jobs.
If there is a maintenance in one week then your **job must have an end time before the reservation starts**.
By this, the job gives a guarantee to the scheduler that it will not interfer with the maintenance reservation.

For example, `scontrol show job JOBID` might report the following

```
JobId=4011580 JobName=snakejob
   UserId=USER(UID) GroupId=GROUP(GID) MCS_label=N/A
   Priority=1722 Nice=0 Account=GROUP QOS=normal
   JobState=PENDING Reason=ReqNodeNotAvail,_Reserved_for_maintenance Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=28-00:00:00 TimeMin=N/A
   SubmitTime=2021-08-30T09:01:01 EligibleTime=2021-08-30T09:01:01
   AccrueTime=2021-08-30T09:01:01
   StartTime=2021-09-09T00:00:00 EndTime=2021-10-07T00:00:00 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2021-08-30T10:20:40
   Partition=long AllocNode:Sid=172.16.35.153:5453
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=(null)
   NumNodes=1-1 NumCPUs=8 NumTasks=8 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=8,mem=4G,node=1,billing=8
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=4G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Power=
   NtasksPerTRES:0
```

Look out for the `Reason` line:

```
Reason=ReqNodeNotAvail,_Reserved_for_maintenance
```

This job is scheduled to run up to 4 weeks and has been submitted on 2021-08-30.

Right now the following reservation is active

```
# scontrol show reservation
ReservationName=root_13 StartTime=2021-09-07T00:00:00 EndTime=2021-09-09T00:00:00 Duration=2-00:00:00
   Nodes=hpc-cpu-[1-36],med[0101-0116,0201-0264,0301-0304,0401-0404,0501-0516,0601-0632,0701-0764]
   NodeCnt=236 CoreCnt=5344 Features=(null) PartitionName=(null)
   Flags=MAINT,IGNORE_JOBS,SPEC_NODES,ALL_NODES TRES=cpu=10176
   Users=root Groups=(null) Accounts=(null) Licenses=(null) State=INACTIVE BurstBuffer=(null) Watts=n/a
   MaxStartDelay=(null)
```

Thus, the scheduler decided to set a `StartTime` of the job to `2021-09-09T00:00:00`, which is the end time of the reservation.
Effectively, the job is forced to run outside the maintenance reservation.

You can resolve this by using a `--time=` parameter to `srun` or `sbatch` such that the job ends before the maintenance reservation starts.
