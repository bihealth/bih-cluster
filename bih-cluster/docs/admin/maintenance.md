# Next Maintenance Window

This page documents the current and known upcoming maintenance windows.

## Login & Compute maintenance, October 22nd, 2024

Kernel and Slurm upgrades. Keep an eye on [our forum](https://hpc-talk.cubi.bihealth.org/t/maintenance-slurm-upgrade-october-22nd/1095) for updates. 

## Login, Compute and Storage Maintenance, December 13-14, 2022

All informationand updates regarding maintenance will be circulated on our forum https://hpc-talk.cubi.bihealth.org/c/announcements/5.

## Login, Compute and Storage Maintenance, March 22-23, 2022

All COMPUTE nodes and STORAGE resources won't be reachable!
  
> All nodes will be running in RESERVATION mode. This means you are still able to schedule new jobs on these nodes if their potential/allowed runtime does not extend into the maintenance window (Tuesday and Wednesday, March 22 and 23, all-day). For example, if you submit a job that can run up to 7 days after March 15 then the job will remain in "pending/PD" state giving the explanation of "all nodes being reserved or unavailable".

Issues of today's maintenance:

- Mounting of storage to `/tmp` on login nodes
- Changing mount options of the root partition on the compute nodes
- Upgrading all nodes kernels and further packages
    - This implies an upgrade of CUDA, Singularity, and further packages
- Cold reboot ("power off, power on") of storage system
- Exchanging `cephfs-2` switches (Tier 2 storage, not relevant for most users)

**IMPORTANT**

- **All nodes will reboot**
- **All running jobs will die**
- **All sessions on login nodes will die**

[Progress Thread on hpc-talk](https://hpc-talk.cubi.bihealth.org/t/maintenance-2022-03-22-progress-thread/90)

## DRMAA Deprecation, March 2, 2022

- The usage of DRMAA on the HPC is deprecated.
- In Snakemake, it has been deprecated in favor of using Snakemake Profiles [as documented](../slurm/snakemake.md#snakemake-and-slurm).
- We will support DRMAA at least until June 30, 2022 but ask all users to migrate away from it as soon as possible.
- Background:
    - With DRMAA, the status of each job is queried for using `scontrol show job JOBID` and `sacct -j JOBID`.
    - This leads to regular remote procedure calls (RPC) to the slurm control daemon.
    - It leads to **a lot of such calls**.
    - It leads to so many calls that it prevents the scheduler from working correctly and leads to service degradation for all users.
- Using Snakemake profiles is easy.
    - Call Snakemake with `snakemake --profile=cubi-v1` instead of `snakemake --drmaa "..."`.
    - In your rules, specify threads, running time and memory as:
      ```
      rule myrule:
        # ...
        threads: 8
        resources:
          time="12:00:00",
          memory="8G",
        # ...
      ```

## Cluster Setting Tuning, March 1, 2022

- We have adjusted the scheduler settings to address high number of jobs by users:
- `SchedulerParameters+=bf_max_job_user=50`: backfill scheduler only considers 50 jobs of each user. This mitigates an issue with some users having too many jobs and thus other users' jobs don't get ahead in the queue
- `EnforcePartLimits=ALL`: jobs that don't fit into their partition are rejected
- `DependencyParameters=kill_invalid_depend`: jobs that have dependencies set that cannot be fulfilled will be killed

## Limiting Global Memory Usage, February 14, 2022

- A global memory allocation limit per user is set per partition.
- The value is set to "max CPU count per user * 7GB".
- Users can allocate up to "max cpu count" CPUs or "max cpu count * 7GB" RAM.
- This is enforced globally (users could allocate 2/3 of their global CPU limit with 3.5 GB RAM and 1/3 with 7GB of RAM, for example).

## Ganglia Fixes & Docs, February 3, 2022

- Reparing GPFS and NVIDIA GPU monitoring in [Ganglia](https://hpc-ganglia.cubi.bihealth.org)
- Root cause was that the Python modules in Ganglia were removed from EPEL.
  We now have a local package build of Ganglia, if you are interested, here is the [patch](https://github.com/bihealth/rpms-ganglia) and [Docker based build instructions](https://github.com/bihealth/rpms-ganglia/issues/1).
- You can find some [documentation about our Ganglia here](../overview/monitoring.md).

## Misc Changes, January 29, 2022

- We have reduced oversubscription to 2x from 4x.
- We have setup the user quota on /tmp on the login nodes to 20MB to improve stability of the nodes.

## Enabling Oversubscription, January 6, 2022

- Many resources remain unused as users allocate too many cores to their jobs.
  Slurm will now oversubscribe jobs in terms of CPUs, i.e., schedule more than one allocated core per physical core/thread.

## Enforcing Usage of `localtmp` Resource, January 31, 2022

- We will enforce using `localtmp` resource for local storage above 100MB.
- See [Slurm: Temporary Files](../slurm/temporary-files.md) for details.

## Temporary File Handling Changes, December 27, 2021

- Each job gets its private `/tmp` using Linux namespaces/cgroups.
  This greatly improves the reliability of cleaning up after jobs.
  (Technically, this is implemented using the Slurm [job_container/tmpfs](https://slurm.schedmd.com/job_container.conf.html)) plugin.
- We are starting to track available local temporary space with Slurm in the general resource (`Gres`) "localtmp".
  In the future this will become a requirement.
  Also see [Slurm: Temporary Files](../slurm/temporary-files.md).

## Cluster Node Upgrades, December 22-23, 2021

- Renaming of cluster head nodes to:
  - `hpc-login-1.cubi.bihealth.org`
  - `hpc-login-2.cubi.bihealth.org`
  - `hpc-portal.cubi.bihealth.org`
  - `hpc-transfer-1.cubi.bihealth.org`
  - `hpc-transfer-2.cubi.bihealth.org`
- Upgraded cluster operating system from CentOS 7.9 to Rocky Linux 8.5.
- Added three more GPU nodes with Tesla V100 GPUS: `hpc-gpu-{5..7}`.  
- Slurm has been upgraded to `28.08.5`.
- Ganglia monitoring generally available at https://hpc-ganglia.cubi.bihealth.org, from internal networks.
- We have applied a number of changes to maximal running times in Slurm configuration.

## GPFS Upgrade, December 20-21, 2021

The GPFS storage system has been upgraded to the latest version to make compatible with Enterprise Linux version 8.

## Slurm upgrade to `21.08.0`, September 8, 2021

Slurm has been upgraded to version `21.08.0`.

## Network re-cabling, September 7-8, 2021

All servers/nodes won't be reachable!

>All nodes will be running in reservation mode. This means you are still able to schedule new jobs on these nodes if their potential/allowed runtime does not extend into the maintenance window (Tuesday and Wednesday, September 7 and 8, all-day). For example, if you submit a job that can run up to 7 days after August 30 then the job will remain in "pending/PD" state giving the explanation of "all nodes being reserved or unavailable".

>If you already have a job running on any nodes that goes beyond September 7, 12:00 am (00:00 Uhr), this job will die.

## Renaming of GPU & High Memory Machines & Scheduler Changes, September 7, 2021

The GPU machines `med030[1-4]` have been renamed to `hpc-gpu-[1-4]`.
The high memory machines `med040[1-4]` have been renamed to `hpc-mem-[1-4]`.
It will probably take us some time to update all places in the documentation.

Further, the `long` partition has been changed to allow jobs with a maximum running time of 14 days.

## New Nodes in the `staging` partition, August 31, 2021

We have installed 36 new nodes (in BETA mode) in the cluster called `hpc-node-[1-36]`.
They have 48 cores (thus 96 hardware threads) each and have 360GiB of main memory available (for the hardware nerds, it's Intel(R) Xeon(R) Gold 6240R CPUs at 2.40GHz, featuring the `cascadelake` architecture).

Right now, they are only available in the `staging` partition.
After some testing we will move them to the other partitions.
We'd like to ask you to test them as well and report any issues to hpc-helpdesk@bih-charite.de.
The nodes have been setup identically to the existing `med0xxx` nodes.
We do not expect big changes but the nodes might not be as stable as other oness.

Here is how you can reach them.

```
hpc-login-1 # srun --immediate=5 --pty --time=24:00:00 --partition=staging bash -i
[...]
hpc-cpu-1 #
```

Note that I'm specifying a maximal running time of 24h so the scheduler will end the job after 24 hours which is before the upcoming maintenance reservation begins.
By default, the scheduler allocates 28 days to the job which means that the job cannot end before the reservation and will be scheduled to start **after** it.
See [Reservations / Maintenances](../slurm/reservations.md) for more information about maintenance reservations.

## Reservation / Maintenance Display on Login, August 30, 2021

User will now be notified on login about maintenance, for example:

```
NOTE: scheduled maintenance(s)

 1: 2021-09-07 00:00:00 to 2021-09-09 00:00:00 ALL nodes

Slurm jobs will only start if they do not overlap with scheduled reservations.
More information:

  - https://bihealth.github.io/bih-cluster/slurm/reservations/
  - https://bihealth.github.io/bih-cluster/admin/maintenance/
```

## Update to Job Sumission Script, August 23, 2021

The `srun` command will now behave as if `--immediate=60` has been specified by default.
It explains how to override this behaviour and possible reasons for job scheduling to fail within 60 seconds (reservations and full cluster).

## Slurm upgrade, August 6, 2021

We upgrade from `20.11.2` to `20.11.8` which contains some fixes for bugs that our users actually stumbled over. The change should be non-intrusive as it's only a patch-level update.

## Networking hardware exchange, August 3, 2021

Following servers won't be reachable:

- GPU nodes (med03xx)
- computing nodes (med0233-0248)

> These nodes are running in reservation mode now. This means you are still able to schedule new jobs on these nodes if their potential/allowed runtime does not extend into the maintenance window (Tuesday, August 3, all-day). For example, if you submit a job that can run up to 7 days after July 26 then the job will remain in "pending/PD" state giving the explanation of "all nodes being reserved or unavailable". If you have a job running on any of the before mentioned nodes that goes beyond August 3, 12:00 am (00:00 Uhr), this job will die. We do not expect the remaining nodes to be affected. However, there remains a minor risk of unexpected downtime of other nodes.

## Server reorganization, July 13, 2021

Affected servers are:

- med02xx
- med07xx

## Server reorganization, June 22 + 23, 2021

> If you have a job running on any of the before mentioned nodes that goes beyond June 22, 6am, this job will die.
> We put a so-called Slurm reservation for the maintenance period.
> Any job that is scheduled before the maintenance and whose end time (start time + max running time) is not before the start of the maintenance will not be scheduled with the message ReqNodeNotAvail, Reserved for maintenance.

Affected servers are:

- med01xx
- med05xx
- med06xx
- med03xx
- med0405

## Memory and PSU exchange, May 31, 2021

- Memory exchange
  - transfer-1.research (OK)
  - med0143 (OK)
  - med0147 (OK)
  - med0206 (FAIL - exchange part broken)
  - med0233 (FAIL - exchange part broken)
  - med0254 (FAIL - exchange part broken)
- PSU exchange
  - med-host024 (OK)

## Moving servers, May 20 + May 25, 2021

- Physically moving proxmox-{2,4} and transfer-2.research (May 20)
- Physically moving proxmox-{1,3} and transfer-1.research (May 25)

## Miscellaneous Maintenances, December 23-25, 2020

**HPC 4 Research**

- Separate HPC 4 Research group GID space from other organization's.
    - **Fully Unavailable**
- Reboot login nodes to increase RAM on hpc-login-2.research
- Update firmwares of transfer-{1,2}.research

## CentOS 8 Migration (in planning)

!!! note

    :construction: This task is currently being planned.
    No schedule has been fixed yet. :construction:

- All nodes will be upgraded to CentOS 8.
- This will be done in a rolling fashion over the course of 1 month.
- The login nodes must be rebooted which we will do with a break of 2 days (one node will remain running).

## Finalize unification of Mass Data Mounts

!!! note

    :construction: This task is currently being planned.
    No schedule has been fixed yet. :construction:

- We will remove the bind mount `/fast` that currently points to `/data/gpfs-1` on HPC 4 Research.
- Users should use `/data` instead of `/fast` everywhere, e.g., `/data/users/$NAME` etc.

# Previous Maintenance Windows

## HPC 4 Research: Miscellaneous Maintenances, December 1, 2020

**Time**: 6am-12am

- Exchange GPFS Controller
    - We need to exchange a central piece of hardware in the storage system.
    - We do not expect a downtime, only a degradation of service.
    - **Access to the GPFS will be degraded**
- Slurm Scheduler
    - Upgrade to the latest and greatest version.
    - Restructure scheduler installation to ease rolling upgrades without future downtimes.
    - Archival of old accounting information to improve schedule performance.
    - **Slurm will be unavailable.**
- Re-Mounting of GPFS
    - The `/fast` file system will be re-mounted to `/data/gpfs-1`.
    - `/fast` becomes a symbolic link to `/data` on all of the cluster.
    - **GPFS access will disappear for some time.**
- Login & Transfer Node Migration
    - The login nodes will be moved from physical machines to virtual machines in high-availability mode.
    - Further, they will be available as `hpc-login-1.cubi.bihealth.org` and `login-2...` instead of  `hpc-login-{1,2}`.
    - The same is true for, `hpc-transfer-{1,2}` which will be replaced by `transfer-1.research.hpc.bihealth.org` and `transfer-2...`.
    - The aim is to improve stability and make everything easier to manage by administration.

### Current Status / Result

- We had to clear the accounting information database to make the update work within an acceptable time (we have 4M+ jobs in there).
  From now on we will only keep the last 31 days in the database (updated nightly).
- The old login and transfer nodes have been made available as nodes `med010[1-3]` and `med012[5-6]`.
- All nodes are available again.
- The maintenance is complete.
## Slurm Scheduler Updates: September 8, 2020

- To improve the scheduling behaviour we will need to restart the Slurm scheduler at ~8am.
- If everything runs well, this will finish after 30minutes (8:30 am).
- Planned Scheduler Changes:
    - Introduce automatic routing of jobs to partitions.
    - Make Slurm scheduler and accounting run more robustly.

## Network Maintenance: June 3, 2020

On June 3, we need to perform a network maintenance at 8 am.

If everything goes well, there might be a short delay in network packages and connections will survive.
In this case, the maintenance will end 8:30 am.

Otherwise, the maintenance will finish by noon.

## Cluster Maintenance with Downtime: June 16

We need to schedule a full cluster downtime on June 16.

## Slurm Migration

We will switch to the [Slurm](https://slurm.schedmd.com/) workload scheduler (from the legacy SGE).
The main reason is that Slurm allows for better scheduling of GPUs (and has loads of improvements over SGE), but the syntax is a bit different.
Currently, our documentation is in an transient state.
We are currently extending our [Slurm-specific documentation](../slurm/quickstart.md).

- **March 7, 2020 (test stage)**:
  Slurm will provide 16 CPU and 3 GPU nodes (with 4 Tesla V100 each), and two high memory nodes, the remaining nodes are available in SGE.
  We ask users to look into scheduling with Slurm.
- **March 31, 2020 (intermediate stage)**:
  Half of the nodes will be migrated to the Slurm cluster (~100), all high memory and GPU nodes will be moved to Slurm.
  New users are advised to use not learn SGE any more but directly use Slurm.
  Support for SGE is limited to bug fixing only (documentation and tips are phased out).
- **May 31, 2020 (sunsetting SGE)**:
  All but 16 nodes will remain in the SGE cluster.
- **June 31, 2020 (the end)**:
  SGE has reached its end of life on hpc4research.

## SSH Key Management

SSH Key Management has switched to using Charite and MDC ActiveDirectory servers.
**You need to upload all keys by the end of April 2020.**

- [MDC Key Upload](../connecting/submit-key/mdc.md)
- [Charite Key Upload](../connecting/submit-key/charite.md)

**Schedule**

- `Feb 4, 2020:` Keys are now also taken from central MDC/Charite servers.
  *You do not need to contact us any more to update your keys (we cannot accelerate the process at MDC).*
- `May 1, 2020:` Keys are now **only** taken from central MDC/Charite servers.
  **You must upload your keys to central servers by then.**

## Switch update, Location Flip of hpc-login-2 and hpc-transfer-1

- Monday, February 23, 9am-15am.

Affected systems:

- `hpc-transfer-1`
- `hpc-transfer-2`
- `hpc-login-2`
- a few compute nodes

The compute nodes are non-critical as we are taking them out of the queues now.

## CentOS 7.6 Upgrade, January 29, February 5

- Wednesday, January 29, 2018: **Reboot hpc-login-1, hpc-transfer-1**
- Wednesday, February 5, 2018: **Reboot hpc-login-2, hpc-transfer-2**

## September 03-30, 2018

Starting monday 03.09.2018 we will be performing rolling update of the cluster from CentOS 7.4 to CentOS 7.5.
Since update will be performed in small bunches of nodes, the only impact you should notice is smaller number of nodes available for computation.

Also, for around two weeks, you can expect that your jobs can hit both CentOS 7.4 & CentOS 7.5 nodes. This should not impact you in any way, but if
you encounter any unexpected behavior of the cluster during this time, please let us know.

At some point we will have to update the transfer, and login nodes. We will do this also in parts, so the you can switch to the other machine.

Key dates are:

**18.09.2018** - hpc-login-1 & hpc-transfer-1 will not be available, and you should switch to hpc-login-2  & hpc-transfer-2 respectively.

**25.09.2018** - hpc-login-2 & hpc-transfer-2 will not be available, and you should switch to hpc-login-1  & hpc-transfer-1 respectively.

Please also be informed that non-invasive maintenance this weekend which we announced has been canceled, so cluster will operate normally.

In case of any concerns, issues, do not hesitate to contact us via hpc-admin@bih-charite.de, or hpc-helpdesk@bih-charite.de.


## June 18, 2018, 0600-1500

Due to tasks we need to perform on BIH cluster, we have planned maintenance:

- Maintenance start: 18.06.2018 06:00 AM
- Maintenance end: 18.06.2018 3:00 PM

During maintenance we will perform several actions:

- GPFS drives re-balancing to improve performance
- ~~OS update on cluster, transfer, and login nodes~~

During maintenance whole cluster will not be usable, this includes:

- you will not be able to run jobs on cluster (SGE queuing system will be shutdown)
- hpc-login-{1,2} nodes will not work reliably during this time
- hpc-transfer-{1-2} nodes, and resources shared by them will be not available

Maintenance window is quite long, since we are dependent on external vendor. However, we will recover services as soon as possible.

We will keep you posted during maintenance with services status.

## March 16-18, 2018 (MDC IT)

**MDC IT** has a network maintenance from **Friday, March 16 18:00 hours until Sunday March 18 18:00 hours.**

This will affect connections to the cluster but no connections within the cluster.

## January 17, 2018 (Complete)

**STATUS:** complete


The first aim of this window is to upgrade the cluster to CentOS 7.4 to patch against the Meltdown/Spectre vulnerabilities. For this, the login and transfer nodes have to be rebooted.

The second aim of this window is to reboot the file server to mitigate some NFS errors. For this, the SGE master has to be stopped for some time.

### Plan/Progress

- [x] reboot med-file1
- [ ] update to CentOS 7.4
    - [x] front nodes
        - [x] hpc-login-1
        - [x] hpc-login-2
        - [x] hpc-login-3 (admin use only)
        - [x] hpc-transfer-1
        - [x] hpc-transfer-2
    - [x] infrastructure nodes
        - [x] qmaster*
        - [ ] ~~install-srv~~
    - [ ] compute nodes
        - [x] med0100 to med0246
        - [ ] med0247 to med0764
    - [x] special purpose compute nodes
        - [x] med0401 (high-memory)
        - [x] med0402 (high-memory)
        - [x] med0403 (high-memory)
        - [x] med0404 (high-memory)
        - [x] med0405 (GPU)

## Previous Maintenance

(since January 2010)

- none
