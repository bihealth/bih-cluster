# Next Maintenance Window

This page documents the current and known upcoming maintenance windows.

## Miscellaneous Maintenances, December 23-25, 2020

**HPC 4 Research**

- Separate HPC 4 Research group GID space from other organization's.
    - **Fully Unavailable**
- Reboot login nodes to increase RAM on login-2.research
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
    - Further, they will be available as `login-1.research.hpc.bihealth.org` and `login-2...` instead of  `med-login{1,2}`.
    - The same is true for, `med-transfer{1,2}` which will be replaced by `transfer-1.research.hpc.bihealth.org` and `transfer-2...`.
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

## Switch update, Location Flip of med-login2 and med-transfer1

- Monday, February 23, 9am-15am.

Affected systems:

- `med-transfer1`
- `med-transfer2`
- `med-login2`
- a few compute nodes

The compute nodes are non-critical as we are taking them out of the queues now.

## CentOS 7.6 Upgrade, January 29, February 5

- Wednesday, January 29, 2018: **Reboot med-login1, med-transfer1**
- Wednesday, February 5, 2018: **Reboot med-login2, med-transfer2**

## September 03-30, 2018

Starting monday 03.09.2018 we will be performing rolling update of the cluster from CentOS 7.4 to CentOS 7.5.
Since update will be performed in small bunches of nodes, the only impact you should notice is smaller number of nodes available for computation.

Also, for around two weeks, you can expect that your jobs can hit both CentOS 7.4 & CentOS 7.5 nodes. This should not impact you in any way, but if
you encounter any unexpected behavior of the cluster during this time, please let us know.

At some point we will have to update the transfer, and login nodes. We will do this also in parts, so the you can switch to the other machine.

Key dates are:

**18.09.2018** - med-login1 & med-transfer1 will not be available, and you should switch to med-login2  & med-transfer2 respectively.

**25.09.2018** - med-login2 & med-transfer2 will not be available, and you should switch to med-login1  & med-transfer1 respectively.

Please also be informed that non-invasive maintenance this weekend which we announced has been canceled, so cluster will operate normally.

In case of any concerns, issues, do not hesitate to contact us via hpc-admin@bihealth.de, or hpc-helpdesk@bihealth.de.


## June 18, 2018, 0600-1500

Due to tasks we need to perform on BIH cluster, we have planned maintenance:

- Maintenance start: 18.06.2018 06:00 AM
- Maintenance end: 18.06.2018 3:00 PM

During maintenance we will perform several actions:

- GPFS drives re-balancing to improve performance
- ~~OS update on cluster, transfer, and login nodes~~

During maintenance whole cluster will not be usable, this includes:

- you will not be able to run jobs on cluster (SGE queuing system will be shutdown)
- med-login{1,2} nodes will not work reliably during this time
- med-transfer{1-2} nodes, and resources shared by them will be not available

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
        - [x] med-login1
        - [x] med-login2
        - [x] med-login3 (admin use only)
        - [x] med-transfer1
        - [x] med-transfer2
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
