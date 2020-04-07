# Nodes and Storage Volumes

!!! info "No mounting **on** the cluster itself."
    For various technical and security-related reasons, it is not possible to mount anything on the cluster nodes by users.
    That is, it is not possible to get file server mounts on the cluster nodes.
    For mounting the cluster storage on your computer, see [Connecting: SSHFS Mounts](../connecting/configure-ssh/linux.md#file-system-mount-via-sshfs).

This document gives an overview of the nodes and volumes on the cluster.

## Cluster Layout

![](figures/Cluster_Layout.png)

## Cluster Nodes

The following groups of nodes are available to cluster users.
There are a number of nodes that are invisible to non-admin staff, hosting the queue master and monitoring tools and providing backup storage for key critical data, but these are not shown here.

- `med-login1..3`
  - available as `med-login{1,2,3}.bihealth.org`
  - login nodes (IPs: `172.16.45.209`, `172.16.45.210`, `172.16.45.211`)
  - do not perform any computation on these nodes!
  - these nodes are not execution nodes
  - each process may at most use 1GB of RAM to increase stability of the node
- `med0101..0124,0127`
  - 25 standard nodes
  - Intel Xeon E5-2650 v2 @2.60Ghz, 16 cores x2 threading
  - 128 GB RAM
- `med0133..0164`
  - 32 standard nodes
  - Intel Xeon E5-2667 v4 @3.20GHz, 16 cores x 2 threading
  - 192 GB RAM
- `med0201..0264`
  - 64 nodes with Infiniband interconnect
  - Intel Xeon E5-2650 v2 @2.60Ghz, 16 cores x2 threading
  - 128 GB RAM
- `med0401..0405` special purpose/high-memory machines
  - Intel Xeon E5-4650 v2 @2.40GHz, 40 cores x2 threading
  - `med0401` and `med0402`
     - 1 TB RAM
  - `med0403` and `med0404`
     - 500 GB RAM
  - `med0405`
     - 2x "Tesla K20Xm" GPU accelleration cards (cluster resource `gpu`)
     - access limited to explicit GPU users
- `med0601..0616`
   - 16 nodes owned by CUBI
   - Intel Xeon E5-2640 v3 @2.60Ghz
   - 192 GB RAM
- `med0618..0633`
   - 16 nodes owned by CUBI
   - Intel Xeon E5-2667 v4 @3.20GHz, 16 cores x 2 threading
   - 192 GB RAM
- `med0701..0764`
   - 64 standard nodes
   - Intel Xeon E5-2667 v4 @3.20GHz, 16 cores x 2 threading
   - 192 GB RAM

If not noted anyway, currently no access restrictions apply per se.

## Cluster Volumes and Locations

The cluster has 2.1 PB of fast storage, currently available at `/fast`.
The storage runs on a DDN appliance using the IBM GPFS file system and is designed for massively parallel access from an HPC system.
In contrast to "single server" NFS systems, the system can provide large bandwidth to all cluster nodes in parallel as long as large data means relatively "few" files are read and written.
The system is not designed for access to many small files. 

The `/fast` volume has the following directories:

- `/fast/users/$USER` --
  For a user's data, the default quota is set to 10GB (see below for the `scratch` sub directory with no quota).
- `/fast/groups/$AG` --
  For the data of a working group lab, there is no default quota.
- `/fast/projects/$PROJECT`
  For sharing data between users that are not in the same lab, there is no default quota.
  :exclamation: Projects can be used for giving subsets of a whole lab access to a certain data set to implement the data privacy concept "Datensparsamkeit".

Inside each user, AG, and project directory, there is a `scratch` directory on which no quotas apply.
This is implemented by making each user, AG, and project directory a distinct "GPFS file set".

:bulb: We recommend that you store temporary files and intermediary/ephemeral results in the `scratch` directory and important files elsewhere.

Further, the cluster has 500 TB of slower memory mounted at `/slow` that is currently only used and available for mirroring Omics data and protecting them from accidental loss.

