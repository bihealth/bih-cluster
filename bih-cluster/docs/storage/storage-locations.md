# Storage and Volumes: Locations
This document describes the forth iteration of the file system structure on the BIH HPC cluster.
It was made necessary because the previous file system was no longer supported by the manufacturer and we since switched to distributed [Ceph](https://ceph.io/en/) storage.

## Organizational Entities
There are the following three entities on the cluster:

1. **Users** *(real people)*
2. **Groups** *(Arbeitsgruppen)* with one leader and an optional delegate
3. **Projects** with one owner and an optional delegate

Each user, group, and project can have storage folders in different locations.

## Data Types and Storage Tiers
Files stored on the HPC fall into one of three categories:

1. **Home** folders store programs, scripts, and user config i.&nbsp;e. long-lived and very important files. 
Loss of this data requires to redo manual work (like programming).

2. **Work** folders store data of potentially large size which has a medium life time and is important.
Examples are raw sequencing data and intermediate results that are to be kept (e.&nbsp;g. sorted and indexed BAM files).
Work data requires time-consuming actions to be restored, such as downloading large amounts of data or long-running computation.

3. **Scratch** folder store temporary files with a short life-time.
Examples are temporary files (e.&nbsp;g. unsorted BAM files).
Scratch data is created to be removed eventually.

Ceph storage comes in two types which differ in their I/O speed, total capacity, and cost.
They are called **Tier 1** and **Tier 2** and sometimes **hot storage** and **warm storage**.
In the HPC filesystem they are mounted in `/data/cephfs-1` and `/data/cephfs-2`.

- Tier 1 storage is fast, relatively small, expensive, and optimized for performance.
- Tier 2 storage is slow, big, cheap, and built for keeping large files for longer times.

Storage quotas are imposed in these locations to restrict the maximum size of folders.
Amount and utilization of quotas is communicated via the [HPC Access](https://hpc-access.cubi.bihealth.org/) web portal.

### Home Directories
Location: `/data/cephfs-1/home/`

Only users have home directories on Tier 1 storage.
This is the starting point when starting a new shell or SSH session.
Important config files are stored here as well as analysis scripts and small user files.
Home folders have a [strict storage quota](./home-quota.md) of 1 GB.

### Work Directories
Location: `/data/cephfs-1/work/`

Groups and projects have work directories on Tier 1 storage.
User home folders contain a symlink to their respective group's work folder.
Files shared within a group/project are stored here as long as they are in active use.
Work folders are generally limited to 1 TB per group.
Project work folders are allocated on an individual basis.

### Scratch Space
Location: `/data/cephfs-1/scratch/`

Groups and projects have scratch space on Tier 1 storage.
User home folders contain a symlink to their respective group's scratch space.
Meant for temporary, potentially large data e. g. intermediate unsorted or unmasked BAM files, data downloaded from the internet etc.
Scratch space is generally limited to 10 TB per group.
Projects are allocated scratch on an individual basis.
Files in scratch will be [automatically removed](scratch-cleanup.md) 2 weeks after their creation.

### Tier 2 Storage
Location: `/data/cephfs-2/`

This is where big files go when they are not in active use.
Groups are allocated 10 TB of Tier 2 storage by default.
File quotas here can be significantly larger as space is much cheaper and more abundant than on Tier 1.

!!! note
    Tier 2 storage is currently not accessible from HPC login nodes.

### Overview

| Tier | Function        | Path                                           | Default Quota |
|:-----|:----------------|:-----------------------------------------------|--------------:|
|    1 | User home       | `/data/cephfs-1/home/users/<user>`             | 1 GB          |
|    1 | Group work      | `/data/cephfs-1/work/groups/<group>`           | 1 TB          |
|    1 | Group scratch   | `/data/cephfs-1/scratch/groups/<group>`        | 10 TB         |
|    1 | Project work    | `/data/cephfs-1/work/projects/<project>`       | On request    |
|    1 | Project scratch | `/data/cephfs-1/scratch/projects/<project>`    | On request    |
|    2 | Group           | `/data/cephfs-2/unmirrored/groups/<group>`     | 10 TB         |
|    2 | Project         | `/data/cephfs-2/unmirrored/projects/<project>` | On request    |
|    2 | Group           | `/data/cephfs-2/mirrored/groups/<group>`       | On request    |
|    2 | Project         | `/data/cephfs-2/mirrored/projects/<project>`   | On request    |

## Snapshots and Mirroring
Snapshots are incremental copies of the state of the data at a particular point in time. 
They provide safety against various "Ops, did I just delete that?" scenarios, meaning they can be used to recover lost or damaged files.
Depending on the location and Tier, CephFS creates snapshots in different frequencies and retention plans.

| Location  | Path  | Retention policy  |
|-----------|-------|-------------------|
| User homes                | /data/cephfs-1/home/users/  | Daily for 14 d |
| Group/project work        | /data/cephfs-1/work/        | Daily for 5 d  |
| Group/project scratch     | /data/cephfs-1/scratch/     | Daily for 3 d  |
| Group/project mirrored    | /data/cephfs-2/mirrored/    | Daily for 7 d, weekly for 8 w  |
| Group/project unmirrored  | /data/cephfs-2/unmirrored/  | Daily for 7 d, weekly for 8 w  |
| SODAR                     | | Daily for 7 d |

Some parts of Tier 1 and Tier 2 snapshots are also mirrored into a separate fire compartment within the data center.
This provides an additional layer of security i. e. physical damage to the servers.

### Accessing Snapshots
To access snapshots simply navigate to the `.snap/` sub-folder of the respective location.
This special folder exists on all levels of the CephFS file hierarchy, so even in your user home directory.
Inside you will find one folder per snapshot created and in those a complete replica of the respective folder at the time of snapshot creation.

For example:

- `/data/cephfs-1/home/.snap/<some_snapshot>/users/<your_user>/` same as:
- `/data/cephfs-1/home/users/<your_user>/.snap/<some_snapshot>`
- `/data/cephfs-1/work/.snap/<some_snapshot>/groups/<your_group>/`
- `/data/cephfs-2/unmirrored/.snap/<some_snapshot>/projects/<your_project>/`

Here is a simple example of how to restore a file:

```sh
$ cd /data/cephfs-2/unmirrored/groups/cubi/.snap/scheduled-2024-03-11-00_00_00_UTC/
$ ls -l
important_file.txt
$ cp important_file.txt /data/cephfs-2/unmirrored/groups/cubi/
```

## Technical Implementation
### Tier 1
- Fast & sparse
- mounted on `/data/cephfs-1`
- Currently 12 Nodes with 10 × 14 TB NVME SSD each
    - 1.68 PB raw storage
    - 1.45 PB erasure coded (EC 8:2)
    - 1.23 PB usable (85 %, Ceph performance limit)
- For typical CUBI use case 3 to 5 times faster I/O then the old DDN
- Hardware costs:
    - One node/chunk: 50.000 € (150 TB)
    - ca. 350 €/TB

### Tier 2
- Slower but much more affordable
- Mounted on /data/cephfs-2, comes in two flavours
  - /data/cephfs-2/mirrored
  - /data/cephfs-2/unmirrored
- Currently 24 nodes in two clusters in two different fire compartments. Each node with ~40 HDDs each 16–18 TB capacity
- Current useable capacity 10+ PB (mirrored nad unmirror combined)
- Hardware costs:
  - ca. 50 € per TB unmirrored
  - ca. 100 € per TB mirrored
  - small chunk extension possible
