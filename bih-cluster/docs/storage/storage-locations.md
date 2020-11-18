# Storage and Volumes: Locations

!!! info "HPC 4 Research Only"

    The following information is only valid for HPC 4 Research.

On the BIH HPC cluster, there are three kinds of entities: users, groups (*Arbeitsgruppen*), and projects.
Each user, group, and project has a central folder for their files to be stored.

## For the Impatient

### Storage Locations

Each user, group, and project directory consists of three locations (using `/fast/users/muster_c` as an example here):

- `/fast/users/muster_c/work`:
  Here, you put your large data that you need to keep.
  Note that there is no backup or snapshots going on.
- `/fast/users/muster_c/scratch`:
  Here, you put your large temporary files that you will delete after a short time anyway.
  **Data placed here will be automatically removed four weeks after creation.**
- `/fast/users/muster_c` (and all other sub directories):
  Here you put your programs and scripts and very important small data.
  By default, you will have a soft quota of 1GB (hard quota of 1.5GB, 7 days grace period).
  However, we create snapshots of this data (every 24 hours) and this data goes to a backup.

You can check your current usage using the command `bih-gpfs-report-quota user $USER`

### Do's and Don'ts

First and foremost:

- **DO NOT place any valuable data in `scratch` as it will be removed within 4 weeks.**

Further:

- **DO** set your `TMPDIR` environment variable to `/fast/users/$USER/scratch/tmp`.
- **DO** add `mkdir -p /fast/users/$USER/scratch/tmp` to your `~/.bashrc` and job script files.
- **DO** try to prefer creating fewer large files over many small files.
- **DO NOT** create multiple copies of large data.
  For sequencing data, in most cases you should not need more than raw times the size of the raw data (raw data + alignments + derived results).

## Introduction

This document describes the third iteration of the file system structure on the BIH HPC cluster.
This iteration was made necessary by problems with second iteration which worked well for about two years but is now reaching its limits.

## Organizational Entities

There are the following three entities on the cluster:

1. normal user accounts ("natural people")
2. groups *(Arbeitsgruppen)* with on leader and an optional delegate
3. projects with one owner and an optional delegate.

Their purpose is described in the document "User and Group Management".

## Storage/Data Tiers

The files fall into one of three categories:

1. **Home** data are programs and scripts of which there is relatively few but which is long-lived and very important.
   Loss of home data requires to redo manual work (like programming).

2. **Work** data is data of potential large size and has a medium life time and important.
   Examples are raw sequencing data and intermediate results that are to be kept (e.g., a final, sorted and indexed BAM file).
   Work data can time-consuming actions to be restored, such as downloading large amounts of data or time-consuming computation.

3. **Scratch** data is data that is temporary by nature and has a short life-time only.
   Examples are temporary files (e.g., unsorted BAM files).
   Scratch data is created to be removed eventually.

## Snapshots, Backups, Archive

- **A snapshot** stores the state of a data volume at a given time.
  File systems like GPFS implement this in a copy-on-write manner, meaning that for a snapshot and the subsequent "live" state, only the differences in data need to be store.d
  Note that there is additional overhead in the meta data storage.

- **A backup** is a copy of a data set on another physical location, i.e., all data from a given date copied to another server.
  Backups are made regularly and only a small number of previous ones is usually kept.

- **An archive** is a single copy of a single state of a data set to be kept for a long time.
  Classically, archives are made by copying data to magnetic tape for long-term storage.

## Storage Locations

This section describes the different storage locations and gives an overview of their properties.

### Home Directories

- **Location** `/fast/{users,groups,projects}/<name>` (except for `work` and `scratch` sub directories)
- the user, group, or project home directory
- meant for documents, scripts, and programs
- default quota for data: default soft quota of 1 GB, hard quota of 1.5 GB, grace period of 7 days
- quota can be increased on request with short reason statement
- default quota for metadata: 10k files soft, 12k files hard
- snapshots are regularly created, see Section \ref{snapshot-details}
- nightly incremental backups are created, the last 5 are kept
- *Long-term strategy:*
    users are expected to manage data life time independently and use best practice for source code and document management best practice (e.g., use Git).
    When users/groups leave the organization or projects ends, they are expected to handle data storage and cleanup on their own.
    Responsibility to enforce this is with the leader of a user's group, the group leader, or the project owner, respectively.

### Work Directories

- **Location** `/fast/{users,groups,projects}/<name>/work`
- the user, group, or project work directory
- meant for larger data that is to be used for a longer time, e.g., raw data, final sorted BAM file
- default quota for data: default soft quota of 1 TB, hard quota of 1.1 TB, grace period of 7 days
- quota can be increased on request with short reason statement
- default quota for metadata: 2 Mfile soft, 2.2M files hard
- no snapshots, no backup
- *Long-term strategy:*
    When users/groups leave the organization or projects ends, they are expected to cleanup unneeded data on their own.
    HPC IT can provide archival services on request.
    Responsibility to enforce this is with the leader of a user's group, the group leader, or the project owner, respectively.

### Scratch Directories

- **Location** `/fast/{users,groups,projects}/<name>/scratch`
- the user, group, or project scratch directory
- **files will be removed 2 weeks after their creation**
- meant for temporary, potentially large data, e.g., intermediate unsorted or unmasked BAM files, data downloaded from the internet for trying out etc.
- default quota for data: default soft quota of 200TB, hard quota of 220TB, grace period of 7 days
- quota can be increased on request with short reason statement
- default quota for metadata: 2M files soft, 2.2M files hard
- no snapshots, no backup
- *Long-term strategy:*
    as data on this volume is not to be kept for longer than 2 weeks, the long term strategy is to delete all files.

## Snapshot Details

Snapshots are made every 24 hours.
Of these snapshots, the last 7 are kept, then one for each day.

## Backup Details

Backups of the snapshots is made nightly.
The backups of the last 7 days are kept.

## Archive Details

BIH HPC IT has some space allocated on the MDC IT tape archive.
User data can be put under archive after agreeing with head of HPC IT.
The process is as describe in Section \ref{sop-data-archival}.

## Technical Implementation

As a quick (very) technical note:

There exists a file system `fast`.
This file system has three independent file sets `home`, `work`, `scratch`.
On each of these file sets, there is a dependent file set for each user, group, and project below directories `users`, `groups`, and `projects`.
`home` is also mounted as `/fast_new/home` and for each user, group, and project, the entry `work` links to the corresponding fileset in `work`, the same for scratch.
Automatic file removal from `scratch` is implemented using GPFS ILM.
Quotas are implemented on the file-set level.
