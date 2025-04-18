# Home
Welcome to the user documentation of the BIH high-performance computing (HPC) cluster, also called HPC 4 Research.
The BIH HPC cluster is managed by [CUBI](https://cubi.bihealth.org) (Core Unit Bioinformatics).
This documentation is maintained by BIH CUBI and the user community.
It is a living document that you can update and add to.
See [How-To: Contribute to this Document](how-to/misc/contribute.md) for details.

:arrow_left: The global table of contents is on the left, the one of the current page is on the right. :arrow_right:

!!! tip "Additional resources"

    - [User discussion forum](https://hpc-talk.cubi.bihealth.org/)
    - [Performance and workload monitoring](https://metrics.cubi.bihealth.org/public-dashboards/dc3e4d5b1ea049429abf39e412c47302)


## Getting Started
Read the following set of pages (in order) to learn how to get access and connect to the cluster.

1. [Getting Access](admin/getting-access.md)
2. [Connecting](connecting/overview.md)
3. [Storage](storage/storage-locations.md)
5. [Slurm](slurm/overview.md)
6. [Getting Help](help/hpc-talk.md) ([Writing Good Tickets](help/good-tickets.md); if no answer found, contact the [HPC Helpdesk](help/helpdesk.md)).
7. [HPC Tutorial](hpc-tutorial/episode-0.md)

!!! note "Acknowledging BIH HPC Usage"
    Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    Please add your publications using the cluster to [this list](misc/publication-list.md).

## News & Maintenance Announcements
- :woman_technologist: February 10th 2025: New GPU nodes `hpc-gpu-[9-11]` with 4× NVIDIA L40 each.
- :snowflake: January 2025: Rolling upgrade to Rocky Linux 9 of all cluster nodes & VMs.
- :jack_o_lantern: October 22nd 2024: Kernel & SLURM upgrades on (almost) all nodes.
- :ram: July 16th: New high-memory node `hpc-mem-5` with 4 TB of RAM.
- :locomotive: Until autumn 2024: Operation Exodus – Migration of all data from GPFS to CephFS storage.
- :maple_leaf: September 30th 2024: Unmounting of `/fast` on all non-transfer nodes.
- :headstone: October 31st 2024: Retirement of GPFS/DDN storage.

See [Maintenance](admin/maintenance.md) for a detailed list of current, planned, and previous maintenance and update work.

## Technical Details
If you are interested in how this HPC cluster is set up on a technical level, we got you covered.
There is [an entire section](./overview/for-the-impatient.md) on this.

## Documentation Structure
The documentation is structured as follows:

- **Administrative** information about administrative processes such as how to get access, register users, work groups, and projects.
- **Connecting** technical help for connecting to the cluster.
- **Storage** describes how and where files are stored.
- **HPC tutorial** a first demo project for getting you started quickly.
- **Cluster Scheduler** technical help for using the Slurm scheduler.
- **OnDemand Portal** introduces web HPC access.
- **Best Practice** guidelines on recommended usage of certain aspects of the system.
- **Static Data (Cubit)** documentation about the static data (files) collection on the cluster.
- **How-To** short(ish) solutions for specific technical problems.
- **Getting Help** explains how you can obtain help in using the BIH HPC.
- **Miscellaneous** contains a growing list of pages that don't fit anywhere else.
