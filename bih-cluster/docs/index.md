# Home

This is the documentation of the BIH high-performance computing (HPC) cluster, also called HPC 4 Research.
The BIH HPC cluster is maintained by CUBI (Core Unit Bioinformatics).

This documentation is maintained by BIH CUBI and the user community.
It is a living document that you can update and add to.
See [How-To: Contribute to this Document](how-to/misc/contribute) for details.

:arrow_left: The global table of contents is on the left, the one of the current page is on the right. :arrow_right:

!!! tip "Additional resources"

    - [User discussion forum](https://hpc-talk.cubi.bihealth.org/)
    - [Performance and workload monitoring](https://metrics.cubi.bihealth.org/public-dashboards/dc3e4d5b1ea049429abf39e412c47302)


## Getting Started

To get started, the following is a suggested (in order) set of pages to read after your first successful connection to the cluster.

1. [Getting Access](admin/getting-access).
2. [Getting Help](help/hpc-talk) ([Writing Good Tickets](help/good-tickets); if no answer found, contact the [HPC Helpdesk](help/helpdesk)).
3. [For the Impatient](overview/for-the-impatient).
4. [The Cluster Tutorial: First Steps](first-steps/episode-0).

Then, continue reading through the manual.


!!! note "Acknowledging BIH HPC Usage"
    Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    Please add your publications using the cluster to [this list](misc/publication-list).

## Maintenance Announcements
- :locomotive: Until mid 2024: Migration of all user data from GPFS to CephFS storage.
- :headstone: Late 2024: Retirement of GPFS/DDN storage.

See [Maintenance](admin/maintenance) for a detailed list of current, planned, and previous maintenance and update work.

## Connecting to the Cluster

You will need to perform some configuration steps after you have been registered with the cluster.
Here are the most important points:

1. [Generating SSH Keys :key: in Linux](connecting/generate-key/linux) or [Windows](connecting/generate-key/windows).
2. [Submitting the key :arrow_up: to Charite](connecting/submit-key/charite) or [to MDC](connecting/submit-key/mdc).
3. [Configuring your SSH client :wrench: on Linux and Mac](connecting/configure-ssh/linux) or [Windows](connecting/configure-ssh/windows).
4. Bonus: [Connecting from external networks :flying_saucer:](connecting/from-external).

There are various other topics covered in the "Connecting" section that might be of interest to you.

!!! tip "tl;dr"

    - Web Access: https://hpc-portal.cubi.bihealth.org
    - SSH-Based Access:

        ```bash
        # Interactive login (choose one)
        ssh username@hpc-login-1.cubi.bihealth.org
        ssh username@hpc-login-2.cubi.bihealth.org

        # File Transfer (choose one)
        sftp username@hpc-transfer-1.cubi.bihealth.org
        sftp username@hpc-transfer-2.cubi.bihealth.org

        # Interactive login into the transfer nodes (choose one)
        ssh username@hpc-transfer-1.cubi.bihealth.org
        ssh username@hpc-transfer-2.cubi.bihealth.org
        ```

## Documentation Structure

The documentation is structured as follows:

- **Administrative** information about administrative processes such as how to get access, register users, work groups, and projects.
- **Getting Help** explains how you can obtain help in using the BIH HPC.
- **Technical details** has detailed information about the cluster setup.
  This includes the description of the hardware, network, software, and policies.
- **Connecting** technical help on connecting to the cluster.
- **First Steps** information for getting you started quickly.
- **Storage** describes how and where files are stored.
- **Cluster Scheduler** technical help on using the Slurm scheduler.
- **Best Practice** guidelines on recommended usage of certain aspects of the system.
- **Static Data** documentation about the static data (files) collection on the cluster.
- **How-To** short(ish) solutions for specific technical problems.
- **Miscellaneous** contains a growing list of pages that don't fit anywhere else.
