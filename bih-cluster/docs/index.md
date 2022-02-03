# Home

This is the documentation of the BIH high-performance compute (HPC).
The BIH HPC cluster is maintained by CUBI (Core Unit Bionformatics).

This documentation is maintained by BIH CUBI, and the user community.

:arrow_left: The global table of contents is on the left, it is on the right for the current page :arrow_right:.

!!! tip "Connecting to the Cluster"

    - Web Access: https://hpc-portal.cubi.bihealth.org
    - SSH-Based Access:

        ```bash
        # Interactive Login as Charite/MDC user
        ssh -l user_c hpc-login-1.cubi.bihealth.org  # OR login-2...
        ssh -l user_m hpc-login-1.cubi.bihealth.org  # OR login-2...
        # File Transfer as Charite/MDC user
        sftp user_c@hpc-transfer-1.cubi.bihealth.org  # OR hpc-transfer-2...
        sftp user_m@hpc-transfer-1.cubi.bihealth.org  # OR hpc-transfer-2...
        # Of course, you can also log into the transfer nodes (no Slurm)
        ssh -l user_c hpc-transfer-1.cubi.bihealth.org  # OR hpc-transfer-2...
        ssh -l user_m hpc-transfer-1.cubi.bihealth.org  # OR hpc-transfer-2...
        ```

## Getting Started

To get started, the following is a suggested (in order) set of pages to read after your first successful connection to the cluster.

1. [Getting Access](admin/getting-access).
2. [Getting Help](help/helpdesk) (and [Writing Good Tickets](help/good-tickets)).
3. [For the Impatient](overview/for-the-impatient).
4. [The Cluster Tutorial: First Steps](first-steps/episode-0).

Then, continue reading through the manual.
This is a living document that you can update and add to.
See [How-To: Contribute to this Document](how-to/misc/contribute) for details.

!!! cite "Acknowledging BIH HPC Usage"
    Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    Please add your publications using the cluster to [this list](misc/publication-list).

## Maintenance Announcements

!!! attention "Current and Future Maintenances"
   
    - :adhesive_bandage: February 3: Ganglia monitoring of GPFS and NVIDIA GPU metrics.
    - :warning: January 31: Enforcing usage of `localtmp` resource above 100MB of node-local storage.
    - :sparkles: January 29: Reduced oversubscription, quotas on `/tmp` on the login nodes.
    - :sparkles: January 6: Oversubscribing physical cores to improve throughput.
    - :sparkles: December 27: Temporary File Handling Improvements
    - :sparkles: December 23-24: Cluster OS Upgrade
    - :sparkles: December 21-22: GPFS Upgrade

See [Maintenance](admin/maintenance) for a detailed list of current, planned, and previous maintenance and update work.

## Connecting to the Cluster

You will need to perform some configuration steps after you have been registered with the cluster (via email to hpc-gatekeeper@bihealth.de by a group leader/PI).
Here are the most important points:

1. [Generating SSH Keys :key: in Linux](connecting/generate-key/linux) or [Windows](connecting/generate-key/windows).
2. [Submitting the key :arrow_up: to Charite](connecting/submit-key/charite) or [to MDC](connecting/submit-key/mdc).
3. [Configuring your SSH client :wrench: on Linux and Mac](connecting/configure-ssh/linux) or [Windows](connecting/configure-ssh/windows).
4. Bonus: [Connecting from external networks :flying_saucer:](connecting/from-external).

There are various other topics covered in the "Connecting" section that might be of interest to you.

## Documentation Structure

The documentation is structured as follows:

- **Administrative** information about administrative processes such as how to get access, register users, work groups, and projects.
- **Getting Help** explains how you can obtain help in using the BIH HPC.
- **Overview** detailed information about the cluster setup.
  This includes the description of the hardware, network, software, and policies.
- **Connecting** technical help on connecting to the cluster.
- **First Steps** information for getting you started quickly.
- **Slurm Scheduler** technical help on using the Slurm scheduler.
- **Best Practice** guidelines on recommended usage of certain aspects of the system.
- **Static Data** documentation about the static data (files) collection on the cluster.
- **How-To** short(ish) solutions for specific technical problems.
- **Miscellaneous** contains a growing list of pages that don't fit anywhere else.
