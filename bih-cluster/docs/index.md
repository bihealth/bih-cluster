# Home

This is the documentation of the BIH high-performance compute (HPC), also called HPC 4 Research.
The BIH HPC cluster is maintained by CUBI (Core Unit Bioinformatics).

This documentation is maintained by BIH CUBI and the user community.
This is a living document that you can update and add to.
See [How-To: Contribute to this Document](how-to/misc/contribute) for details.

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
2. [Getting Help](help/hpc-talk) ([Writing Good Tickets](help/good-tickets); if no answer found, contact the [HPC Helpdesk](help/helpdesk)).
3. [For the Impatient](overview/for-the-impatient).
4. [The Cluster Tutorial: First Steps](first-steps/episode-0).

Then, continue reading through the manual.


!!! cite "Acknowledging BIH HPC Usage"
    Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    Please add your publications using the cluster to [this list](misc/publication-list).

## Maintenance Announcements

!!! attention "Current and Future Maintenances"

    - :warning: 2023, May 09-10: [BIH HPC / CUBI services downtime: maintenance at MDC data center](https://hpc-talk.cubi.bihealth.org/t/bih-hpc-cubi-services-downtime-maintenance-at-mdc-data-center-09-10-05-2023/525).
    - :nerd: August 30: Replace defect Nvidia V100 and make hpc-gpu-4 available to slurm users again.
    - :nerd: March 9: [Updated TensorFlow How-To to Slurm and Tensorflow 2.0](how-to/software/tensorflow/)
    - :warning: March 22-23: Compute and Storage Downtime.
    - :warning: March: Deprecation of using DRMAA.
    - :sparkles: March 1: New scheduler settings to address high job per user count.
    - :warning: February 14: Limiting allocateable memory per user.
    - :book: February 3: Adding Ganglia documentation.
    - :adhesive_bandage: February 3: Ganglia monitoring of GPFS and NVIDIA GPU metrics.
    - :warning: January 31: Enforcing usage of `localtmp` resource above 100MB of node-local storage.

See [Maintenance](admin/maintenance) for a detailed list of current, planned, and previous maintenance and update work.

## Connecting to the Cluster

You will need to perform some configuration steps after you have been registered with the cluster (via email to hpc-helpdesk@bih-charite.de by a group leader/PI).
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
