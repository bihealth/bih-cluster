# Home

This is the documentation of the BIH high-performance compute (HPC) cluster.
This documentation is maintained bih BIH HPC IT, BIH CUBI (Core Unit Bioinformatics), and the user community.

:arrow_left: The global table of contents is on the left, it is on the right for the current page :arrow_right:.

## Getting Started

To get started, the following is a good reading order after your first successful connection to the cluster.

1. [Getting Access](administrative/getting-access).
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
    - :test_tube: GPU upgrade - testing
    - :test_tube: SLURM scheduler - testing
    - :test_tube: Documentation migration - April 15, 2020
    - :test_tube: Head node migration -- June 15, 2020
    - :hourglass: SGE scheduler end of life - June 31, 2020
    - :calendar: CentOS 8 migration -- summer 2020
    - :calendar: `/fast` => `/data` -- autumn 2020

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

- **Administrative** provides information about administrative processes such as how to get access, register users, work groups, and projects.
- **Getting Help** explains how you can obtain for using the BIH HPC.
- **Overview** provides detailed information about the cluster setup.
  This includes the description of the hardware, network, software, and policies.
- **Connecting** provides technical help on performing connections to the cluster.
- **First Steps** provides information for getting you started quickly.
- **Slurm Scheduler** provides technical help on using the Slurm scheduler.
- **Best Practice** provides guidelines on recommended usage of certain aspects of the system.
- **Static Data** provides documentation about the static data (files) collection on the cluster.
- **How-To** provides short(ish) solutions for specific technical problems.
- **Miscellaneous** contains a growing list of pages not fitting anywhere else.
