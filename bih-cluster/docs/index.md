# Home

This is the documentation of the BIH high-performance compute (HPC) cluster.
This documentation is maintained bih BIH HPC IT, BIH CUBI (Core Unit Bioinformatics), and the user community.
*This is a living document, click the pencil icon next to the page header (e.g., "Home") to propose a change via GitHub.*

!!! tip "Get Access"
    1. **Register Users.** Group leaders register their members via  hpc-gatekeeper@bihealth.de.
    2. **Upload Key.** Upload your SSH key through the Charite and MDC Infrastructure.
    3. **Connect.** `ssh -l <user>_<c or m>@med-login1.bihealth.org`

!!! faq "Getting Help"
    Our helpdesk can be reached via email to hpc-helpdesk@bihealth.de.
    Please read our guide on [how to write good tickets](help/good-tickets) first.

!!! cite "Acknowledging BIH HPC Usage"
    Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    Please add your publications using the cluster to [this list](misc/publication-list).

!!! example "How to Edit these Pages"

    Click on the edit link at the top of each page as shown below.
    More details in [How-To: Contribute to Docs](how-to/contribute-to-docs).

    ![Click on the edit link of the page.](how-to/misc/figures/how-to-contribute.png)

## Maintenance Announcements

!!! attention "Current and Future Maintenances"
    - :test_tube: GPU upgrade - testing
    - :test_tube: SLURM scheduler - testing
    - :hourglass: Documentation migration - April 15, 2020
    - :hourglass: SGE scheduler end of life - June 31, 2020
    - :calendar: Head node migration -- ~May 31, 2020
    - :calendar: CentOS 8 migration -- summer 2020

See [Maintenance](admin/maintenance) for a detailed list of current, planned, and previous maintenance and update work.

## Connecting to the Cluster

You will need to perform some configuration steps after you have been registered with the cluster (via email to hpc-gatekeeper@bihealth.de by a group leader/PI).
Here are the most important points:

1. [Generating SSH Keys :key:](connecting/generate-key).
2. [Submitting the key :arrow_up: to Charite or to MDC](connecting/submit-key).
3. [Configuring your SSH client :wrench: on Linux and Mac or Windows](connecting/configure-ssh).
4. Bonus: [Connecting from external networks :flying_saucer:](connecting/from-external).

There are various other topics covered in the "Connecting" section that might be of interest to you.

## Getting Started

To get started, the following is a good reading order after your first successful connection to the cluster.

1. [Getting Help](help/helpdesk) (and [Writing Good Tickets](help/good-tickets)).
2. [Manual For the Impatient](manual/for-the-impatient).
3. [The Cluster Tutorial: First Steps](tutorial/first-steps-episode-0).

Then, continue reading through the manual.

## Documentation Structure

:arrow_left: The global table of contents is on the left, it is on the right for the current page :arrow_right:.

The documentation is structured as follows:

- **Administrative** provides information about administrative processes such as how to get access, register users, work groups, and projects.
- **Getting Help** explains how you can obtain for using the BIH HPC.
- **Manual** provides detailed information about the cluster setup.
  This includes the description of the hardware, network, software, and policies.
- **Connecting** provides technical help on performing connections to the cluster.
- **Slurm Scheduler** provides technical help on using the Slurm scheduler.
- **Best Practice** provides guidelines on recommended usage of certain aspects of the system.
- **Cubit Static Data** provides documentation about the static data (files) collection on the cluster.
- **How-To** provides short(ish) solutions for specific technical problems.
- **Miscellaneous** contains a growing list of pages not fitting anywhere else.
