# Getting Access

=== "HPC 4 Research"

    !!! tip "Get Access"

        A Charite or MDC account is required for accessing HPC 4 Research.

        1. **Register Users.** Group leaders register their members via  hpc-gatekeeper@bihealth.de.
        2. **Upload Key.** Upload your SSH key through the Charite and MDC Infrastructure.
        3. **Connect.** `ssh -l <user>_<c or m>@login-1.research.bihealth.org`

=== "HPC 4 Clinic"

    !!! tip "Get Access"

        For using HPC 4 Clinic, a Charite account is required.

        1. **Register Users.** Group leaders register their members via  hpc-gatekeeper@bihealth.de.
        2. **Upload Key.** Upload your SSH key through the Charite and MDC Infrastructure.
        3. **Connect.** `ssh -l <user>_<c or m>@login-1.clinic.bihealth.org`

Access to the BIH HPC clusters is based on work groups (also known as labs, units).
Each group is headed by a leader (also known as principle investigator/PI).
Data can also be managed in **project** which allow cross-group collaboration but also providing a limited access space, e.g., for controlled data access data where only a few group members may access the data.


=== "HPC 4 Research"

    HPC 4 Research is available for Charite, MDC, and BIH networks.

    ```bash
    # Charite Users
    host:~$ ssh -l user_c login-1.research.hpc.bihealth.org
    # MDC Users
    host:~$ ssh -l user_m login-1.research.hpc.bihealth.org
    ```

    !!! important "Accounts and Email Adresses"

        - All users on the cluster must already have an account with either Charite/BIH or MDC.
        - Please only use email addresses from the institutions Charite, BIH, MDC.

=== "HPC 4 Clinic"

    HPC 4 Clinic is only from the Charite network.

    ```bash
    host:~$ ssh -l user login-1.clinic.hpc.bihealth.org
    ```

    !!! important "Accounts and Email Adresses"

        - All users on the cluster must already have an account either Charite/BIH.
        - You must use your Charite/BIH email address.


## Work Groups

The process to create a new group is as follows.

1. The group leader sends an email to hpc-gatekeeper@bihealth.de and fills out the form below.
   Please consider the notes below on this page.
2. hpc-gatekeeper decides on the request and the corresponding objects are created on the cluster (users, groups, directories).
3. All new users are notified and further instructions are sent to them via email.

Subsequently, both owner and delegate can initiate changes (new users, resource changes etc.) to the group.

## Form

Example values are given in curly braces.

```
# Group "ag-doe"

Cluster: {HPC 4 Clinic or HPC 4 Research}
Group leader/PI: {John Doe}
Delegate [optional]: {Max Mustermann}
Purpose of cluster usage [short]: {RNA-seq analysis in colorectal cancer}

Required resources:
- storage in TB: {1 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}
- Number of files [if >> 1M]: {less than 1M}

Users for each member:

# User 1

- cluster: {HPC 4 Clinic or HPC 4 Research}
- first name: John
- last name: Doe
- affiliation: Charite, Department of Oncology
- institute email: john.doe@charite.de
- institute phone: 030-8445-0
- user has account with
    - [ ] BIH
    - [ ] Charite
    - [ ] MDC
- BIH/Charite/MDC user name: doej
- duration of cluster access (max 1 year): 2020-03-30 to 2021-03-30

[etc.]
```

## Notes

- Work groups on the cluster must have an owner (its leader, principal investigator, etc.)
- Group ownership implies control but also accountability for their work group and members.
- Through the delegate mechanism, control can be delegated to up to one person (e.g., post-doc in the lab).
- **Users can be members of one work group only.**
  For multi-group collaborations, please use the project mechanism described below.

## Projects

Projects are very similar to work groups with the main distinction that

- users can be a member of more than one project, but
- project membership does not grant cluster access (group membership is still required).

Project creation can be initiated by group leaders and group delegates with the following process.

1. The initiator sends an email to hpc-gatekeeper@bihealth.de and fills out the following form.
   Please consider the notes below on this page.
2. hpc-gatekeeper decides on the request and the corresponding objects are created on the cluster (users, groups, directories).
3. All new users are notified and further instructions are sent to them via email.

Subsequently, both owner and delegate can initiate changes (new users, resource changes etc.) to the project.

## Form

Example values are given in curly braces.

```
# Project "doe-dbgap-rna"

Cluster: {HPC 4 Research}
Project owner: {John Doe}
Delegate [optional]: {Max Mustermann}
Purpose of cluster usage [short]: {RNA-seq data from dbGAP}

Required resources:
- storage in TB: {1 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}
- Number of files [if >> 1M]: {less than 1M}

Additional members:
- Susi Sorglos <Susi.Sorglos@charite.de>
```
