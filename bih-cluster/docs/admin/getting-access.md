# Getting Access

!!! tip "Get Access"

    A Charite or MDC account is required for accessing HPC 4 Research.

    1. **Register Users.** **Group leaders/PIs** register their members via hpc-gatekeeper@bihealth.de using the forms below.
    2. **Upload Key.** Upload your SSH key through the Charite and MDC infrastructure.
    3. **Connect.** `ssh -l <user>_<c or m>@hpc-login-1.cubi.bihealth.org` or using `hpc-login-2.cubi.bihealth.org`.

Access to the BIH HPC clusters is based on the **work groups** (also known as labs, units) and **projects** concepts. 
The work groups [data/folder structure](../../storage/storage-locations/) can only be accessed by the work group members and is not accessible to other cluster users.
Collaborative projects involving multiple PIs/groups should be realized using the **project** mechanism described below.
Please note, the **hot** near cluster fast storage is rather expensive and sparse resource which is not intended for long term storage. 

Independent group leaders at BIH/Charité/MDC can request a **work group** on the cluster and name **group members**. 
The work group **leader** (the group PI) has to take responsibility and undertake the necessary measures to ensure that all **group members** follow the [cluster policies and etiquette](../policies/) for fair usage and do not abuse the cluster resources for unethical or illegal purposes.
Major and/or continued violations may lead to exclusion of the work group form the cluster.  

The **leader** may also name one **delegate** (typically an IT savvy Post-Doc) who is allowed to take decision about cluster usage and work group on behalf of the leader. 
The above mentioned responsibilities stay with the work group leader.  

Data can also be managed in **projects** which allow 

1. cross-group collaboration spaces and
2. substructured access spaces, e.g., data where only selected group members may access data (e.g. not each internship student may get access to all valuable and potentially sensitive data).

Projects are also an excellent way to partition your data into set with different life 

HPC 4 Research is available for Charite, MDC, and BIH networks.

```bash
# Charite Users
host:~$ ssh -l user_c hpc-login-1.cubi.bihealth.org
# MDC Users
host:~$ ssh -l user_m hpc-login-1.cubi.bihealth.org
```

!!! important "Accounts and Email Adresses"

    - All users on the cluster must already have an account with either Charite/BIH or MDC.
    - Please only use email addresses from the institutions Charite, BIH, MDC.

## Work Groups

The process to create a new group is as follows.

1. The group leader sends an email to hpc-gatekeeper@bihealth.de and fills out the form below.
   Please consider the notes below on this page.
2. hpc-gatekeeper decides on the request and the corresponding objects are created on the cluster (users, groups, directories).
3. All new users are notified and further instructions are sent to them via email.

Subsequently, both owner and delegate (but **only owner and delegate**) can initiate changes (new users, resource changes etc.) to the group.

## Form

Example values are given in curly braces.

```
# Group "ag-doe"

Cluster: HPC 4 Research
Group leader/PI: {John Doe}
Delegate [optional]: {Max Mustermann}
Purpose of cluster usage [short]: {RNA-seq analysis in colorectal cancer}

Required resources:
- storage in TB: {1 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}
- Number of files [if >> 1M]: {less than 1M}

# Users

## User 1

- cluster: HPC 4 Research
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

We **strongly** dis-encourage on boarding non lab members into your group. 
This cause biases in usage accounting, may raise concerns in IT security and data privacy audits and also puts unfair responsibilities on the group leader. 
Please use the **project** mechanism described below. 

## Projects

Projects are very similar to work groups with the main distinction that
- users can be a member of multiple projects (no upper limit) 
- projects can be accessible for member of different groups 

Please note, project membership does not grant cluster access (a primary group membership is still required).

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

Cluster: HPC 4 Research
Project owner: {John Doe}, {doej_c}
Delegate [optional]: {Max Mustermann}, {musterm_c}
Purpose of cluster usage [short]: {RNA-seq data from dbGAP}

Required resources:
- storage in TB: {1 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}
- Number of files [if >> 1M]: {less than 1M}

Additional members:
- {Susi Sorglos}, {sorgls_c}
```

## Users

If you wish to add users to your AG, please use the following form. Note that users you want to add to a project need to be associated with an AG first.

The inquiry has to be send to hpc-gatekeeper@bih-charite.de, either by the PI or the delegate.

```
# User 1

- cluster: HPC 4 Research
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
- AG: ag-abcd

[etc]
```
