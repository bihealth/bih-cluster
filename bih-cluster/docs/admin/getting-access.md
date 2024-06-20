# Getting Access
Access to the BIH HPC cluster is conceptually based on **user groups** (also known as labs or units) and **projects**.
Users have a relatively limited storage quota within their private home folder and store big data primarily within their group's work space or in project folders.
**Projects** are collaborative efforts involving multiple PIs/groups and are allocated separate storage space on the cluster.

Independent group leaders at BIH/Charité/MDC can request a **group** on the cluster and name group **members**. 
The work group **leader** (the group PI) bears the responsibility for the group's **members** and ensures that [cluster policies and etiquette](./policies.md) are followed.
In brief: Fair usage rules apply and the cluster ist not to be abused for unethical or illegal purposes.
Major and/or continued violations may lead to exclusion of the entire group.  

The group leader may also name one **delegate** (typically an IT-savvy Post-Doc) who is thereby allowed to take decision about cluster usage and work group management on behalf of the group leader. 
The above mentioned responsibilities stay with the group leader.  

!!! note

    - A Charité or MDC user account is required for accessing HPC 4 Research.
    - Please only use email addresses from the institutions Charite, BIH, or MDC in the forms below.

## Work Groups and Users
All cluster users are member of exactly one primary work group.
This affiliation is usually defined by real life organisational structures within Charité/BIH/MDC.
Leaders of independent research groups (PIs) can apply for a new cluster work group as follows:

1. The group leader sends an email to hpc-helpdesk@bih-charite.de and includes the filled-out form below.
   Please read the notes box before sending.
2. The HPC helpdesk decides on the request and creates corresponding objects on the cluster (users, groups, directories).
3. New users are notified and sent further instructions via email.

!!! warning "Important"
    Changes to an existing group (adding new users, changes in resources, etc.) can only be requested by group leaders and delegates.

### Form: New Group
Example values are given in curly braces.

```
# Group "ag-doe"
Group leader/PI: {John Doe}
Delegate [optional]: {Max Mustermann}
Purpose of cluster usage [short]: {RNA-seq analysis in colorectal cancer}

Required resources:
- Tier 1 storage: {1 TB}
- Tier 2 storage: {0 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}

# Users
## User 1
- first name: John
- last name: Doe
- affiliation: Charité, Department of Oncology
- institute email: john.doe@charite.de
- user has account with
    - [ ] BIH
    - [x] Charite
    - [ ] MDC
- BIH/Charité/MDC user name: doej

[etc.]
```

### Form: Add User to Group
Example values are given in curly braces.

```
# New user of AG Doe
- first name: Mia
- last name: Smith
- affiliation: Charité, Department of Oncology
- institute email: mia.smith@charite.de
- user has account with
    - [ ] BIH
    - [x] Charite
    - [ ] MDC
- BIH/Charité/MDC user name: smithm
```

!!! note "Notes"
    - All cluster groups must have an owner and may have one delegate.
    - Group ownership implies control but also accountability for their group files and members.
    - Users can only be members of one primary work group.
    - We **strongly** dis-encourage on-boarding non-lab members into your group.
      This cause biases in usage accounting, may raise concerns in IT security and data privacy audits and also puts unfair responsibilities on the group leader.

## Projects

Projects are secondary user groups to enable:

- collaboration and data sharing across different work groups,
- fine-grained allocation of additional storage resources,
- organising data in a fine-grained manner for better data lifecycle management. 

Project creation can be initiated by group leaders and group delegates as follows:

1. Send an email to hpc-helpdesk@bih-charite.de and includes the filled-out form below.
   Please read the notes box before sending.
2. The HPC helpdesk decides on the request and creates corresponding objects on the cluster (groups, directories).

!!! warning "Important"
    Changes to an existing project (adding new users, changes in resources, etc.) can only be requested by project owners and delegates.
    Please send us cluster user names for adding new project members.

### Form

Example values are given in curly braces.

```
# Project "doe-dbgap-rna"
Project owner: {John Doe}, {doej_c}
Delegate [optional]: {Max Mustermann}, {musterm_c}
Purpose of cluster usage [short]: {RNA-seq data from dbGAP}

Required resources:
- Tier 1 storage: {0 TB}
- Tier 2 storage: {1 TB}
- CPU hours per year: {approx. 1000}
- GPU hours per year: {none}

Additional members (cluster user names):
- {sorgls_c}
- ...
```

!!! note "Notes"
    - All projects must have one owner and may have one delegate.
    - Please note that we will enforce [kebab case]([url](https://en.wikipedia.org/wiki/Letter_case#Kebab_case)) for all project names and folders.
    - Tier 1 project storage will be supplemented with 10 TB of T1 scratch by default.
    - Users can be associated with multiple projects.
    - Project membership does not grant cluster access. A primary group affiliation is still required.
