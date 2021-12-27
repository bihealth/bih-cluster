
This page describes strictly enforced policies valid on the BIH HPC clusters.

The aim of the HPC systems is to support the users in their scientific work and relies on their cooperation.
First and foremost, the administration team enforces state of the art IT security and reliability practices through their organizational and operational processes and actions. We kindly ask user to follow the Cluster Etiquette describe below to allow for fair use and flexible access to the shared resources. Beyond this, policies are introduced or enforced only when required to ensure non-restrictive access to the resources themselves. Major or recurrent breaches of policies may lead to exclusion from service. 

We will update this list of policies over time.
Larger changes will be announced through the mailing list.

## Cluster Etiquette

1. The clusters are soft-partitioned shared resources that are made available under a "fair use" policy as far as possible.
2. The general assumption that if a user interferes with the work of others (e.g., by blocking compute slots) then this happens accidentally.
    - Please do not do this.
    - If you see this happening try to contact the user yourself (use `getent passswd $USER` to find out the user's office contact details).
    - Send an email to hpc-helpdesk@bihealth.de if you need administrative intervention.
3. All users must be subscribed to the cluster mailing list (they are subscribed automatically when the account is created).
4. When leaving please send an email to hpc-helpdesk@bihealth.de such that we can shutdown your account in an organized fashion.
   We also need to arrange for cleaning up your data.
5. The cluster mailing list bih-cluster@bihealth.de is the primary contact channel for administration to user.
   Users must be subscribed to the mailing list.
   Users must follow the announcements, failure to do so can lead to missing important policy changes and thus losing access to the cluster or data.
6. **Do not perform any computation on the login nodes.**
   This includes: running `conda`, archive management tools such as `tar`, `(un)zip`, or `gzip`.
   You should probably only run `screen`/`tmux` and maybe a text editor there.

## Cluster Policies

### File System Policies

In the case of violations marked with a shield (:shield:) administration reserves the right to remove write and possibly read permission to the given locations.
Policies marked with a robot (:robot:) are automatically enforced.

1. Storage on the GPFS file system is a sparse resource try to use both data volume and file sparingly.
   Note well that small files above ~4KB take up at least 8MB of space.
2. Default quotas are as follows (each user, group, project has a `home`, `work`, and `scratch` volume).
   You can request an increase by an email to hpc-gatekeeper@bihealth.de.
    - `home` 10k files, 1GB space
    - `work` 2M files, 1TB space
    - `scratch` 20M files, 200TB space
3. The overall throughput limit is 10GB/sec.
   Try not to overload the cluster I/O wise.
4. :shield: :robot: User home/work/group file sets have to be owned by the user, group is `hpc-users` and mode is `u=rwx,go=`; POSIX ACLs are prohibited.
    This policy is automatically enforced every 5 minutes.
5. :shield: :robot: Group and project home/work/group file sets have to be owned by the owner, group set to the corresponding unix group and mode is `u=rwx,g=rwxs,o=`; POSIX ACLs are prohibited.
    This policy is automatically enforced every 5 minutes.
6. :robot: All files in scratch will be moved into a read-only "trash can" inside `scratch/BIH_TRASH` after 14 days (by `mtime`) over night.
   Trash directories will be removed after 14 further days.
    - Users can arrange with hpc-gatekeeper@bihealth.de to keep files longer by using `touch` on files in `scratch` and subsequently bumping the `mtime`.
    - :warning: In the case of abuse of this mechanism / failure to communicate with hpc-gatekeeper, administration reserves the right to drastically reduce scratch quota of affected users and employ other measures to ensure stability of operations.
    - You can learn more in the [Automated Scratch Cleanup](../storage/scratch-cleanup.md) section.
7. :shield: Administration will not delete any files (outside of `/tmp`).
   In the case that users need to delete files that they can access but not update/delete, administration will either give write permissions to the Unix group of the work group or project or change the owner to the owner/delegate of this group.
   This can occur in a group/project directory of a user who has left the organization.
   In the case that a user laeaves the organiation, the owner/delegate of the hosting group can request getting access to the user's files with the express agreement of this user.
8. Only use `/tmp` in Slurm-controlled jobs.
   This will enforce that Slurm can clean up after you.

### Connections

Network connections are a topic important in security.
In the case of violations marked with a shield (:shield:) administration reserves the right to terminate connections without notice and perform other actions.

1. Data transfers should happen through the transfer nodes (HPC 4 Research) and/or the compute nodes themselves.
2. :shield: The cluster is not meant as a "hop node".
   Do not use it to connect to the login node first and then jump to another host outside of the cluster network. Doing so is a breach of cluster policies and quite possibly your organization's IT security policies
3. :shield: As a corollary, SSH reverse tunnels are strictly prohibited.
4. Outgoing connections are meant for data transfers only (in other words: using SSH/SCP to download file is fine).
5. :shield: Do not leave outgoing connections open longer than necessary.
6. Sessions of `screen` and `tmux` are only allowed to run on the head nodes.
   They will be terminated automatically on the compute nodes.

### Interactive Use

1. Interactive sessions block resources to the scheduler.
   Reduce interactive use to the minimal time and resources possible.
2. The cluster is optimized for batch processing.
   Interactive use is a secondary aim.
   Administration attempts to strike a good balance here but batch usage is most important.
   Consider using our Open on Demand service for interactive use. 
3. Interactive use should happen through the Slurm scheduler (`srun`).
4. SSH connections to the nodes are allowed for monitoring purposes but not meant for computation.
   Administration enforces this by restricting all jobs outside of Slurm to use at most 1 core and 128 MB of RAM.
   This limit is enforced per node per user with Linux cgroups.

### GPU Use

1. Interactive sessions block resources to the scheduler.
   Interactive GPU use is discouraged.
2. Accessing GPUs outside of the Slurm scheduler has been disabled by administration.

### Account Policies

1. :shield: Sharing accounts and/or credentials is strictly prohibited.
   Doing so is a breach of cluster policies and certainly also of your organization's IT security policies.
2. :shield: Hosting shared services on the cluster is also strictly prohibited.
        - This includes Jupyter servers that shall only be used by the user starting them, this also includes work schedulers such as Dask.
        - You can assume that the cluster internal network is secure and you do not have to encrypt connections between nodes.
        - Connections towards outside of the clsuter must be encrypted (e.g., via SSH tunnels; incoming ones as reverse tunneling is prohibited, see above).
        - Access to any service must be protected by appropriate means, e.g., passwords, tokens or client certifictes.

### Maintenances

1. Maintenance that are expected to cause major service interruptions (the whole system becomes unusable and/or jobs might be prevented to run etc.) are announced 14 days in advance.
2. Maintenance of login nodes (e.g., reboot one node while the other is still available) are announced 7 days in advance.
3. Maintenance of transfer nodes are announced 1 day in advance.
   Rationale: transfer nodes expected to not have any interactive sessions running.

### Credentials Policies

1. Login on HPC 4 Research is currently based on SSH keys only.
   HPC 4 Clinic currently allows login with password or SSH keys.
2. SSH keys must deposited with the host organizations (Charite/MDC) as documented.
3. For technical reasons, the compute nodes also use the `~/.ssh/authorized_keys` file but their usage is discouraged.
