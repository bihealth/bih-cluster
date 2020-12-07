
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
3. The overall throughput limit is 20GB/sec.
   Try not to overload the cluster I/O wise.
4. :shield: :robot: User home/work/group file sets have to be owned by the user, group is `hpc-users` and mode is `u=rwx,go=`; POSIX ACLs are prohibited.
    This policy is automatically enforced every 5 minutes.
5. :shield: :robot: Group and project home/work/group file sets have to be owned by the owner, group set to the corresponding unix group and mode is `u=rwx,g=rwxs,o=`; POSIX ACLs are prohibited.
    This policy is automatically enforced every 5 minutes.

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
