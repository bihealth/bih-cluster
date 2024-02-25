# Job Scheduler

Once logged into the cluster through the login nodes, users use the Slurm scheduler for job submission.
In Slurm nomenclature, cluster compute nodes are assigned to one or more **partitions**.
Submitted jobs are assigned to nodes according to the partition's configuration.

## Partitions

The BIH HPC has the partitions described below.
The cluster focuses on life science applications and not "classic HPC" with numerical computations using MPI.
Thus, all partitions except for `mpi` only allow to reserve resources on one node.
This makes the cluster easier to use as users don't have to explicitely specify this limit when submitting their jobs.

### `standard`

Jobs are submitted to the `standard` partition by default.
From the, the scheduler will route the jobs to their actual partition using the routing rule set described below.
You can override this routing by explicitely assigning a partition (but this is discouraged).

1. Jobs requesting a GPU resources are routed to the `gpu` queue.
2. Else, jobs requesting more than 200 GB of RAM are routed to the `highmem` queue.
3. Else, jobs are assigned to the partitions `debug`, `short`, `medium`, and `long` long depending on their configured maximal running time.
   The partitions are evaluated in the order given above and the first fitting partition will be used.

### `debug`

This partition is for very short jobs that should be executed quickly, e.g., for tests.
The job running time is limited to one hour and at most 128 cores can be used per user but the jobs are submitted with highest priority.

* **maximum run time:** 1 hour
* **maximum cores:** 128 cores per user
* **partition name:** `debug`
* **argument string:** maximum run time: `--time 01:00:00`

### `short`

This partition is for jobs running only few hours.
The priority of short jobs is high and many cores can be used at once to reward users for splitting their jobs into smaller parts.

* **maximum run time:** 4 hours
* **maximum cores:** 2000 cores
* **partition name:** `short`
* **argument string:** maximum run time: `--time 04:00:00`

### `medium`

This partition is for jobs running for multiple days.
Users can only allocate the equivalent of 4 nodes.

* **maximum run time:** 7 days
* **maximum cores:** 128 cores/slots (4 nodes)
* **partition name:** `medium`
* **argument string:** maximum run time: `--time 7-00:00:00`

### `long`

This partition is for long-running tasks.
Only one node can be reserved for so long to discourage really long-running jobs and encourage users for splitting their jobs into smaller parts.

* **maximum run time:** 14 days
* **maximum cores:** 32 cores/slots (1 node)
* **partition name:** `long`
* **argument string:** maximum run time: `--time 14-00:00:00`

### `gpu`

Jobs requesting GPU resources are automatically assigned to the `gpu` partition.

The GPU nodes are only part of the `gpu` partition so they are not blocked by normal compute jobs.
The maximum run time is relatively high (14 days) to allow for longer training jobs.
Contact [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) if you have longer running jobs that you really cannot make run any shorter for assistance.

For access to it you have register [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#gpu-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `gpu`
* **argument string:** select `$count` GPUs: `-p gpu --gres=gpu:$card:$count` (`card=tesla` or `card=a40`), maximum run time: `--time 14-00:00:00`

### `highmem`

Jobs requesting more than 200 GB of RAM are automatically routed to the `highmem` partition.

The high memory nodes are only part of the `highmem` partition so they are not blocked by normal compute jobs.
The maximum run time is relatively high (14 days) to allow for longer jobs.
Contact [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) for assistance if you have longer running jobs that you really cannot make run any shorter.

For access you must first contact [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#high-memory-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `highmem`
* **argument string:** `-p highmem`, maximum run time: `--time 14-00:00:00`

### `mpi`

Jobs are not routed automatically to the `mpi` partition but you have to explitely request the partition.
This is the only partition in which more than one node can be allocated to a job.

You can submit multi-node jobs into the `mpi` partition.
The maximum run time is relatively high (14 days) to allow for longer jobs.
Don't abuse this.
Contact [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) for assistance if you have longer running jobs that you really cannot make run any shorter.

For access you must first contact [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#gpu-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `highmem`
* **argument string:** `-p mpi`, maximum run time: `--time 14-00:00:00`

### `critical`

Jobs are not routed into `critial` automatically, and the partition has to be selected manually.

This partition is for time-critical jobs with deadlines.
For access to it you have to first ask [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de).
See [Resource Registration: Critical Partition](../admin/resource-registration.md#critical-partition) for details.

As long as the cluster is not very busy, requests for critical jobs will be granted most of the time.
However, do not use this partition without arranging with hpc-helpdesk as killing jobs will be used as the *ultima ratio* in case of such policy violations.

* **maximum run time:** 7 days
* **maximum cores:** 2000 cores/slots (48 nodes)
* **partition name:** `critical`
* **argument string:** maximum run time: `--time 7-00:00:00`
