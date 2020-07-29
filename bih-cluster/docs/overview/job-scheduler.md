# Job Scheduler

Once logged into the cluster through the login nodes, users use the Slurm scheduler for job submission.
In Slurm nomenclature, cluster compute nodes are assigned to one or more **partitions**. Submitted jobs are assigned to nodes according to the partition's configuration.

!!! note "Selecting Partitions"
    Users should only select partitions explicitly if they need special resources such as GPUs, high-memory nodes, Infiniband, or the critical partition.
    Otherwise, they should let Slurm decide on the optimal assignment.

!!! todo "MPI Partition"
    We recently introduced an `mpi` partition that allows you to allocate more than one node.
    We will update this soon.

## Partitions

The BIH HPC has the following partitions.

!!! important "Use `mpi` Partition for Multi-Node Jobs"
    All partitions except for `mpi` have their maximum number of available nodes set to `1`.
    This makes the system easier to use for single-node/multi-core programs as users don't have to specify `--min-nodes=1 --max-nodes=1`.
    Simply use the `mpi` partition for multi-node jobs.

### `debug`

!!! todo "debug/default partition"
    Actually, the default partition is the "debug" partition for now.
    We will rename things soon.

### `default`

This partition is for normal jobs running for less than 4 hours.
The number of overall jobs allowed per user is high to reward users for splitting their computational needs into smaller parts (this makes them easier for the scheduler to deal with).

* **maximum run time:** 4 hours
* **maximum cores:** >6000 cores (all nodes, **no limit**)
* **partition name:** `default`
* **priority:** **`default`** < `medium` < `long` < `critical` < `mini`
* **argument string:** maximum run time: `--time 04:00:00`

### `medium`

This partition is for jobs running for multiple days.

* **maximum run time:** 7 days
* **maximum cores:** 128 cores/slots (4 nodes)
* **partition name:** `medium`
* **priority:** `default` < **`medium`** < `long` < `critical` < `mini`
* **argument string:** maximum run time: `--time 7-00:00:00`

### `long`

This partition is for long-running tasks.

* **maximum run time:** 28 days
* **maximum cores:** 32 cores/slots (1 node)
* **partition name:** `long`
* **priority:** `default` < `medium` < **`long`** < `critical` < `mini`
* **argument string:** maximum run time: `--time 28-00:00:00`

### `critical`

This partition is for time-critical jobs with deadlines.
For access to it you have to first ask [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de).
See [Resource Registration: Critical Partition](../admin/resource-registration.md#critical-partition) for details.

As long as the cluster is not very busy, requests for critical jobs will be granted most of the time.
However, do not use this partition without arranging with hpc-gatekeeper as killing jobs will be used as the *ultima ratio* in case of such policy violations.

* **maximum run time:** 7 days
* **maximum cores:** 1536 cores/slots (48 nodes)
* **partition name:** `critical`
* **priority:** `default` < `medium` < `long` < **`critical`** < `mini`
* **argument string:** maximum run time: `--time 7-00:00:00`

### `gpu`

The GPU nodes are only part of the `gpu` partition so they are not blocked by normal compute jobs.
The maximum run time is relatively high (14 days) to allow for longer training jobs.
Contact [hpc-helpdesk@bihealth.de](mailto:hpc-helpdesk@bihealth.de) if you have longer running jobs that you really cannot make run any shorter for assistance.

For access to it you have register [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#gpu-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `gpu`
* **argument string:** select `$count` nodes: `-p gpu --gres=gpu:tesla:$count`, maximum run time: `--time 14-00:00:00`

### `highmem`

The high memory nodes are only part of the `highmem` partition so they are not blocked by normal compute jobs.
The maximum run time is relatively high (14 days) to allow for longer jobs.
Contact [hpc-helpdesk@bihealth.de](mailto:hpc-helpdesk@bihealth.de) for assistance if you have longer running jobs that you really cannot make run any shorter.

For access you must first contact [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#high-memory-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `highmem`
* **argument string:** `-p highmem`, maximum run time: `--time 14-00:00:00`

### `mpi`

You can submit multi-node jobs into the `mpi` partition.
The maximum run time is relatively high (14 days) to allow for longer jobs.
Don't abuse this.
Contact [hpc-helpdesk@bihealth.de](mailto:hpc-helpdesk@bihealth.de) for assistance if you have longer running jobs that you really cannot make run any shorter.

For access you must first contact [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#gpu-nodes) for details.

* **maximum run time:** 14 days
* **partition name:** `highmem`
* **argument string:** `-p mpi`, maximum run time: `--time 14-00:00:00`
