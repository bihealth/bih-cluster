# Job Scheduler

Once logged into the cluster through the login nodes, users use the Slurm scheduler for job submission.
In Slurm nomenclature, the nodes are assigned to one or more **partitions** and jobs are then assigned to nodes according to the partition's configuration.

!!! note "Selecting Partitions"
    Users should only select partitions explictely if they need special resources such as GPUs, high-memory nodes, Infiniband, or the critical partition.
    Otherwise, they should let Slurm decide about the optimal execution.

!!! todo "MPI Partition"
    We recently introduced an `mpi` partition that allows you to allocte more than one node.
    We will update this soon.

## Partitions

The BIH HPC has the following partitions.

!!! important "Use `mpi` Partition for Multi-Node Jobs"
    All partitions except for the `mpi` partition have the maximal number of availble nodes set to `1`.
    This makes the system easier to use for single-node/multi-core apps as users don't have to specify `--min-nodes=1 --max-nodes=1`.
    Simply use the `mpi` partition for multi-node jobs.

### `debug`

!!! todo "debug/default partition"
    Actually, the default partition is the "debug" partition for now.
    We will rename things soon.

### `default`

This partition is for normal jobs running in less than 4 hours.
The number of overall allows jobs by one user is high to reward users for chopping down their tasks into smaller parts and making them scheduler friendly.

* **maximal running time:** 4 hours
* **maximal cores:** >6000 cores (all nodes, **no limit**)
* **priority:** **`default`** < `medium` < `long` < `critical` < `mini`
* **argument string:** maximal running time: `--time 04:00:00`

### `medium`

This partition is for jobs running for multiple days.

* **maximal running time:** 7 days
* **maximal cores:** 128 cores/slots (4 nodes)
* **priority:** `default` < **`medium`** < `long` < `critical` < `mini`
* **argument string:** maximal running time: `-t 7-00:00:00`

### `long`

This partition is for long-running tasks.

* **maximal running time:** 28 days
* **maximal cores:** 32 cores/slots (1 node)
* **project name:** `long`
* **queue name:** `long.q`
* **priority:** `default` < `medium` < **`long`** < `critical` < `mini`
* **argument string:** maximal running time: `-t 28-00:00:00`

### `critical`

This partition is for time-critical jobs with deadlines.
For access to it you have to first ask [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de).
See [Resource Registration: Critical Partition](../admin/resource-registration.md#critical-partition) for details.

As long as the cluster is not very busy, requests for critical jobs will be granted most of the time.
However, do not use this queue without arranging with hpc-gatekeeper as killing jobs will be used as the *ultima ratio* in case of such policy violations.

* **maximal running time:** 7 days of maximal running time
* **maximal cores:** 1536 cores/slots (48 nodes)
* **project name:** `critical`
* **priority:** `default` < `medium` < `long` < **`critical`** < `mini`
* **argument string:** maximal running time: `-t 7-00:00:00`

### `gpu`

The GPU nodes are only part of the `gpu` partition so they are not blocked by normal compute jobs.
The maximal running time is relatively high (14 days) to allow for longer training jobs.
Contact hpc-helpdesk@bihealth.de if you have longer running jobs that you really cannot make run any shorter for assistance.

For access to it you have register [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#gpu-nodes) for details.

* **maximal running time:** 14 days of maximal running time
* **partition name:** `gpu`
* **argument string:** select `$count` nodes: `-p gpu --gres=gpu:tesla:$count`, maximal running time: `-t 14-00:00:00`

### `highmem`

The high membory nodes are only part of the `highmem` partition so they are not blocked by normal compute jobs.
The maximal running time is relatively high (14 days) to allow for longer jobs.
Contact hpc-helpdesk@bihealth.de if you have longer running jobs that you really cannot make run any shorter for assistance.

For access to it you have register [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#high-memory-nodes) for details.

* **maximal running time:** 14 days of maximal running time
* **partition name:** `highmem`
* **argument string:** `-p highmem`, maximal running time: `-t 14-00:00:00`

### `mpi`

You can submit multi-node jobs into the `mpi` partition.
The maximal running time is relatively high (14 days) to allow for longer jobs.
Don't abuse this.
Contact hpc-helpdesk@bihealth.de if you have longer running jobs that you really cannot make run any shorter for assistance.

For access to it you have register [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) (who will grant all requests).
See [Resource Registration: GPU Nodes](../admin/resource-registration.md#high-memory-nodes) for details.

* **maximal running time:** 14 days of maximal running time
* **partition name:** `highmem`
* **argument string:** `-p mpi`, maximal running time: `-t 14-00:00:00`
