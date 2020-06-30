# Background In: Job Schedulers

Assuming that you have read [Background In: Computers](bg-computers.md) and [Background In: Systems](bg-systems.md), how does this relate to the BIH HPC?
What is "a cluster" in general and "the BIH HPC" specifically?

Computers can be tied together ("clustered") for various reasons.
The aim of high-performance clusters is to provide computing power from multiple machines to provide an advantage over "multiple individual machines".
Generally, the components are:

- individual CPUs with their memory and local disks,
- a shared file system for exchanging data, and
- a scheduling system to allow for sharing the system fairly between multiple users.

## Nodes

Each individual system (that is tied together with other systems) is called a **node**.
Such nodes are separate systems, having one or multiple (individual) CPUs, memory (RANM), and disks, that are physically separate from other node's hardware.
Each node connnects to other nodes (and the storage system) via the network.

## Storage

While each node provides (limited) local storage, the centralized storage offers terabytes, if not petabytes, of storage instead of gigabytes.
This centralized storage is **shared** which has positive implications (larger than local storage, can be accessed by multiple nodes at high throughput) and negative implications (latency and single-node throughput might be lower than for local storage and is dependent on other nodes).
In the BIH HPC, `/fast` (which is the same as `/data/gpfs-1`) is the shared storage while `/tmp` allows access to the local storage.
Also, depending on the hardware, the local storage might differ between media (spinning disk vs. SSH) while the central storage is independent of the local storage system.

## Jobs

The individual machines (that are connected by network and also via network by the central storage) are connected by the **scheduling system (Slurm)**.
Jobs specifications are submitted to the system. 
These specifications consist of resource requirements (e.g., wall-clock time, number of CPUs, amount of RAM, GPUs, if any).

The scheduler (Slurm) will then allocate "chunks" in the given "resource space" of requirements.
These jobs will be placed on a number of criteria, e.g.:

- Small jobs (e.g., few CPUs, GPUs, and memory) will be favoured over large jobs.
- Jobs from users consuming few resources overall will be favoured over users with high cluster usage.
- Jobs from high-priority partitions (e.g., `critical`) will be favoured over jobs from low-priority partitions (e.g., `debug` or `medium`).

## Scheduling

The Slurm scheduler uses refined/advanced algorithms for deciding which job to schedule.
You can read the [Priority Multifactor](https://slurm.schedmd.com/priority_multifactor.html) documentation for details.
Most importantly, read about [`sprio`](https://slurm.schedmd.com/priority_multifactor.html#sprio) on how to inspect your jobs for their current prioritiy rating.

Please note that job scheduling and prioritization is a bit of a "black art" and subject to constant optimization.
Don't try to play the admins as they are more powerful than you (they have root ;).
On the other hand, the admins promise not "to play users" and aim for providing "fair share" of the system (as the admins consider "the user in general"/"the user collective") to be smarter that them (we also consider you to be **collaborative** ;).
