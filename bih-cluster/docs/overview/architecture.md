# Cluster Architecture

BIH HPC IT provides acess to high-performance compute (HPC) cluster systems.
A cluster system bundles a high number of nodes and in the case of HPC, the focus is on performance (with contrast to high availability clusters).

## HPC 4 Research

=== "HPC 4 Research"

    ### Cluster Hardware

    - approx. 256 nodes (from three generations),
    - 4 high-memory nodes (2 nodes with 512 GB RAM, 2 nodes with 1 TB RAM),
    - 4 GPU nodes (with 4 Tesla GPUs each), and
    - a high-perfomance parallel GPFS files system.

    ### Network Interconnect

    - All nodes are connected with 2x10GbE,
    - 32 nodes provide Infiniband interconnect (lower latency, but MPI library required).

=== "HPC 4 Clinic"

    ### Cluster Hardware

    - 24 general purpose nodes,
    - 3 GPU nodes (with 4 Tesla GPUs each), and
    - a scale-out NAS system based on Isilon

    ### Network Interconnect

    - All nodes are connected with 25GbE.

## Cluster Management

Users don't connect to nodes directly but rather create interactive or batch jobs to be executed by the cluster job scheduler *Slurm*.

- **Interactive jobs** open interactive sessions on compute nodes (e.g., R or iPython sessions).
  These jobs are run directly in the user's terminal.
- **Batch jobs** consist a job script with execution instructions (a name, resource requirements etc.)
  These are submitted to the cluster and then assigned to compute hosts by the job scheduler.
  Users can configure the scheduler to send them an email upon completion.
  Users can submit many batch jobs at the same time and the scheduler will execute them once the cluster offers sufficient resources.
- **Web-based access** can be achieved using the [OnDemand Portal](/ondemand/overview/)

## Head vs. Compute Nodes

As common with HPC systems, users cannot directly access the compute nodes but rather connect to so-called *head nodes*.
The BIH HPC system provides the following head nodes:

- `login-1` and `login-2` that accept SSH connections and are meant for **low intensity**, interactive work such as editing files, running screen/tmux sessions, and logging into the compute nodes.
  Users should run **no computational tasks** and **no large-scale data transfer** on these nodes.
- `transfer-1` and `transfer-2` also accept SSH connections.
  Users should **run all large-scale data transfer through these nodes.**

## Common Use Case

After registration and client configurations, users with typically connect to the HPC system through the login nodes:

```bash
local:~$ ssh -l jdoe_c login-1.research.hpc.bihealth.org
res-login-1:~$
```

Subsequently, they might submit batch jobs to the cluster for execution through the Slurm scheduling system or open interactive sessions:

```bash
res-login-1:~$ sbatch job_script.sh
res-login-1:~$ srun --pty bash -i
med0104:~$
```
