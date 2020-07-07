# Introduction to Scheduling

As explained elsewhere in more detail, an HPC cluster consists of multiple computers that are connected via network and work together.
Multiple users can use the system to do their work.
This means that the system needs to **join** multiple computers (*nodes*) to provide a coherent view on them and the same time **partition** the system to allow for multiple users working concurrently.

```

    user 1         user 2          ...
  
  .---. .---.    .---. .---.
  | J | | J |    | J | | J |
  | o | | o |    | o | | o |       ...
  | b | | b |    | b | | b |
  | 1 | | 2 |    | 3 | | 4 |
  '---' '---'    '---' '---'

.------------------------------------------.
|            Cluster Scheduler             |
'------------------------------------------'

.----------.  .------------.  .------------.
| multiple |  |  separate  |  | computers  |
'----------'  '------------'  '------------'

```

## Interlude: Partitioning Single Computers

Overall, this partitioning is not so much different from how your workstation or laptop works.
Most likely, your computer (or even your smartphone) has multiple processors (or cores).
You can run multiple programs on the same computer and the fact that (a) there is more than one core and (b) there is more than one program running is opaque to the running programs (except if they explicitely communicate with each other).
The programs can explicitely take advantage of the multiple processor cores
The main difference is that you use your computer in an interactive fashion (you perform an action and expect an immediate reaction).

Even with a single processor (and core), your computer manages to run more than program at the same time.
This is done with the so-called *time-slicing* approach where the operating system lets each programs run in turn for a short time (a few milliseconds).
A program with a higher priority might get more time slices than one with a lower (e.g., your audio player has soft real-time requirements and you would hear artifacts if it is starved for compute resources).
Your operating system protects each program from the next by creating an *address space* for each program.
When two programs are running, the value of the memory at any given position in one program is independent from the value in the other program.
Your operating system offers explicit functionality for sharing certain memory areas that two programs can use to exchange data efficiently.

Similarly, file permissions with Unix users/groups or Unix/Windows ACLs (access control lists) are used to isolate users from each other.
Programs can share data by accessing the same file if they can both access it.
There are special files called *sockets* that allow for network-like inter-process communication but of course two programs on the same computer can also connect (virtually) via the computer network (no data will actually go through a cable).

## Interlude: Resource Types

As another diversion, let us consider how Unix manages its resources.
This fact is important to understand how to request resources from the scheduler later on.

First of all, a computer might offer a certain **feature** such as a specific hardware platform or special network connection.
Examples for this on the BIH HPC are specific Intel processor generations such as `haswell` or the availability of Infiniband networking.
You can request these with so-called constraints; they are not allocated to specific jobs.

Second, there are resources that are allocated to specific jobs.
The most important resources here are:

- computing resources (processors/CPUs (central progressing units) and cores, details are explained below),
- main memory / RAM, and
- special hardware such as GPUs, and
- (wall-clock) time that a job wants to run as an upper bound.

Generally, once a resource has been allocated to a job, they are not available to another.
This means if you allocating more resources to your job that you actually need (*overallocation*) then these resources are not available to other jobs (whether they are from you or other users).
This will be explained further below.

Another example for resource allocation are *licenses*.
The BIH HPC has a few Matlab 2016b licenses that users can request.
As long as a license is allocated to one job, they are unavailable to another.

## Nodes, Sockets, Processors, Cores, Threads

Regarding compute resources, Slurm differentiates between:

- nodes: a compute server,
- sockets: a socket in the compute server that hosts one physical processor,
- processor: a CPU or a *CPU core* in a multi-core computer (all CPUs in the BIH HPC are multi-core), and
- (hardware) threads: most Intel CPUs feature hardware threads (also known as "hyperthreading") where each each core looks like two cores.

In most cases, you will use one compute node only.
For using two nodes, you would have to use message passing, e.g., MPI, for processes on more than one node to communicate.
Thus, you would mostly use single-threaded and multi-threaded processes or multiple processes on the same node.

![Slurm's Multiprocessing Definition](figures/slurm-multicore.gif){: .center}

*Above: Slurm's nomenclature for sockets, processors, cores, and threads (from [Slurm Documentation](https://slurm.schedmd.com/mc_support.html#defs)).*{: .center}

Co-locating processes/threads on the same socket has certain implications that are mostly useful for numerical applications.
We will not further go into detail here.
Slurm provides many different features of ways to specify allocation of "pinning" to specific process locations.
If you need this feature, we trust that you find sufficient explanation in the Slurm documentation.

Usually, you would allocate multiple cores (which Slurm uses synonymously with processors) on a single node (allocation on a single node is the default).

## How Scheduling Works

Slurm is an acronym for "Simple Linux Unix Resource Manager" and the word "scheduler" does not occur.
Actually, one classically differentiate managing of resources and scheduling.
The **resource manager** would allocate the user resources according to request for their job and ensure that there are no conflicts.
In the case of no resources being available, the **scheduler** will put the user's job into a queue and then distribute compute sources to them as they become available.
In the following, both is simplified as the term *scheduler*.

Also, the interesting case occurs if there currently not enough resources available for at least two jobs that are submitted to the scheduler.
The scheduler has to decide how to proceed.
When simplifying to only scheduling cores, each job will request a number of cores and a number of cores.
The scheduler will then generate a schedule that might look as follows.

```
core
  ^
4 |   |---job2---|
3 |   |---job2---|
2 |   |---job2---|
1 | |--job1--|
  +--------------------------> t time
       5    1    1    2
	        0    5    0
```

`job1` has allocated one core and `job2` has allocated two cores.
When a `job3` for one core is submitted at *t = 5* then it has to wait at least as long until `job1` is finished.
If `job3` would allocate two or more cores, it would have to wait at least until `job2` has finished as well.

We can now ask several questions, including the following:

1. *What if a job runs shorter than the allocated time?*
    -- In this case, resources become free and the scheduler will attempt to select the next job(s) to run.
2. *What if a job runs longer than the allocated time?*
    -- In this case, the scheduler will send a [Unix signal](https://en.wikipedia.org/wiki/Signal_(IPC)) to the process first.
    The job will be given a bit more time and then it will be terminated.
    You will find a note about this at the end of your job log file.
3. *What if multiple jobs compete for resources?*
    -- The scheduler will prefer certain jobs over certain others using the [Slurm Multifactor Priority Plugin](https://slurm.schedmd.com/priority_multifactor.html).
    In practice, small jobs will be preferred over large jobs, users with few used resources in the last month will be favored over heavy consumers, long-waiting jobs will be favored over recently submitted jobs, and many other factors.
    You can use the [sprio utility](https://slurm.schedmd.com/priority_multifactor.html#sprio) for inspecting the factors in real-time.
4. *How does the scheduler handle new requests?*
    -- Generally, the scheduler will add new jobs into the waiting queue.
    The scheduler runs priority and scheduler recomputations with different computational effort.
    Slurm is configured to perform computationally simple schedule recomputations quite regularly and larger recomputations will be done more infrequently.

Also see the [Slurm Frequently Asked Questions](https://slurm.schedmd.com/faq.html).

Please note that even if all jobs were known from the beginning then scheduling is a so-called *NP-complete* problem and whole computer science journals and books are dedicated only to scheduling.
Things get more complex in the case of online scheduling where jobs can appear at any time.
In practice, Slurm does a fantastic job with its heuristics but it heavily relies on parameter tuning.
HPC administration is constantly working on optimizing the scheduler settings.
In the case that you observe inexplicable behaviour, please notify us at hpc-helpdesk@bihealth.de.

## Slurm Partitions

In Slurm, the nodes of a cluster are split into **partitions**.
Nodes are assigned to one or more partition (see the [Job Scheduler](../overview/job-scheduler.md) section for details).
Jobs can also be assigned to one or more partitions and are executed in nodes of the given partition.

In the BIH HPC, partitions are used to stratify jobs of certain running times and provide different quality of service (e.g., maximal number of CPU cores available to a user for jobs of a certain running time and size).
The partitions `gpu` and `highmem` provide special hardware (the nodes are not assigned to other partitions) and the partition `mpi` allows MPI-parallelism and the allocation of jobs to more than one node.
The [Job Scheduler](../overview/job-scheduler.md) provides further details.
