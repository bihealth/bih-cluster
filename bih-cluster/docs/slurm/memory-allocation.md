# Memory Allocation

Memory allocation is one of the topics that users find confusing most often.
This section first gives some technical background and then explains how to implement this properly with Slurm on the BIH HPC.

## Technical Background

!!! note "Technical Background Summary"

    - **virtual memory** is what your programs tells the operating system it wants to use
    - **resident set size** is the amount of memory that your program *actually* uses
    - most memory will be allocated on the **heap**

Main memory used to be one of the most important topics when programming, as computers had so little.
There is the infamous quote "640KB ought ot be enough for anybody" wrongly attribute to Bill Gates which refers to the fact that early computers could only address that amount of memory. In MS DOS, one had to use special libraries for a program to use more memory.
Today, computers are very fast and memory is plentiful and people can (rightfully) forget about memory allocation ... as long as they don't use "much" memory by today's standards.

The Linux operating system differentiates between the following types of memory:

- **virtual memory size (vsize)**, the amount of memory that a process (virtually) allocates,
- **resident set size (rss)**, the amount of memory actually used and currently in the computer's main memory,
- the **swap memory usage**, the amount of active memory that is not present in main memory but on the computer's disk,
- sometimes, the **shared memory** is also interesting, and
- it might be interesting to know about **heap** and **stack** size.

Note that above we are talking about processes, not Slurm jobs yet.
Let us look at this in detail:

Each program uses some kind of memory management.
For example, in C the `malloc` and `free` functions manually allocate and free memory while in Java, R, and Python, memory allocation and release is done automatically using a concept called garbage collection.
Each program starts with a certain **virtual memory size**, that is the amount of memory it can address, say 128MB.
When the program allocates memory, the memory allocation mechanism will check whether it has sufficient space left.
If not, it will request an increase in virtual memory from the operating system, e.g., to 256MB.
If this fails then the program can try to handle the error, e.g., terminate gracefully, but many programs will just panic and stop.
Otherwise, the program will get access to more memory and happily continue to run.

However, programs can allocate humonguous amounts of virtual memory and only use a little.
Memory is organized in "pages" (classically these are 4096 bytes each, but can be larger using so-called "huge page" features).
The operating system tracks which memory pages are actually used by a process.
The total size of these pages is called the **resident set size**: the amount of memory that is actually currently used by a program.
Programs can also mark pages as unused again, thus freeing resident memory and can also decrease their virtual memory.

In some cases it is still interesting to use **swap memory**.
Here, the contents of resident memory are copied to disk by the operating system.
This process is completely transparent to the program; the data remains available at the original positions in the *virtual memory*!
However, accessing it will take some time as it must be read back into main memory from the disk.
In this way, it was possible for a computer with 4MB of RAM and a disk of 100MB to run programs that used 8MB.
Of course, this was only really useable for programs that ran in the background. One could really feel the latency if a graphical program was using swapped memory (you could actually hear the hard drive working).
Today, swap storage is normally only relevant when put your computer into hibernation.
Given the large main memory on the cluster nodes, their small local hard drives (just used for loading the operating system), and the extreme slowness involved in using swapped memory, the BIH HPC nodes have no swap memory allocated.

Most HPC users will also use **shared memory**, at least implicitly.
Whenever a program uses `fork` to create a subprocess (BTW, this is *not* a thread), the program can chose to "copy" its current address space.
The second process then has access to the same memory than the parent process in a *copy-on-write* fashion.
This allows, for example, pre-loading a database, and also allows the use of already loaded library code by the child process as well.
If the child process writes to the copy-on-write memory of the parent, the relevant memory page will be copied and attributed to the child.

Two or more processes can share the same memory explicitly. This is usually used for inter-process communication but the Bowtie program uses it for sharing the memory of indices.
For example, the Python `multiprocessing` module will use this, including if you have two MPI processes running on the same host.

Memory is also separated into segments, the most interesting ones are **heap** and **stack** memory.
For compiled languages, memory can be allocated on either.
For C, an `int` variable will be allocated on the stack.
Every time you call a function, a *stack frame* is created in memory to hold the local variables and other information for the duration of the function execution.
The stack thus grows through function calls made by your program and shrinks when the functions return.
The stack size for a process is limited (by `ulimit -s`) and a program that goes too deep (e.g., via infinite recursion) will be terminated by the operating system if it exceeds this limit.
Again in C, `int * ptr = (int *)malloc(10 * sizeof(int));` will allocate memory for one variable (an integer pointer) on the stack and memory for 10 integers on the heap.
When the function returns, the `ptr` variable on the stack will be freed but to free the array of integers, you'd have to call `free(ptr)`.
If the memory is not freed then this constitutes a *memory leak*, but that is another topic.

Other relevant segments are **code**, where the compiled code lives, and **data**, where static data such as strings displayed to the user are stored.
As a side node, in interpreted languages such as R or Python, the code and data segments will refer to the code and data of Python while the actual program text will be on the heap.

## Interlude: Memory in Java

!!! note "Memory in Java Summary"

    - set `-XX:MaxHeapSize=<size>` (e.g., `<size>=2G`) for your program and only tune the other parameters if needed
    - also consider the amount of memory that Java needs for heap management in your Slurm allocations

Java's memory management provides for some interesting artifacts.
When running simple Java programs, you will never run into this but if you need to use gigabytes of memory in Java then you will have to learn a bit about Java memory management.
This is the case when running [GATK](https://gatk.broadinstitute.org/hc/en-us) programs, for example.

As different operating systems handle memory management differently, the Java virtual machine does its own memory management to provide a consistent interface.
The following three settings are important in governing memory usage of Java:

- `-Xmx<size>`/`-XX:MaxHeapSize=<size>` -- the maximal *Java heap* size
- `-Xms<size>`/`-XX:InitialHeapSize=<size>` -- the initial *Java heap* size
- `-Xss<size>`/`-XX:ThreadStackSize=<size>` -- maximal stack size available to a Java thread (e.g., the main thread)

Above, `<size>` is a memory specification, either in bytes or with a suffix, e.g., `80M`, or `1G`.

On startup, Java does roughly the following:

- Setup the core virtual machine, load libraries, etc. and allocate (vsize) consume (rss) memory on the *OS (operating system) heap*.
- Setup the *Java heap* allocate (vsize) and consume (rss) memory on the *OS heap*.
  In particular, Java will need to setup data structures for the memory management of each individual object.
- Run the program where Java data and Java threads will lead to memory allocation (vsize) and consumption (rss) of memory.

Memory freed by the Java garbage collector can be re-used by other Java objects (rss remains the same) or be freed in the operating system (rss decreases).
The Java VM program itself will also consume memory on the *OS stack* but that is negligible.

Overall, the Java VM needs to store in main memory:

- The Java VM, program code, Java thread stacks etc. (very little memory).
- The *Java heap* (potentially a lot of memory).
- The *Java heap management* data structures (so-called "off-heap", but of course on the *OS heap*) (potentially also considerable memory).

In the BIH HPC context, the following is recommended to:

- Set the Java heap to an appropriate size (use trial-and-error to determine the correct size or look through internet forums).
- Only tune initial heap size in the case of performance issues (unlikely in batch processing).
- Only bump the stack size when problems occur.
- Consider "off-heap" memory when performing Slurm allocations.

## Memory Allocation in Slurm

!!! note "Memory Allocation in Slurm Summary"

    - most user will simply use `--memory=<size>` (e.g., `<size>=3G`) to allocate memory per node
    - both interactive `srun` and batch `sbatch` jobs are governed by Slurm memory allocation
    - the sum of all memory of all processes started by your job may not exceed the job reservation.
    - please don't over-allocate memory, see "Memory Accounting in Slurm" below for details

Our Slurm configuration uses Linux cgroups to enforce a maximum amount of resident memory.
You simply specify it using `--memory=<size>` in your `srun` and `sbatch` command.

In the (rare) case that you provide more flexible number of threads (Slurm tasks) or GPUs, you could also look into `--mem-per-cpu` and `--mem-per-gpu`.
The [official Slurm sbatch](https://slurm.schedmd.com/sbatch.html) manual is quite helpful, as is `man sbatch` on the cluster command line.

Slurm (or rather Linux via cgroups) will track all memory started by all jobs by your process.
If each process works independently (e.g., you put the output through a pipe `prog1 | prog2`) then the amount of memory consumed will at any given time be the sum of the RSS of both processes *at that time*.
If your program uses `fork`, which uses memory in a copy-on-write fashion, the shared memory is of course only counted once.
Note that Python's multiprocessing does **not** use copy on write: its data will be explicitly copied and consume additional memory.
Refer to the Scipy/Numpy/Pandas etc. documentation on how to achieve parallelism without copying too much data.

The amount of virtual memory that your program can reserve is only "virtually" unlimited (pun not intended).
However, in practice, the operating system will not like you allocating more than physically available.
If your program attempts to allocate more memory than requested via Slurm, your program will be killed.

This is reported to you in the Slurm job output log as something like:

```
slurmstepd: error: Detected 1 oom-kill event(s) in step <JOB ID>.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
```

You can inspect the amount of memory available on each node in total with `sinfo --format "%.10P %.10l %.6D %.6m %N"`, as shown below.

```bash
$ sinfo --format "%.10P %.10l %.6D %.6m %N"
 PARTITION  TIMELIMIT  NODES MEMORY NODELIST
    debug*    8:00:00    240 128722 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
    medium 7-00:00:00    240 128722 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
      long 28-00:00:0    240 128722 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
  critical 7-00:00:00    176 128722 med[0101-0164,0501-0516,0601-0632,0701-0764]
   highmem 14-00:00:0      4 515762 med[0401-0404]
       gpu 14-00:00:0      4 385215 med[0301-0304]
       mpi 14-00:00:0    240 128722 med[0101-0164,0201-0264,0501-0516,0601-0632,0701-0764]
```

## Memory/CPU Accounting in Slurm

!!! note "Memory Accounting in Slurm Summary"

    - you can use Slurm accounting to see memory and CPU usage of your program
    - use `sacct -j JOBID --format=JobID,MaxRSS` to display the RSS usage of your program
    - use `sacct -j JOBID --format=Elapsed,AllocCPUs,TotalCPU` to display information about CPU usage
    - consider using the helpful script below to compute overallocated memory

While Slurm runs your job, it collects information about the job such as the running time, exit status, and memory usage.
This information is available through the *scheduling* system via the `squeue` and `scontrol` commands, but only while the job is pending execution, executing, or currently completing.
After job completion, the information is only available through the **Slurm accounting** system.

You can query information about jobs, e.g., using `sacct`:

```bash
$ sacct -j 1607166
       JobID    JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
1607166      snakejob.+   critical                    16  COMPLETED      0:0
1607166.bat+      batch                               16  COMPLETED      0:0
1607166.ext+     extern                               16  COMPLETED      0:0
```

This shows that the job with ID `1607166` with a job ID starting with `snakejob.` has been run in the `critical` partition, been allocated 16 cores and had an exit code of `0:0`.
For technical reasons, there is a `batch` and an `extern` sub step.
Actually, Slurm makes it possible to run various steps in one batch [as documented in the Slurm documentation](https://slurm.schedmd.com/job_launch.html).

The `sacct` command has various command-line options that you can read about via `man sacct` or [in the Slurm documentation](https://slurm.schedmd.com/sacct.html).
We can use `--brief`/`-b` to show only a brief summary.

```bash
$ sacct -j 1607166 --brief
       JobID      State ExitCode
------------ ---------- --------
1607166       COMPLETED      0:0
1607166.bat+  COMPLETED      0:0
1607166.ext+  COMPLETED      0:0
```

Similarly, you can use `--long` to display extended information (see the manual for the displayed columns).
Very long report lines can be piped into `less -S` for easier display.
You can fine-tune the information to display with a format string to `--format`:

```bash
$ sacct -j 1607166 --format=JobID,ReqMem,MaxRSS,Elapsed,TotalCPU,AllocCPUS
       JobID   ReqMem     MaxRSS    Elapsed   TotalCPU  AllocCPUS
------------ --------- ---------- ---------- ---------- ----------
1607166          60Gn              13:07:31 7-16:21:29         16
1607166.bat+     60Gn   4314560K   13:07:31 7-16:21:29         16
1607166.ext+     60Gn          0   13:07:31  00:00.001         16
```

From this command, we can read that we allocate 60GB memory of memory per node (suffix `n`, here `Gn` for gigabytes per node) and the maximum RSS is reported as 4.3GB.
You can use this information to fine-tune your memory allocations.
As a side-remark, a suffic `c` indicates the memory per core (e.g., that could be`60Gc`)

Further, the program ran for 13 hours and 7 minutes with allocated 16 CPU cores and consumed a total of 7 days, 16 hours, and 21 minutes of CPU time.
Thus, a total of 10,061 CPU minutes were spent in 787 minutes wall-clock time.
This yields an overall empirical degree of parallelism of about 10061 / 787 = 14, and a parallel efficiency of 14 / 16 = 88%.
The discussion of parallel efficiency is a topic not covered here.

However, you can use the `awk` script below to compute the empirical parallelism (`EmpPar`) and the parallel efficiency (`ParEff`).
The script also displays the difference I requested, and used RSS (`DiffRSS`).
The script [can be found here](assets/quick-sacct.awk).

```bash
$ sacct -j 1607166 --format=JobID,ReqMem,MaxRSS,Elapsed,TotalCPU,AllocCPUS \
                | awk -f quick-sacct.awk
       JobID     ReqMem     MaxRSS    Elapsed   TotalCPU  AllocCPUS     EmpPar   ParEff  DiffMEM
------------ ---------- ---------- ---------- ---------- ----------  --------- -------- --------
1607166            60Gn              13:07:31 7-16:21:29         16      0.00     0.00        -
1607166.bat+       60Gn   4314560K   13:07:31 7-16:21:29         16     14.05     0.88    55.89
1607166.ext+       60Gn          0   13:07:31  00:00.001         16      0.00     0.00        -
```
