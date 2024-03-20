# How-To: Build and Run OpenMPI Program

This article describes how to build an run an OpenMPI program.
We will build a simple C program that uses the OpenMPI message passing interface and run it in parallel.
You should be able to go from here with other languages and more complex programs.
We will use a simple Makefile for building the software.

## Loading OpenMPI Environment

First, load the OpenMPI package.

```bash
hpc-login-1:~$ srun --pty bash -i
med0127:~$ module load openmpi/4.3.0-0
```

Then, check that the installation works

```bash
med0127:~$ ompi_info | head
                 Package: Open MPI root@med0127 Distribution
                Open MPI: 4.0.3
  Open MPI repo revision: v4.0.3
   Open MPI release date: Mar 03, 2020
                Open RTE: 4.0.3
  Open RTE repo revision: v4.0.3
   Open RTE release date: Mar 03, 2020
                    OPAL: 4.0.3
      OPAL repo revision: v4.0.3
       OPAL release date: Mar 03, 2020
```

## Building the example

Next, clone the OpenMPI example project from Gitlab.

```terminal
med0127:~$ git clone git@github.com:bihealth/bih-cluster-openmpi-example.git
med0127:~$ cd bih-cluster-openmpi-example/src
```

`Makefile`

```Makefile
.PHONY: default clean

# configure compilers
CC=mpicc
CXX=mpicxx
# configure flags
CCFLAGS += $(shell mpicc --showme:compile)
LDFLAGS += $(shell mpicc --showme:link)

default: openmpi_example

openmpi_example: openmpi_example.o

clean:
	rm -f openmpi_example.o openmpi_example
```

`openmpi_example.c`

```c++
#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
```

`run_mpi.sh`

```shell
#!/bin/bash

# Example job script for (single-threaded) MPI programs.

# Generic arguments

# Job name
#SBATCH --job-name openmpi_example
# Maximal running time of 10 min
#SBATCH --time 00:10:00
# Allocate 1GB of memory per node
#SBATCH --mem 1G
# Write logs to directory "slurm_log"
#SBATCH -o slurm_log/slurm-%x-%J.log

# MPI-specific parameters

# Run 64 tasks (threads/on virtual cores)
#SBATCH --nodes 64

# Make sure to source the profile.d file (not available on head nodes).
/etc/profile.d/modules.sh

# Load the OpenMPI environment module to get the runtime environment.
module load openmpi/3.1.0-0

# Launch the program.
mpirun -np 64 ./openmpi_example
```

The next step is building the software

```bash
med0127:~$ make
mpicc    -c -o openmpi_example.o openmpi_example.c
mpicc -pthread -Wl,-rpath -Wl,/opt/local/openmpi-4.0.3-0/lib -Wl,--enable-new-dtags -L/opt/local/openmpi-4.0.3-0/lib -lmpi  openmpi_example.o   -o openmpi_example
med0127:~$ ls -lh
total 259K
-rw-rw---- 1 holtgrem_c hpc-ag-cubi  287 Apr  7 23:29 Makefile
-rwxrwx--- 1 holtgrem_c hpc-ag-cubi 8.5K Apr  8 00:15 openmpi_example
-rw-rw---- 1 holtgrem_c hpc-ag-cubi  760 Apr  7 23:29 openmpi_example.c
-rw-rw---- 1 holtgrem_c hpc-ag-cubi 2.1K Apr  8 00:15 openmpi_example.o
-rwxrwx--- 1 holtgrem_c hpc-ag-cubi 1.3K Apr  7 23:29 run_hybrid.sh
-rwxrwx--- 1 holtgrem_c hpc-ag-cubi  663 Apr  7 23:35 run_mpi.sh
drwxrwx--- 2 holtgrem_c hpc-ag-cubi 4.0K Apr  7 23:29 sge_log
```

The software will run outside of the MPI environment -- but in a single process only, of course.

```bash
med0127:~$ ./openmpi_example
Hello world from processor med0127, rank 0 out of 1 processors
```

## Running OpenMPI Software

All of the arguments are already in the `run_mpi.sh` script.

```
med01247:~# sbatch run_mpi.sh
```

Explanation of the OpenMPI-specific arguments

- `--ntasks 64`: run 64 processes in the MPI environment.

Let's look at the slurm log file, e.g., in `slurm_log/slurm-openmpi_example-3181.log`.

```
med0124:~$  cat slurm_log/slurm-openmpi_example-*.log
Hello world from processor med0133, rank 6 out of 64 processors
Hello world from processor med0133, rank 25 out of 64 processors
Hello world from processor med0133, rank 1 out of 64 processors
Hello world from processor med0133, rank 2 out of 64 processors
Hello world from processor med0133, rank 3 out of 64 processors
Hello world from processor med0133, rank 7 out of 64 processors
Hello world from processor med0133, rank 9 out of 64 processors
Hello world from processor med0133, rank 12 out of 64 processors
Hello world from processor med0133, rank 13 out of 64 processors
Hello world from processor med0133, rank 15 out of 64 processors
Hello world from processor med0133, rank 16 out of 64 processors
Hello world from processor med0133, rank 17 out of 64 processors
Hello world from processor med0133, rank 18 out of 64 processors
Hello world from processor med0133, rank 23 out of 64 processors
Hello world from processor med0133, rank 24 out of 64 processors
Hello world from processor med0133, rank 26 out of 64 processors
Hello world from processor med0133, rank 27 out of 64 processors
Hello world from processor med0133, rank 31 out of 64 processors
Hello world from processor med0133, rank 0 out of 64 processors
Hello world from processor med0133, rank 4 out of 64 processors
Hello world from processor med0133, rank 5 out of 64 processors
Hello world from processor med0133, rank 8 out of 64 processors
Hello world from processor med0133, rank 10 out of 64 processors
Hello world from processor med0133, rank 11 out of 64 processors
Hello world from processor med0133, rank 14 out of 64 processors
Hello world from processor med0133, rank 19 out of 64 processors
Hello world from processor med0133, rank 20 out of 64 processors
Hello world from processor med0133, rank 21 out of 64 processors
Hello world from processor med0133, rank 22 out of 64 processors
Hello world from processor med0133, rank 28 out of 64 processors
Hello world from processor med0133, rank 29 out of 64 processors
Hello world from processor med0133, rank 30 out of 64 processors
Hello world from processor med0134, rank 32 out of 64 processors
Hello world from processor med0134, rank 33 out of 64 processors
Hello world from processor med0134, rank 34 out of 64 processors
Hello world from processor med0134, rank 38 out of 64 processors
Hello world from processor med0134, rank 39 out of 64 processors
Hello world from processor med0134, rank 42 out of 64 processors
Hello world from processor med0134, rank 44 out of 64 processors
Hello world from processor med0134, rank 45 out of 64 processors
Hello world from processor med0134, rank 46 out of 64 processors
Hello world from processor med0134, rank 53 out of 64 processors
Hello world from processor med0134, rank 54 out of 64 processors
Hello world from processor med0134, rank 55 out of 64 processors
Hello world from processor med0134, rank 60 out of 64 processors
Hello world from processor med0134, rank 62 out of 64 processors
Hello world from processor med0134, rank 35 out of 64 processors
Hello world from processor med0134, rank 36 out of 64 processors
Hello world from processor med0134, rank 37 out of 64 processors
Hello world from processor med0134, rank 40 out of 64 processors
Hello world from processor med0134, rank 41 out of 64 processors
Hello world from processor med0134, rank 43 out of 64 processors
Hello world from processor med0134, rank 47 out of 64 processors
Hello world from processor med0134, rank 48 out of 64 processors
Hello world from processor med0134, rank 49 out of 64 processors
Hello world from processor med0134, rank 50 out of 64 processors
Hello world from processor med0134, rank 51 out of 64 processors
Hello world from processor med0134, rank 52 out of 64 processors
Hello world from processor med0134, rank 56 out of 64 processors
Hello world from processor med0134, rank 57 out of 64 processors
Hello world from processor med0134, rank 59 out of 64 processors
Hello world from processor med0134, rank 61 out of 64 processors
Hello world from processor med0134, rank 63 out of 64 processors
Hello world from processor med0134, rank 58 out of 64 processors
```

## Running Hybrid Software (MPI+Multithreading)

In some cases, you want to mix multithreading (e.g., via OpenMP) with MPI to run one process with multiple threads that then can communicate via shared memory.
Note that OpenMPI will let processes on the same node communicate via shared memory anyway, so this might not be necessary in all cases.

The file `run_hybrid.sh` shows how to run an MPI job with 8 threads each.

> **Note well that memory is allocated on a per-slot (thus per-thread) base!**

`run_hybrid.sh`

```
#!/bin/bash

# Example job script for multi-threaded MPI programs, sometimes
# called "hybrid" MPI computing.

# Generic arguments

# Job name
#SBATCH --job-name openmpi_example
# Maximal running time of 10 min
#SBATCH --time 00:10:00
# Allocate 1GB of memory per node
#SBATCH --mem 1G
# Write logs to directory "slurm_log"
#SBATCH -o slurm_log/slurm-%x-%J.log

# MPI-specific parameters

# Run 8 tasks (threads/on virtual cores)
#SBATCH --ntasks 8
# Allocate 4 CPUs per task (cores/threads)
#SBATCH --cpus-per-task 4

# Make sure to source the profile.d file (not available on head nodes).
source /etc/profile.d/modules.sh

# Load the OpenMPI environment module to get the runtime environment.
module load openmpi/4.0.3-0

# Launch the program.
mpirun -n 8 ./openmpi_example
```

We changed the following

- run 8 tasks ("processes")
- allocate 4 threads each

Let's look at the log output:

```terminal
# cat slurm_log/slurm-openmpi_example-3193.log
Hello world from processor med0133, rank 1 out of 8 processors
Hello world from processor med0133, rank 3 out of 8 processors
Hello world from processor med0133, rank 2 out of 8 processors
Hello world from processor med0133, rank 6 out of 8 processors
Hello world from processor med0133, rank 0 out of 8 processors
Hello world from processor med0133, rank 4 out of 8 processors
Hello world from processor med0133, rank 5 out of 8 processors
Hello world from processor med0133, rank 7 out of 8 processors
```

Each process can now launch 4 threads (e.g., by defining `export OMP_NUM_THREADS=4` before the program call).