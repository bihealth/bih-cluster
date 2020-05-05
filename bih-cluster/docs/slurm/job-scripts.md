# Slurm Job Scripts

This page describes how to create SLURM job scripts.

SLURM job scripts look as follows.
On the top you have lines starting with `#SBATCH`.
These appear as comments to bash scripts.
These lines are interpreted by `sbatch` in the same way as command line arguments.
That is, when later submitting the script with `sbatch my-job.sh` you can either have the parameter to the `sbatch` call or in the file.

!!! note "Multi-Node Allocation in Slurm"

    Classically, jobs on HPC systems are written in a way that they can run on multiple nodes at once, using the network to communicate.
    Slurm comes from this world and when allocating more than one CPU/core, it might allocate them on different nodes.
    Please use `--nodes=1` to force Slurm to allocate them on a single node.

**Creating the Script**

```bash
host:example$ cat >my-job.sh <<"EOF"
#!/bin/bash
#
#SBATCH --job-name=this-is-my-job
#SBATCH --output=output.txt
#
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100M

date

hostname
>&2 echo "Hello World"

sleep 1m

date
EOF
```

Also see the [SLURM Rosetta Stone](rosetta-stone.md) for more options.

**Submit, Look at Queue & Result**

```
host:example$ sbatch script.sh 
Submitted batch job 315
host:example$ squeue  -u holtgrem_c
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
               315     debug this-is- holtgrem  R       0:40      1 med0127 
host:example$ sleep 2m
host:example$ squeue  -u holtgrem_c
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
host:example$ cat output.txt 
Wed Mar 25 13:30:56 CET 2020
med0127
Hello World
Wed Mar 25 13:31:56 CET 2020
```