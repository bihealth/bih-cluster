# How-To: Connect to High-Memory Nodes

!!! info "HPC 4 Research Only"

    High memory nodes are available on the HPC 4 Research cluster only.
    However, the HPC 4 Clinic cluster nodes have 350GB+ of RAM which is likely enough for many applications.

## Prequisites

You have to register with [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) for requesting access.

Afterwards, you can connect to the High-Memory using the `highmem` SLURM partition (see below).
Jobs allocating more than 200GB of RAM will be routed automatically to the `highmem` nodes.

## How-To

In the cluster there are four High-memory used which can be used:

```
med-login1:~$ sinfo -p highmem
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
highmem      up 14-00:00:0      3   idle med040[1-4] 
```

To connect to one of them, simply allocate more than 200GB of RAM in your job.

```
med-login1:~$ srun --pty --memory=300GB bash -i
med0401:~$
```

You can also pick one of the hostnames:

```
med-login1:~$ srun --pty --memory=300GB --nodelist=med0403 bash -i
med0403:~$
```

After successfull login, you can see that you are in "highmem" queue:

```
med0403:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
[...]
               270   highmem     bash holtgrem  R       1:25      1 med0403 

```
