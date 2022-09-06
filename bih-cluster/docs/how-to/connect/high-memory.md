# How-To: Connect to High-Memory Nodes

## Prequisites

You have to register with [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) for requesting access.

Afterwards, you can connect to the High-Memory using the `highmem` SLURM partition (see below).
Jobs allocating more than 200GB of RAM should be routed automatically to the `highmem` nodes.

## How-To

In the cluster there are four High-memory used which can be used:

```
res-login-1:~$ sinfo -p highmem
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
highmem      up 14-00:00:0      3   idle med040[1-4] 
```

To connect to one of them, simply allocate more than 200GB of RAM in your job.

```
res-login-1:~$ srun --pty --memory=300GB bash -i
med0401:~$
```

You can also pick one of the hostnames:

```
res-login-1:~$ srun --pty --memory=300GB --nodelist=med0403 bash -i
med0403:~$
```

After successfull login, you can see that you are in "highmem" queue:

```
med0403:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
[...]
               270   highmem     bash holtgrem  R       1:25      1 med0403 

```
