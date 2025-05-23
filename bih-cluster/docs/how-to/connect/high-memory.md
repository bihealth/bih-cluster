# How-To: Connect to High-Memory Nodes

The cluster has 4 high-memory nodes with 0.75 (2x), 1.5 (2x), and 1x4 (1x) TB of RAM.
You can connect to these nodes using the `highmem` SLURM partition (see below).
Jobs allocating more than 200 GB of RAM are automatically routed to the `highmem` nodes.

!!! info
    Fair use rules apply.
    As high-memory nodes are a limited resource, excessive use by single users is prohibited and can lead to mitigating actions.
    Be nice and cooperative with other users.
    Tip: `getent passwd USER_NAME` will give you a user's contact details.

## How-To

In the cluster there are four High-memory used which can be used:

```
hpc-login-1:~$ sinfo -p highmem
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST 
highmem      up 14-00:00:0      5   idle hpc-mem-[1-5] 
```

To connect to one of them, simply allocate more than 200GB of RAM in your job.

```
hpc-login-1:~$ srun --pty --mem=300GB bash -i
hpc-mem-1:~$
```

You can also pick one of the hostnames:

```
hpc-login-1:~$ srun --pty --mem=300GB --nodelist=hpc-mem-2 bash -i
hpc-mem-2:~$
```

After successfull login, you can see that you are in "highmem" queue:

```
med0403:~$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
[...]
               270   highmem     bash holtgrem  R       1:25      1 hpc-mem-2 

```
