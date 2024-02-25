# How-To: Connect to GPU Nodes

The cluster has seven nodes with four Tesla V100 GPUs each: `hpc-gpu-{1..7}` and one node with 10 A400 GPUs: `hpc-gpu-8`.

Connecting to a node with GPUs is easy.
You simply request a GPU using the `--gres=gpu:$CARD:COUNT` (for `CARD=tesla` or `CARD=a40`) argument to `srun` and `batch`.
This will automatically place your job in the `gpu` partition (which is where the GPU nodes live) and allocate a number of `COUNT` GPUs to your job.

!!! note

    Recently, `--gres=gpu:tesla:COUNT` was often not able to allocate the right partion on it's own.
    If scheduling a GPU fails, consider additionally indicating the GPU partion explicitely with `--partition gpu` (or `#SBATCH --partition gpu` in batch file).

!!! hint

    Make sure to read the FAQ entry "[I have problems connecting to the GPU node! What's wrong?](../../help/faq.md#i-have-problems-connecting-to-the-gpu-node-whats-wrong)".

!!! important "Interactive Use of GPU Nodes is Discouraged"

    While interactive computation on the GPU nodes is convenient, it makes it very easy to forget a job after your computation is complete and let it run idle.
    While your job is allocated, it blocks the **allocated** GPUs and other users cannot use them although you might not be actually using them.
    Please prefer batch jobs for your GPU jobs over interactive jobs.

    Further, interactive GPU jobs are currently limited to 24 hours.
    We will monitor the situation and adjust that limit to optimize GPU usage and usability.
    
!!! important "Allocation of GPUs through Slurm is mandatory"

    In other word: using GPUs from SSH sessions is prohibited.
    The scheduler is not aware of manually allocated GPUs and this interferes with other users' jobs.

## Prequisites

You have to register with [hpc-helpdesk@bih-charite.de](mailto:hpc-helpdesk@bih-charite.de) for requesting access.
Afterwards, you can connect to the GPU nodes as shown below.

## Preparation

We will setup a miniconda installation with `pytorch` testing the GPU.
If you already have this setup then you can skip this step

```bash
res-login-1:~$ srun --pty bash
med0703:~$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
med0703:~$ bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
med0703:~$ source ~/miniconda3/bin/activate
med0703:~$ conda create -y -n gpu-test pytorch cudatoolkit=10.2 -c pytorch
med0703:~$ conda activate gpu-test
med0703:~$ python -c 'import torch; print(torch.cuda.is_available())'
False
med0703:~$ exit
res-login-1:~$
```

The `False` shows that CUDA is not available on the node but that is to be expected.
We're only warming up!

## Allocating GPUs

Let us now allocate a GPU.
The Slurm schedule will properly allocate GPUs for you and setup the environment variable that tell CUDA which devices are available.
The following dry run shows these environment variables (and that they are not available on the login node).

```bash
res-login-1:~$ export | grep CUDA_VISIBLE_DEVICES
res-login-1:~$ srun --gres=gpu:tesla:1 --pty bash
med0303:~$ export | grep CUDA_VISIBLE_DEVICES
declare -x CUDA_VISIBLE_DEVICES="0"
med0303:~$ exit
res-login-1:~$ srun --gres=gpu:tesla:2 --pty bash
med0303:~$ export | grep CUDA_VISIBLE_DEVICES
declare -x CUDA_VISIBLE_DEVICES="0,1"
```

You can see that you can also reserve multiple GPUs.
If we were to open two concurrent connections (e.g., in a `screen`) to the same node when allocating one GPU each, the allocated GPUs would be non-overlapping.
Note that any two jobs are isolated using Linux cgroups ("container" technology 🚀) so you cannot accidentally use a GPU of another job.

Now to the somewhat boring part where we show that CUDA actually works.

```bash
res-login-1:~$ srun --gres=gpu:tesla:1 --pty bash
hpc-gpu-1:~$ nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2019 NVIDIA Corporation
Built on Wed_Oct_23_19:24:38_PDT_2019
Cuda compilation tools, release 10.2, V10.2.89
hpc-gpu-1:~$ source ~/miniconda3/bin/activate
hpc-gpu-1:~$ conda activate gpu-test
hpc-gpu-1:~$ python -c 'import torch; print(torch.cuda.is_available())'
True
```

!!! note

    Recently, `--gres=gpu:tesla:COUNT` was often not able to allocate the right partion on it's own.
    If scheduling a GPU fails, consider additionally indicating the GPU partion explicitely with `--partition gpu` (or `#SBATCH --partition gpu` in batch file).

## Bonus #1: Who is using the GPUs?

Use `squeue` to find out about currently queued jobs (the `egrep` only keeps the header and entries in the `gpu` partition).

```bash
res-login-1:~$ squeue | egrep -iw 'JOBID|gpu'
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                33       gpu     bash holtgrem  R       2:26      1 hpc-gpu-1
```

## Bonus #2: Is the GPU running?

To find out how active the GPU nodes actually are, you can connect to the nodes (without allocating a GPU you can do this even if the node is full) and then use `nvidia-smi`.

```bash
res-login-1:~$ ssh hpc-gpu-1 bash
hpc-gpu-1:~$ nvidia-smi
Fri Mar  6 11:10:08 2020
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 440.33.01    Driver Version: 440.33.01    CUDA Version: 10.2     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla V100-SXM2...  Off  | 00000000:18:00.0 Off |                    0 |
| N/A   62C    P0   246W / 300W |  16604MiB / 32510MiB |     99%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla V100-SXM2...  Off  | 00000000:3B:00.0 Off |                    0 |
| N/A   61C    P0   270W / 300W |  16604MiB / 32510MiB |    100%      Default |
+-------------------------------+----------------------+----------------------+
|   2  Tesla V100-SXM2...  Off  | 00000000:86:00.0 Off |                    0 |
| N/A   39C    P0    55W / 300W |      0MiB / 32510MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   3  Tesla V100-SXM2...  Off  | 00000000:AF:00.0 Off |                    0 |
| N/A   44C    P0    60W / 300W |      0MiB / 32510MiB |      4%      Default |
+-------------------------------+----------------------+----------------------+
                                                                               
+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID   Type   Process name                             Usage      |
|=============================================================================|
|    0     43461      C   python                                     16593MiB |
|    1     43373      C   python                                     16593MiB |
+-----------------------------------------------------------------------------+
```

## Fair Share / Fair Use

**Note that allocating a GPU makes it unavailable for everyone else.**
**There is currently no technical restriction in place on the GPU nodes, so please behave nicely and cooperatively.**
If you see someone blocking the GPU nodes for long time, then first find out who it is.
You can type `getent passwd USER_NAME` on any cluster node to see the email address (and work phone number if added) of the user.
Send a friendly email to the other user, most likely they blocked the node accidentally.
If you cannot resolve the issue (e.g., the user is not reachable) then please contact hpc-helpdesk@bih-charite.de with your issue.
