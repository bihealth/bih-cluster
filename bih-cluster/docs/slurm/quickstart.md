# Slurm Quickstart

**Create an interactive bash session** (`srun` will run bash in real-time, `--pty` connects its `stdout` and `stderr` to your current session).

```bash
res-login-1:~$ srun --pty bash -i
med0740:~$ echo "Hello World"
Hello World
med0740:~$ exit
res-login-1:~$
```

Note you probably want to **longer running time for your interactive jobs**.
This way, your jobs can run for up to 28 days.
This will make your job be routed automatically into the `long` partition as it is the only one that can fit your job.

```bash
res-login-1:~$ srun --pty --time 28-00 bash -i
med0740:~$
```

**Pro-Tip:** Using Bash aliases for quick access.

```
res-login-1:~$ alias slogin="srun --pty bash -i"
res-login-1:~$ slogin
med0740:~$ exit
res-login-1:~$ cat >>~/.bashrc <<"EOF"
# Useful aliases for logging in via Slurm
alias slogin="srun --pty bash -i"
alias slogin-x11="srun --pty --x11 bash -i"
EOF
```

**Create an interactive R session on the cluster** (assuming conda is active and the environment `my-r` is created, e.g., with `conda create -n my-r r`).

```bash
res-login-1:~$ conda activate my-r
res-login-1:~$ srun --pty R
R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night"
Copyright (C) 2019 The R Foundation for Statistical Computing
[...]
Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.


> Sys.info()["nodename"]
 nodename
"med0740"
> q()
Save workspace image? [y/n/c]:
res-login-1:~$
```

**Create an interactive iPython session** on the cluster (assuming conda is active and the environment `my-python` is created, e.g., with `conda create -n my-python python=3 ipython`).

```bash
res-login-1:~$ conda activate my-python
res-login-1:~$ srun --pty ipython
Python 3.8.2 | packaged by conda-forge | (default, Mar  5 2020, 17:11:00)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.13.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import socket; socket.gethostname()
Out[1]: 'med0740'

In [2]: exit
res-login-1:~$
```

**Allocate 4 cores (default is 1 core), and a total of 4GB of RAM on one node** (alternatively use `--mem-per-cpu` to set RAM per CPU); `sbatch` accepts the same argument.

```bash
res-login-1:~$ srun --cpus-per-task=4 --nodes=1 --mem=4G --pty bash
med0740:~$ export | grep SLURM_CPUS_ON_NODE
4
med0740:~$ your-parallel-script --threads 4
```

**Submit an R script to the cluster in batch mode** (`sbatch` schedules the job for later execution).

```bash
res-login-1:~$ cat >job-script.sh <<"EOF"
#!/bin/bash
echo "Hello, I'm running on $(hostname) and it's $(date)"
EOF
res-login-1:~$ sbatch job-script.sh
Submitted batch job 7

# Some time later:
res-login-1:~$ cat slurm-7.out
Hello, I'm running on med0740 and it's Fri Mar  6 07:36:42 CET 2020
res-login-1:~$
```
