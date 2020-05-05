# Slurm Quickstart

**Create an interactive bash session** (`srun` will run bash in real-time, `--pty` makes `stdout` and `stderr` to your current session).

```bash
med-login1:~$ srun --pty bash -i
med0740:~$ echo "Hello World"
Hello World
med0740:~$ exit
med-login1:~$
```

Note you probably want to **specify the long partition and a longer running time for your interactive jobs**.
This way, your jobs can run for up to 28 days.

```bash
med-login1:~$ srun --pty --partition long --time 28-00 bash -i
med0740:~$
```

**Pro-Tip:** Using Bash aliases for quick access.

```
med-login1:~$ alias slogin="srun --pty bash -i"
med-login1:~$ slogin
med0740:~$ exit
med-login1:~$ cat >>~/.bashrc <<"EOF"
# Useful aliases for logging in via Slurm
alias slogin="srun --pty bash -i"
alias slogin-x11="srun --pty --x11 bash -i"
EOF
```

**Create an interactive R session on the cluster** (assuming conda is active and the environment `my-r` is created, e.g., with `conda create -n my-r r`).

```bash
med-login1:~$ conda activate my-r
med-login1:~$ srun --pty R
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
med-login1:~$
```

**Create an interactive iPython session** on the cluster (assuming conda is active and the environment `my-python` is created, e.g., with `conda create -n my-python python=3 ipython`).

```bash
med-login1:~$ conda activate my-python
med-login1:~$ srun --pty ipython
Python 3.8.2 | packaged by conda-forge | (default, Mar  5 2020, 17:11:00)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.13.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import socket; socket.gethostname()
Out[1]: 'med0740'

In [2]: exit
med-login1:~$
```

**Allocate 4 cores (default is 1 core), and a total of 4GB of RAM on one node** (alternatively use `--mem-per-cpu` to set RAM per CPU); `sbatch` accepts the same argument.

```bash
med-login1:~$ srun --cpus-per-task=4 --nodes=1 --mem=4G --pty bash
med0740:~$ export | grep SLURM_CPUS_ON_NODE
4
med0740:~$ your-parallel-script 4
```

**Submit an R script to the cluster in batch mode** (`sbatch` schedules the job for later execution).

```bash
med-login1:~$ cat >job-script.sh <<"EOF"
#!/bin/bash
echo "Hello, I'm running on $(hostname) and it's $(date)"
EOF
med-login1:~$ sbatch job-script.sh
Submitted batch job 7
med-login1:~$ cat slurm-7.out
Hello, I'm running on med0740 and it's Fri Mar  6 07:36:42 CET 2020
med-login1:~$
```
