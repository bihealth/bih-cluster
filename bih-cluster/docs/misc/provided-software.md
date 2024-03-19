# Administration-Provided Software

Some software is provided by HPC Administration based on the criteria that it is:

- system-near or system-level,
- very commonly used.

Currently, this includes:

- GCC v7.2.0
- CMake v3.11.0
- LLVM v6.0.0
- OpenMPI v4.0.3

On the GPU node, this also includes a recent NVIDIA CUDA version.

To see which software is available, use `module avail` on a compute node (this will not work on login nodes):

```terminal
$ module avail
--------------------- /opt/local/modules ---------------------
cmake/3.11.0-0  llvm/6.0.0-0
gcc/7.2.0-0     openmpi/4.0.3-0
```

To load software, use `module load`.
This will adjust the environment variables accordingly, in particular update `PATH` such that the executable are available.

```terminal
$ which gcc
/bin/gcc
$ module load gcc/7.2.0-0
$ which gcc
/opt/local/gcc-7.2.0-0/bin/gcc
```

!!! hint "Problems with executing `module`?"

    See the [corresponding FAQ entry](../../help/faq/#environment-modules-dont-work-and-i-get-module-command-not-found) in the case that you get a `-bash: module: command not found` when calling `module`.
