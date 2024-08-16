# How-To: Install Custom Scientific Software

This page gives an end-to-end example how to build and install [Gromacs](http://www.gromacs.org/) as an example for managing complex scientific software installs in user land.
You don't have to learn or understand the specifics of Gromacs.
We use it as an example as there are some actual users on the BIH cluster.
However, installing it is out of scope of BIH HPC administration.

Gromacs is a good example as it is a sufficiently complex piece of software.
Quite some configuration is done on the command line and there is no current software package of it in the common RPM repositories.
However, it is quite well-documented and easy to install for scientific software so there is a lot to be learned.

## Related Documents

- [How-To: Build and Run OpenMPI Programs](openmpi.md)

## Steps for Installing Scientific Software

We will perform the following step:

1. Download and extract the source of the software
2. Configure the software (i.e., create the actual build system `Makefile`s)
3. Compile the software
4. Install the software
5. Create [environment module](http://modules.sourceforge.net/) files so the software is easy to use

Many scientific software packages will have more dependencies.
If the dependencies are available as CentOS Core or EPEL packages (such as zlib), HPC IT administration can install them.
However, otherwise you will have to install them on their own.

!!! warning

    Do not perform the compilation on the login nodes but go to a compute node instead.

## Downloading and Extracting Software

This is best done in your `scratch` directory as we don't have to keep these files around for long.
Note that the files in your `scratch` directory will automatically be removed after 2 weeks.
You can also use your `work` directory here.

```bash
hpc-login-1:~$ srun --pty bash -i
med0127:~$ mkdir $HOME/scratch/gromacs-install
med0127:~$ cd $HOME/scratch/gromacs-install
med0127:~$ wget http://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz
med0127:~$ tar xf gromacs-2018.3.tar.gz
med0127:~$ ls gromacs-2018.3
admin    cmake           COPYING          CTestConfig.cmake  INSTALL  scripts  src
AUTHORS  CMakeLists.txt  CPackInit.cmake  docs               README   share    tests
```

So far so good!

## Perform the Configure Step

This is the most critical step.
Most scientific C/C++ software has a build step and allows for, e.g., disabling and enabling features or setting installation paths.
Here, you can configure the software depending on your needs and environment.
Also, it is the easiest step to mess up.

Gromac's documentation is actually quite good but the author had problems to follow it to the letter.
Gromacs recommends to create an MPI and a non-MPI build but the precise way did not work.
This installation creates two flavours for Gromacs 2018.3, but in a different way than the Gromacs documentation proposes.

First, here is how to configure the non-MPI flavour
Gromacs wants a modern compiler, so we load `gcc`.
We will need to note down the precise version we used so later we can load it for running Gromacs with the appropriate libraries.
We will install gromacs into `$HOME/work/software`, which is appropriate for user-installed software, but it could also go into a group or project directory.
Note that we install the software into your work directory as software installations are quite large and might go above your home quota.
Also, software installations are usually not precious enough to waste resources on snapshots and backups.
Also that we force Gromacs to use `AVX_256` for SIMD support (Intel *sandy bridge* architecture) to not get unsupported CPU instruction errors.

```bash
med0127:~$ module load gcc/7.2.0-0 cmake/3.11.0-0
med0127:~$ module list
Currently Loaded Modulefiles:
  1) gcc/7.2.0-0      2) cmake/3.11.0-0
med0127:~$ mkdir gromacs-2018.3-build-nompi
med0127:~$ cd gromacs-2018.3-build-nompi
med0127:~$ cmake ../gromacs-2018.3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_MPI=OFF \
    -DGMX_SIMD=AVX_256 \
    -DCMAKE_INSTALL_PREFIX=$HOME/work/software/gromacs/2018.3
```

Second, here is how to configure the MPI flavour.
Note that we are also enabling the `openmpi` module.
We will also need the precise version here so we can later load the correct libraries.
Note that we install the software into the directory `gromacs-mpi` but switch off shared library building as recommended by the Gromacs documentation.

```bash
med0127:~$ module load openmpi/3.1.0-0
med0127:~$ module list
Currently Loaded Modulefiles:
  1) gcc/7.2.0-0       2) cmake/3.11.0-0    3) openmpi/4.0.3-0
med0127:~$ mkdir gromacs-2018.3-build-mpi
med0127:~$ cd gromacs-2018.3-build-mpi
med0127:~$ cmake ../gromacs-2018.3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_MPI=ON \
    -DGMX_SIMD=AVX_256 \
    -DCMAKE_INSTALL_PREFIX=$HOME/work/software/gromacs-mpi/2018.3 \
    -DCMAKE_C_COMPILER=$(which mpicc) \
    -DCMAKE_CXX_COMPILER=$(which mpicxx) \
    -DBUILD_SHARED_LIBS=off
```

## Perform the Build and Install Steps

This is simple, using `-j 32` allows us to build with 32 threads.
If something goes wrong: meh, the "joys" of compilling C software.

!!! note "Getting Support for Building Software"

    BIH HPC IT cannot provide support for compiling scientific software.
    Please contact the appropriate mailing lists or forums for your scientific software.
    You should only contact the BIH HPC IT helpdesk only if you are sure that the problem is with the BIH HPC cluster.
    You should try to resolve the issue on your own and with the developers of the software that you are trying to build/use.

For the no-MPI version:

```bash
med0127:~$ cd ../cd gromacs-2018.3-build-nompi
med0127:~$ make -j 32
[...]
med0127:~$ make install
```

For the MPI version:

```bash
med0127:~$ cd ../cd gromacs-2018.3-build-mpi
med0127:~$ make -j 32
[...]
med0127:~$ make install
```

## Create Environment Modules Files

For Gromacs 2018.3, the following is appropriate.
You should be able to use this as a template for your environment module files:

```bash
med0127:~$ mkdir -p $HOME/local/modules/gromacs
med0127:~$ cat >$HOME/local/modules/gromacs/2018.3 <<"EOF"
#%Module
proc ModulesHelp { } {
    puts stderr {
      Gromacs molecular simulation toolkit (non-MPI version)

      - http://www.gromacs.org
    }
}

module-whatis {Gromacs molecular simulation toolkit (non-MPI)}

set root /data/cephfs-1/home/users/YOURUSER/work/software/gromacs-mpi/2018.3

prereq gcc/7.2.0-0

conflict gromacs
conflict gromacs-mpi

prepend-path    LD_LIBRARY_PATH         $root/lib64
prepend-path    LIBRARY_PATH            $root/lib64
prepend-path    MANPATH                 $root/share/man
prepend-path    PATH                    $root/bin
setenv          GMXRC                   $root/bin/GMXRC
EOF
```

```bash
med0127:~$ mkdir -p $HOME/local/modules/gromacs-mpi
med0127:~$ cat >$HOME/local/modules/gromacs-mpi/2018.3 <<"EOF"
#%Module
proc ModulesHelp { } {
    puts stderr {
      Gromacs molecular simulation toolkit (MPI version)

      - http://www.gromacs.org
    }
}

module-whatis {Gromacs molecular simulation toolkit (MPI)}

set root /data/cephfs-1/home/users/YOURUSER/work/software/gromacs-mpi/2018.3

prereq openmpi/4.0.3-0
prereq gcc/7.2.0-0

conflict gromacs
conflict gromacs-mpi

prepend-path    LD_LIBRARY_PATH         $root/lib64
prepend-path    LIBRARY_PATH            $root/lib64
prepend-path    MANPATH                 $root/share/man
prepend-path    PATH                    $root/bin
setenv          GMXRC                   $root/bin/GMXRC
EOF
```

With the next command, make your local modules files path known to the environemtn modules system.

```bash
med0127:~$ module use $HOME/local/modules
```

You can verify the result:

```bash
med0127:~$ module avail

------------------ /data/cephfs-1/home/users/YOURUSER/local/modules ------------------
gromacs/2018.3     gromacs-mpi/2018.3

-------------------- /usr/share/Modules/modulefiles --------------------
dot         module-info null
module-git  modules     use.own

-------------------------- /opt/local/modules --------------------------
cmake/3.11.0-0  llvm/6.0.0-0    openmpi/3.1.0-0
gcc/7.2.0-0     matlab/r2016b-0 openmpi/4.0.3-0
```

### Interlude: Convenient `module use`

You can add this to your `~/.bashrc` file to always execute the `module use` after login.
Note that `module` is not available on the login or transfer nodes, the following should work fine:

```bash
med0127:~$ cat >>~/.bashrc <<"EOF"
case "${HOSTNAME}" in
  login-*|transfer-*)
    ;;
  *)
    module use $HOME/local/modules
    ;;
esac
EOF
```

Note that the paths chosen above are sensible but arbitrary.
You can install any software anywhere you have permission to -- somewhere in your user and group home, maybe a project home makes most sense on the BIH HPC, no root permissions required.
You can also place the module files anywhere, as long as the `module use` line is appropriate.

As a best practice, you could use the following location:

- User-specific installation:
    - `$HOME/work/software` as a root to install software to
    - `$HOME/work/software/$PKG/$VERSION` for installing a given software package in a given version
    - `$HOME/work/software/modules` as the root for modules to install
    - `$HOME/work/software/$PKG/$VERSION` for the module file to load the software in a given version
    - `$HOME/work/software/modules.sh` as a Bash script to contain the line `module use $HOME/work/software/modules`
- Group/project specific installation for a shared setup.
  Don't forget to give the group and yourself read permission only so you don't accidentally damage files after instalation (`chmod ug=rX,o= $GROUP/work/software`, the upper case `X` is essential to only set `+x` on directories and not files):
    - `$GROUP/work/software` as a root to install software to
    - `$GROUP/work/software/$PKG/$VERSION` for installing a given software package in a given version
    - `$GROUP/work/software/modules` as the root for modules to install
    - `$GROUP/work/software/$PKG/$VERSION` for the module file to load the software in a given version
    - `$GROUP/work/software/modules.sh` as a Bash script to contain the `case` Bash snippet from above but with `module use $GROUP/work/software/modules`
    - This setup allows multiple users to provide software installations and share it with others.

### Going on with Gromacs

Every time you want to use Gromacs, you can now do

```bash
med0127:~$ module load gcc/7.2.0-0 gromacs/2018.3
```

or, if you want to have the MPI version:

```bash
med0127:~$ module load gcc/7.2.0-0 openmpi/4.0.3-0 gromacs-mpi/2018.3
```

## Launching Gromacs

Something along the lines of the following job script should be appropriate.
See [How-To: Build Run OpenMPI Programs](openmpi.md) for more information.

```bash
#!/bin/bash

# Example job script for (single-threaded) MPI programs.

# Generic arguments

# Job name
#SBATCH --job-name gromacs
# Maximal running time of 10 min
#SBATCH --time 00:10:00
# Allocate 1GB of memory per CPU
#SBATCH --mem 1G
# Write logs to directory "slurm_log/<name>-<job id>.log" (dir must exist)
#SBATCH --output slurm_log/%x-%J.log

# MPI-specific parameters

# Launch on 8 nodes (== 8 tasks)
#SBATCH --ntasks 8
# Allocate 4 CPUs per task (== per node)
#SBATCH --cpus-per-task 4

# Load the OpenMPI and GCC environment module to get the runtime environment.
module load gcc/4.7.0-0
module load openmpi/4.0.3-0

# Make custom environment modules known. Alternative, you can "module use"
# them in the session you use for submitting the job.
module use $HOME/local/modules
module load gromacs-mpi/2018.3

# Launch the program on 8 nodes and tell Gromacs to use 4 threads for each
# invocation.
export OMP_NUM_THREADS=4
mpirun -n 8 gmx_mpi mdrun -deffnm npt_1000
```

```bash
med0127:~$ mkdir slurm_log
med0127:~$ sbatch job_script.sh
Submitted batch job 3229
```
