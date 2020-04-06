# Custom Environment Modules

This document contains a few tips for helping you using environment modules more effectively.
As the general online documentation is lacking a bit, we also give the most popular commands here.

## How does it Work?

Environment modules are descriptions of software packages.
The `module` command is provided which allows the manipulation of environment variables such as `PATH`, `MANPATH`, etc., such that programs are available without passing the full path.
Environment modules also allow specifying dependencies between packages and conflicting packages (e.g., when the same binary is available in two packages).
Further, environment variables allow the parallel installation of different software versions in parallel and then using software "a la carte" in your projects.

## Popular Commands

### Querying

List currently loaded modules:

```terminal
$ module list
```

Show all available modules

```terminal
$ module avail
```

### Loading/Unloading Modules

Load one module, make sure to use a specific version to avoid ambiguities.

```terminal
$ module load Jannovar/0.16-Java-1.7.0_80
```

Unload one module

```terminal
$ module unload Jannovar
```

Unload all modules

```terminal
$ module purge
```

### Getting Help

Get help for environment modules

```terminal
$ module help
```

Get help for a particular environment module

```terminal
$ module help Jannovar/0.16-Java-1.7.0_80
```

### Using your own Module Files

You can also create your own environment modules.
Simply create a directory with module files and then use `module use` for using the modules from the directory tree.

```terminal
$ module use path/to/modules
```

## FAQ: Why `-bash: module: command not found`?

On the login nodes, the `module` command is not installed.
You should not run any computations there, so why would you need environment modules there? ;)

```terminal
meg-login2$ module
-bash: module: command not found
```

Simply use `qrsh` to get to one of the compute nodes.

## Auto-loading a set of Modules

You will certainly finding yourself using a set of programs regularly without it being part of the core cluster installation, e.g., SAMtools, or Python 3.
Just putting the appropriate `module load` lines in your `~/.bashrc` will generate warnings when logging into the login node.
It is thus recommended to use the following snippet for loading modules automatically on logging into a compute node:

```bash
case "${HOSTNAME}" in
    med-login*)
        ;;
    *)
        # load Python3 environment module
        module load Python/3.4.3-foss-2015a

        # Define path for temporary directories, don't forget to cleanup!
        # Also, this will only work after /fast is available.
        export TMPDIR=/fast/users/$USER/scratch/tmp

        # Make DRMAA library from central installation available
        export DRMAA_LIBRARY_PATH=/opt/sge/lib/lx-amd64/libdrmaa.so
        ;;
esac
```
