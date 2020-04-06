# Temporary Files

Often, it is necessary to use temporary files, i.e., write something out in the  middle of your program, read it in again later, and then discard these files.
For example, `samtools sort` has to write out chunks of sorted read alignments for allowing to sort files larger than main memory.

## Environment Variable `TMPDIR`

Traditionally, in Unix, the environment variables `TMPDIR` is used for storing the location of the temporary directory.
When undefined, usually `/tmp` is used.

## Temporary Directories on the BIH Cluster

Generally, there are two locations where you could put temporary files:

- `/fast/users/$USER/scratch/tmp` -- inside your scratch folder on the fast GPFS file system; this location is available from all cluster nodes
- `/tmp` -- on the local node's temporary folder; this location is only available on the node itself

### Best Practice:  Use `/fast/users/$USER/scratch/tmp`

!!! warning "Use GPFS-based TMPDIR"
    Generally setup your environment to use `/fast/users/$USER/scratch/tmp` as filling the local disk of a node with forgotten files can cause a lot of problems.

Ideally, you append the following to your `~/.bashrc` to use `/fast/users/$USER/scratch/tmp` as the temporary directory.
This will also create the directory if it does not exist.

```bash
export TMPDIR=/fast/users/$USER/scratch/tmp
mkdir -p $TMPDIR
```

**Prepending this to your job scripts is also recommended as it will ensure that the temporary directory exists.**

## `TMPDIR` and the SGE scheduler

In the current configuration, jobs submitted through `qsub` will have `TMPDIR` set to a sub directory in `/tmp` on the local machine that is automatically cleaned up after a job finishes.
In the older nodes, the local disk is a relatively slow spinning disk, in the newer nodes, the local disk is a relatively fast SSD.
Further, the local disk is independent from the GPFS file system, so I/O volume to it does not affect the network or any other job on other nodes.

!!! imporant "SGE and TMPDIR"
    Note that the SGE scheduler will override any value given to `TMPDIR` from the outside. If you want to override `TMPDIR` you have to do this within your job script.

For smaller I/O volumes, it is probably sensible to just use `TMPDIR` from the scheduler.
For larger ones, you probably want to fall back to `/fast/users/$USER/scratch/tmp`.
However, this is then not auto-cleaned any more.
Here is how you get this behaviour.

### Use Bash Trapcs

You can use the following code at the top of your job script to set `TMPDIR` to the location in your home directory and get the directory automatically cleaned when the job is done (regardless of successful or erroneous completion):

```bash
# First, point TMPDIR to the scratch in your home as mktemp will use thi
export TMPDIR=$HOME/scratch/tmp
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap
trap "rm -rf $TMPDIR" EXIT
```
