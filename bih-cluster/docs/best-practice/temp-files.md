# Temporary Files

!!! info "Temporary Files and Slurm"

    See [Slurm: Temporary Files](../slurm/temporary-files.md) for information how Slurm controls access to local temporary storage.

Often, it is necessary to use temporary files, i.e., write something out in the  middle of your program, read it in again later, and then discard these files.
For example, `samtools sort` has to write out chunks of sorted read alignments for allowing to sort files larger than main memory.

## Environment Variable `TMPDIR`

Traditionally, in Unix, the environment variables `TMPDIR` is used for storing the location of the temporary directory.
When undefined, usually `/tmp` is used.

## Temporary Directories on the BIH Cluster

Generally, there are two locations where you could put temporary files:

- `/data/cephfs-1/home/users/$USER/scratch/tmp` -- inside your scratch folder on the CephFS file system; this location is available from all cluster nodes
- `/tmp` -- on the local node's temporary folder; this location is only available on the node itself.
  The slurm scheduler uses Linux namespaces such that every **job** gets its private `/tmp` even when run on the same node.

### Best Practice:  Use `scratch/tmp`

!!! warning "Use CephFS-based TMPDIR"

    Generally setup your environment to use `/data/cephfs-1/home/users/$USER/scratch/tmp` as filling the local disk of a node with forgotten files can cause a lot of problems.

Ideally, you append the following to your `~/.bashrc` to use `/data/cephfs-1/home/users/$USER/scratch/tmp` as the temporary directory.
This will also create the directory if it does not exist.
Further, it will create one directory per host name which prevents too many entries in the temporary directory.

```bash
export TMPDIR=$HOME/scratch/tmp/$(hostname)
mkdir -p $TMPDIR
```

**Prepending this to your job scripts is also recommended as it will ensure that the temporary directory exists.**

## `TMPDIR` and the scheduler

In the older nodes, the local disk is a relatively slow spinning disk, in the newer nodes, the local disk is a relatively fast SSD.
Further, the local disk is independent from the CephFS file system, so I/O volume to it does not affect the network or any other job on other nodes.
Please note that by default, Slurm will not change your environment variables.
This includes the environment variable `TMPDIR`.

Slurm will automatically update temporary files in a job's `/tmp` on the local file system when the job terminates.
To automatically clean up temporary directories on the shared file system, use the following tip.

### Use Bash Traps

You can use the following code at the top of your job script to set `TMPDIR` to the location in your home directory and get the directory automatically cleaned when the job is done (regardless of successful or erroneous completion):

```bash
# First, point TMPDIR to the scratch in your home as mktemp will use thi
export TMPDIR=$HOME/scratch/tmp
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap
trap "rm -rf $TMPDIR" EXIT
```
