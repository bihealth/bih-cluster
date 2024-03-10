# Migration from old GPFS to new CephFS
## What is going to happen?
Files on the cluster's main storage `/data/gpfs-1` aka. `/fast` will move to a new file system.
That includes users' home directories, work directories, and work-group directories.
Once files have been moved to their new locations, `/fast` will be retired.

## Why is this happening?
`/fast` is based on a high performance proprietary hardware (DDN) & file system (GPFS).
The company selling it has terminated support which also means buying replacement parts will become increasingly difficult.

## The new storage
There are *two* file systems set up to replace `/fast`, named *Tier 1* and *Tier 2* after their difference in I/O speed:

- **Tier 1** is faster than `/fast` ever was, but it only has about 75 % of its usable capacity.
- **Tier 2** is not as fast, but much larger, almost 3 times the current usable capacity.

The **Hot storage** Tier 1 is reserved for files requiring frequent random access, user homes, and scratch.
Tier 2 (**Warm storage**) should be used for everything else.
Both file systems are based on the open-source, software-defined [Ceph](https://ceph.io/en/) storage platform and differ in the type of drives used.
Tier 1 or Cephfs-1 uses NVME SSDs and is optimized for performance, Tier 2 or Cephfs-2 used traditional hard drives and is optimized for cost.

So these are the three terminologies in use right now:
- Cephfs-1 = Tier 1 = Hot storage = `/data/cephfs-1`
- Cephfs-2 = Tier 2 = Warm storage = `/data/cephfs-2`

There are no more quotas on the number of files.

## New file locations
Naturally, paths are going to change after files move to their new location.
Due to the increase in storage quality options, there will be some more more folders to consider.

### Users
- Home on Tier 1: `/data/cephfs-1/home/users/<user>`
- Work on Tier 1: `/data/cephfs-1/work/groups/<ag-doe>/users/<user>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<ag-doe>/users/<user>`

!!! warning
    User work & scratch spaces are now part of the user's group folder.
    This means, groups should coordinate internally to distribute their allotted quota evenly among users.

The implementation is done _via_ symlinks created by default when the user account is moved to its new destination.

| Symlink named | Points to |
|:--------------|:----------|
| `/data/cephfs-1/home/users/<user>/work`    | `/data/cephfs-1/work/groups/<group>/users/<user>`    |
| `/data/cephfs-1/home/users/<user>/scratch` | `/data/cephfs-1/scratch/groups/<group>/users/<user>` |

Additional symlinks are created from the user's home directory to avoid storing large files (R packages for example) in their home.
The full list of symlinks is:

- HPC web portal cache: `ondemand`
- General purpose caching: `.cache` & `.local`
- Containers: `.singularity` & `.apptainer`
- Programming languages libraries & registry: `R`, `.theano`, `.cargo`, `.cpan`, `.cpanm`, & `.npm`
- Others: `.ncbi`, `.vs`

!!! warning
    Automatic symlink creation will not create a symlink to any conda installation.

### Groups
- Work on Tier 1: `/data/cephfs-1/work/groups/<group>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<group>`
- Mirrored work on Tier 2: `/data/cephfs-2/mirrored/groups/<group>`

!!! note
    Un-mirrored work space on Tier 2 is available on request.

### Projects
- Work on Tier 1: `/data/cephfs-1/work/projects/<project>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/projects/<project>`

!!! note
    Tier 2 work space (mirrored & un-mirrored) is available on request.

## Recommended practices
### Data locations
#### Tiers
- Tier 1: Designed for many I/O operations. Store files here which are actively used by your compute jobs.
- Tier 2: Big, cheap storage. Fill with files not in active use.
- Tier 2 mirrored: Extra layer of security. Longer term storage of invaluable data.

#### Folders
- Home: Persistent storage for configuration files, templates, generic scripts, & small documents.
- Work: Persistent storage for conda environments, R packages, data actively processed or analyzed.
- Scratch: Non-persistent storage for temporary or intermediate files, caches, etc. Automatically deleted after 14 days.

### Project life cycle
1. Import the raw data on Tier 2 for validation (checksums, …)
2. Stage raw data on Tier 1 for QC & processing.
3. Save processing results to Tier 2 and validate copies.
4. Continue analysis on Tier 1.
5. Save analysis results on Tier 2 and validate copies.
6. Reports & publications can remain on Tier 2.
7. After publication (or the end of the project), files on Tier 1 can be deleted.

### Example use cases

Space on Tier 1 is limited.
Your colleagues, other cluster users, and admins will be very grateful if you use it only for files you actively need to perform read/write operations on.
This means main project storage should probably always be on Tier 2 with workflows to stage subsets of data onto Tier 1 for analysis.

These examples are based on our experience of processing diverse NGS datasets.
Your mileage may vary but there is a basic principle that remains true for all projects.

#### DNA sequencing (WES, WGS)

Typical Whole Genome Sequencing data of a human sample at 100x coverage requires about 150  GB of storage, Whole Exome Sequencing between 6 and 30  GB.
These large files require considerable I/O resources for processing, in particular for the mapping step.
A prudent workflow for these kind of analysis would therefore be the following:

1. For one sample in the cohort, subsample its raw data files (`fastqs`) from the Tier 2 location to Tier 1. [`seqtk`](https://github.com/lh3/seqtk) is your friend!
2. Test, improve & check your processing scripts on those smaller files.
3. Once you are happy with the scripts, copy the complete `fastq` files from Tier 2 to Tier 1. Run the your scripts on the whole dataset, and copy the results (`bam` or `cram` files) back to Tier 2.
4. **Remove raw data & bam/cram files from Tier 1**, unless the downstream processing of mapped files (variant calling, structural variants, ...) can be done immediatly.

!!! tip
    Don't forget to use your `scratch` area for transient operations, for example to sort your `bam` file after mapping.
	More information on how to efficiently set up your temporary directory [here](../best-practice/temp-files.md).

#### bulk RNA-seq

Analysis of RNA expression datasets are typically a long and iterative process, where the data must remain accessible for a significant period.
However, there is usually no need to keep raw data files and mapping results available once the gene & transcripts counts have been generated.
The count files are much smaller than the raw data or the mapped data, so they can live longer on Tier 1. 

A typical workflow would be:

1. Copy your `fastq` files from Tier 2 to Tier 1.
2. Get expression levels, for example using `salmon` or `STAR`, and store the results on Tier 1.
3. Import the expression levels into `R`, using `tximport` and `DESeq2` or `featureCounts` & `edgeR`, for example.
4. Save expression levels (`R` objects) and the output of `salmon`, `STAR`, or any mapper/aligner of your choice to Tier 2.
5. **Remove raw data, bam & count files from Tier 1.**

!!! tip
    If using `STAR`, don't forget to use your `scratch` area for transient operations.
	More information on how to efficiently set up your temporary directory [here](../best-practice/temp-files.md)

#### scRNA-seq

The analysis workflow of bulk RNA & single cell dataset is conceptually similar:
Large raw files need to be processed once and only the outcome of the processing (gene counts matrices) are required for downstream analysis.
Therefore, a typical workflow would be:

1. Copy your `fastq` files from Tier 2 to Tier 1.
2. Perform raw data QC (for example with `fastqc`).
3. Get the count matrix, e.  g. using `Cell Ranger` or `alevin-fry`, perform count matrix QC and store the results on Tier 1.
4. **Remove raw data, bam & count files from Tier 1.**
5. Downstream analysis with `seurat`, `scanpy`, or `Loupe Browser`.

#### Machine learning

## Data migration process from old `/fast` to CephFS 
1. Administrative preparations  
    1. HPC-Access registration (PIs will receive in invite mail)
    2. PIs need to assign a delegate.
    3. PI/Delegate needs to add group and the projects once.
    4. New Tier 1 & 2 resources will be allocated.
2. Users & group home directories are moved by HPC admin (big bang).
3. All directories on `/fast` set to read-only, that is:
    - `/fast/home/users/<user>` & `/fast/work/users/<user>`
    - `/fast/home/groups/<group>` & `/fast/work/groups/<group>`
4. Work data migration is done by the users (Tier 2 is primary target, Tier 1 staging when needed).  

Best practice and/or tools will be provided.

!!! note
    The users' `work` space will be moved to the group's `work` space.
