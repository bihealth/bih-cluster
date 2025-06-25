# Migration from old GPFS to new CephFS
!!! warning "Important"
    Access to `/fast` was removed September 30th 2024.
    The following is kept for documenting the process for the curious reader,
    but no action needs to be taken.

## What is going to happen?
Files on the cluster's main storage `/data/gpfs-1` aka. `/fast` will move to a new file system.
That includes users' home directories, work directories, and work-group directories.
Once files have been moved to their new locations, `/fast` will be retired.

Simultaneously we will move towards a more unified naming scheme for project and group folder names.
From now on, all such folders names shall be in [kebab-case](https://en.wikipedia.org/wiki/Letter_case#Kebab_case).
This is Berlin after all.
Group folders will also be renamed, removing the "ag_" prefix.

Detailed communication about the move will be communicated via the cluster mailinglist and the [user forum](https://hpc-talk.cubi.bihealth.org/).
For technical help, please consult the [Data Migration Tips and tricks](migration-faq.md).

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

More information about CephFS can be found [here](storage-locations.md).

## New file locations
Naturally, paths are going to change after files move to their new location.
Due to the increase in storage quality options, there will be some more folders to consider.

### Users
- Home on Tier 1: `/data/cephfs-1/home/users/<user>`
- Work on Tier 1: `/data/cephfs-1/work/groups/<doe>/users/<user>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<doe>/users/<user>`

!!! warning "Important"
    User `work` & `scratch` spaces are now part of the user's group folder.
    This means, groups need to coordinate internally to distribute their allotted quota according to each user's needs.

The implementation is done _via_ symlinks created by default when the user account is moved to its new destination:

- `~/work -> /data/cephfs-1/work/groups/<group>/users/<user>`
- `~/scratch -> /data/cephfs-1/scratch/groups/<group>/users/<user>`

### Groups
- Work on Tier 1: `/data/cephfs-1/work/groups/<group>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<group>`
- Tier 2 storage: `/data/cephfs-2/unmirrored/groups/<group>`
- Mirrored space on Tier 2 is available on request.

### Projects
- Work on Tier 1: `/data/cephfs-1/work/projects/<project>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/projects/<project>`
- Tier 2 storage is available on request.

## Recommended practices
### Data locations
#### Tiers
- Tier 1: Designed for many I/O operations. Store files here which are actively used by your compute jobs.
- Tier 2: Big, cheap storage. Fill with files not in active use.
- Tier 2 mirrored: Extra layer of security. Longer term storage of invaluable data.

#### Folders
- Home: Configuration files, templates, generic scripts, & small documents.
- Work: Conda environments, R packages, data actively processed or analyzed.
- Scratch: Non-persistent storage for temporary or intermediate files, caches, etc.

### Project life cycle
1. Import raw data on Tier 2 for validation (checksums, …)
2. Stage raw data on Tier 1 for QC & processing.
3. Save processing results to Tier 2.
4. Continue analysis on Tier 1.
5. Save analysis results on Tier 2.
6. Reports & publications can remain on Tier 2.
7. After publication (or the end of the project), files on Tier 1 should be deleted.

## Example use cases

Space on Tier 1 is limited.
Your colleagues, other cluster users, and admins will be very grateful if you use it only for files you actively need to perform read/write operations on.
This means main project storage should probably always be on Tier 2 with workflows to stage subsets of data onto Tier 1 for analysis.

These examples are based on our experience of processing diverse NGS datasets.
Your mileage may vary but there is a basic principle that remains true for all projects.

### DNA sequencing (WES, WGS)

Typical Whole Genome Sequencing data of a human sample at 100x coverage requires about 150  GB of storage, Whole Exome Sequencing files occupy between 6 and 30  GB.
These large files require considerable I/O resources for processing, in particular for the mapping step.
A prudent workflow for these kind of analysis would therefore be the following:

1. For one sample in the cohort, subsample its raw data files (`fastqs`) from the Tier 2 location to Tier 1. [`seqtk`](https://github.com/lh3/seqtk) is your friend!
2. Test, improve & check your processing scripts on those smaller files.
3. Once you are happy with the scripts, copy the complete `fastq` files from Tier 2 to Tier 1. Run the your scripts on the whole dataset, and copy the results (`bam` or `cram` files) back to Tier 2.
4. **Remove raw data & bam/cram files from Tier 1**, unless the downstream processing of mapped files (variant calling, structural variants, ...) can be done immediatly.

!!! tip
    Don't forget to use your `scratch` area for transient operations, for example to sort your `bam` file after mapping.
	More information on how to efficiently set up your temporary directory [here](../best-practice/temp-files.md).

### bulk RNA-seq

Analysis of RNA expression datasets are typically a long and iterative process, where the data must remain accessible for a significant period.
However, there is usually no need to keep raw data files and mapping results available once the gene & transcripts counts have been generated.
The count files are much smaller than the raw data or the mapped data, so they can live longer on Tier 1. 

A typical workflow would be:

1. Copy your `fastq` files from Tier 2 to Tier 1.
2. Perform raw data quality control, and store the outcome on Tier 2.
3. Get expression levels, for example using `salmon` or `STAR`, and store the results on Tier 2.
4. Import the expression levels into `R`, using `tximport` and `DESeq2` or `featureCounts` & `edgeR`, for example.
5. Save expression levels (`R` objects) and the output of `salmon`, `STAR`, or any mapper/aligner of your choice to Tier 2.
6. **Remove raw data, bam & count files from Tier 1.**

!!! tip
    If using `STAR`, don't forget to use your `scratch` area for transient operations.
	More information on how to efficiently set up your temporary directory [here](../best-practice/temp-files.md)

### scRNA-seq

The analysis workflow of bulk RNA & single cell dataset is conceptually similar:
Large raw files need to be processed once and only the outcome of the processing (gene counts matrices) are required for downstream analysis.
Therefore, a typical workflow would be:

1. Copy your `fastq` files from Tier 2 to Tier 1.
2. Perform raw data QC, and store the results on Tier 2. 
3. Get the count matrix, e.&nbsp;g. using `Cell Ranger` or `alevin-fry`, perform count matrix QC and store the results on Tier 2.
4. **Remove raw data, bam & count files from Tier 1.**
5. Downstream analysis with `seurat`, `scanpy`, or `Loupe Browser`.

### Machine learning

There is no obvious workflow that covers most used cases for machine learning. 
However,

- Training might be done on scratch where data access is quick and data size not as constrained as on work space. But files will disappear after 14 days.
- Some models can be updated with new data, without needing to keep the whole dataset on Tier 1.

## Data migration process from old `/fast` to CephFS 
1. After being contacted by HPC admins, delegates move project folders to Tier 2. Additional Tier 1 storage is granted on request.
2. User homes and group folders are moved by HPC admins to Tier 1 and 2 as appropriate. This is done on a group-by-group basis.
3. Users move contents of their work directories into the new shared group work space.

Best practices and tools will be provided.
