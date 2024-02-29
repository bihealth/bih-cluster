## What is going to happen?
Files on the cluster's main storage (`/data/gpfs-1` aka. `/fast`) will move to a new filesystem.
That includes users' home directories, work directories, and workgroup directories.
Once files have been moved to their new locations, `/fast` will be retired.

## Why is this happening?
`/fast` is based on a high performance proprietary hardware (DDN) & filesystem (GPFS).
The company selling it has terminated support which also means buying replacement parts will become increasingly difficult.

## The new storage
There are *two* filesystems set up to replace `/fast`, named *Tier 1* and *Tier 2* after their difference in I/O speed:

- **Tier 1** is faster than `/fast` ever was, but it only has about 75  % of its usable capacity.
- **Tier 2** is not as fast, but much larger, almost 3 times the current usable capacity.

The **Hot storage** Tier 1 is reserved for large files, requiring frequent random access.
Tier 2 (**Warm storage**) should be used for everything else.
Both filesystems are based on the open-source, software-defined [Ceph](https://ceph.io/en/) storage platform and differ in the type of drives used.
Tier 1 or Cephfs-1 uses NVME SSDs and is optimized for performance, Tier 2 or Cephfs-2 used traditional hard drives and is optimised for cost.

So these are the three terminologies in use right now:
- Cephfs-1 = Tier 1 = Hot storage
- Cephfs-2 = Tier 2 = Warm storage

### Snapshots and Mirroring
Snapshots are incremental copies of the state of the data at a particular point in time. 
They provide safety against various "Oups, did I just delete that?" scenarios, meaning they can be used to recover lost or damaged files.

Depending on the location and Tier, Cephfs utilises snapshots in differ differently.
Some parts of Tier 1 and Tier 2 snapshots are also mirrored into a separate fire compartment within the datacenter to provide an additional layer of security.

| Tier | Location                 | Path                         | Retention policy                | Mirrored |
|-----:|:-------------------------|:-----------------------------|:--------------------------------|:---------|
|    1 | User homes               | `/data/cephfs-1/home/users/` | Hourly for 48 h, daily for 14 d | yes      |
|    1 | Group/project work       | `/data/cephfs-1/work/`       | Four times a day, daily for 5 d | yes      |
|    1 | Group/project scratch    | `/data/cephfs-1/scratch/`    | Daily for 3 d                   | no       |
|    2 | Group/project mirrored   | `/data/cephfs-2/mirrored/`   | Daily for 30 d, weekly for 16 w | yes      |
|    2 | Group/project unmirrored | `/data/cephfs-2/unmirrored/` | Daily for 30 d, weekly for 16 w | no       |

User access to the snapshots is documented here: https://hpc-docs.cubi.bihealth.org/storage/accessing-snapshots

### Quotas

| Tier | Function | Path | Default Quota |
|-----:|:---------|:-----|--------------:|
|    1 | User home          | `/data/cephfs-1/home/users/<user>`              | 1 GB        |
|    1 | Group work         | `/data/cephfs-1/work/groups/<group>`           | 1 TB        |
|    1 | Group scratch      | `/data/cephfs-1/scratch/groups/<group>`        | 10 TB       |
|    1 | Projects work      | `/data/cephfs-1/work/projects/<project>`       | individual  |
|    1 | Projects scratch   | `/data/cephfs-1/scratch/projects/<project>`    | individual  |
|    2 | Group mirrored     | `/data/cephfs-2/mirrored/groups/<group>`       | 4 TB        |
|    2 | Group unmirrored   | `/data/cephfs-2/unmirrored/groups/<group>`     | On request  |
|    2 | Project mirrored   | `/data/cephfs-2/mirrored/projects/<project>`   | On request  |
|    2 | Project unmirrored | `/data/cephfs-2/unmirrored/projects/<project>` | On request  |

There are no quotas on the number of files.

## New file locations
Naturally, paths are going to change after files move to their new location.
Due to the increase in storage quality options, there will be some more more folders to consider.

### Users
- Home on Tier 1: `/data/cephfs-1/home/users/<user>`
- Work on Tier 1: `/data/cephfs-1/work/groups/<ag-doe>/users/<user>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<ag-doe>/users/<user>`

[!IMPORTANT]
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
- Programming languanges libraries & registry: `R`, `.theano`, `.cargo`, `.cpan`, `.cpanm`, & `.npm`
- Others: `.ncbi`, `.vs`

[!IMPORTANT]
Automatic symlink creation will **not create a symlink to any conda installation**.

### Groups
- Work on Tier 1: `/data/cephfs-1/work/groups/<group>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/groups/<group>`
- Mirrored work on Tier 2: `/data/cephfs-2/mirrored/groups/<group>`

[!NOTE]
Un-mirrored work space on Tier 2 is available on request.

### Projects
- Work on Tier 1: `/data/cephfs-1/work/projects/<project>`
- Scratch on Tier 1: `/data/cephfs-1/scratch/projects/<project>`

[!NOTE]
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
#### DNA sequencing (WES, WGS)

#### bulk RNA-seq

#### scRNA-seq

#### Machine learning

## Data migration process from old `fast/` to CephFS 
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

[!NOTE]
The users' `work` space will be moved to the group's `work` space.

## Technical details on the new infrastructure
### Tier 1
- Fast & expensive (flash drive based), mounted on `/data/cephfs-1`
- Currently 12 Nodes with 10 × 14 TB NVME/SSD each installed
    - 1.68 PB raw storage
    - 1.45 PB erasure coded (EC 8:2)
    - 1.23 PB usable (85 %, ceph performance limit)
- For typical CUBI use case 3 to 5 times faster I/O then the old DDN
- Two more nodes in purchasing process
- Example of flexible extension:
    - Chunk size: 45 kE for one node with 150 TB, i.e. ca. 300 E/TB

### Tier 2
- Slower but more affordable (spinning HDDs), mounted on `/data/cephfs-2`
- Currently ten nodes with 52 HDDs slots plus SSD cache installed, per node ca. 40 HDDs with 16 to 18 TB filled, i.e.
    - 6.6 PB raw
    - 5.3 PB erasure coded (EC 8:2)
    - 4.5 PB usable (85 %; Ceph performance limit)
- Nine more nodes in purchasing process with 5+ PB
- Very Flexible Extension possible:
    - ca. 50 Euro per TB, 100 Euro mirrored, starting at small chunk sizes
    
### Tier 2 mirror
Similar hardware and size duplicate (another 10 nodes, 6+ PB) in separate fire compartment
