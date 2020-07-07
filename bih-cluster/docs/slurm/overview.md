# Scheduling Overview

The BIH HPC uses the [Slurm](https://slurm.schedmd.com/overview.html) scheduling system.
This section of the manual attempts to give an overview of what scheduling is and how you can use the Slurm scheduler.
For more detailed information, you will have to refer to the [Slurm website](https://slurm.schedmd.com/overview.html) and the Slurm man pages (e.g., by entering `man sbatch` or `man srun` on the HPC terminal's command line).

For a quick introduction and hands-on examples, please see the manual sections

- Overview, starting with [For the Impatient](/overview/for-the-impatient), and
- First Steps/Tutorial, starting with [Episode 0](/first-steps/episode-0).

## Annotated Contents

- [Background on Scheduling](background.md) -- some background on scheduling and the terminology used
- [Quickstart](quickstart.md) -- the most important Slurm commands, explained, with examples
- [Cheat Sheet](cheat-sheet.md) -- for quick reference
- [Job Scripts](job-scripts.md) -- how to setup job scripts with Slurm
- [Slurm and Snakemake](snakemake.md) -- how to use Snakemake with Slurm
- Introduction to Slurm Commands
    - [`srun`](commands-srun.md) -- running parallel jobs **now**
    - [`sbatch`](commands-sbatch.md) -- submission of batch jobs
    - [`scancel`](commands-scancel.md) -- stop/kill jobs
    - [`sinfo`](commands-sinfo.md) -- display information about the Slurm cluster
    - [`squeue`](commands-squeue.md) -- information about pending and running jbos
    - [`scontrol`](commands-scontrol.md) -- detailed information (and control)
    - [`sacct`](commands-sacct.md) -- access Slurm accounting information (pending, running, and past jobs)
- [Format Strings in Slurm](format-strings.md) -- format strings allow to display extended information about Slurm scheduler objects
- [X11 Forwarding](x11.md) -- X11 forwarding in Slurm (simple; short)
- [Rosetta Stone](rosetta-stone.md) -- lookup table for SGE <-> Slurm
- [Migrating from SGE](migrating.md) -- hints for migrating from SGE to Slurm (:spider_web: deprecated, will be removed)

## A Word on "Elsewhere"

Many facilities run Slurm clusters and have their documentation available on the internet and we will list some that we found useful below.
However, beware that Slurm is a highly configurable and extensible system.
Other sites may have different configurations and plugins enabled than we have (or might even have written custom plugins that are not available at BIH).
In any case, it's always useful to look "über den Tellerrand".

- [Quick Start User Guide](https://slurm.schedmd.com/quickstart.html) - the official guide from the Slurm creators.
- [Slurm `man` Pages](https://slurm.schedmd.com/man_index.html) - web versions of `man <slurm command>`.
- [TU Dresden Slurm Compendium](https://doc.zih.tu-dresden.de/hpc-wiki/bin/view/Compendium/Slurm) - nice documentation from the installation in Dresden.
  Beware that their installation is highly customized, in particular partition selection is automatized for them (and it is not for us).
- [Slurm at CECI](https://support.ceci-hpc.be/doc/_contents/QuickStart/SubmittingJobs/SlurmTutorial.html) - CECI is a HPC consortium from Belgium.
- [Slurm at the Arctic University of Norway](https://hpc-uit.readthedocs.io/en/latest/jobs/examples.html)
- [Slurm at Technical University of Denmark](https://wiki.fysik.dtu.dk/niflheim/SLURM) - if you want to get an insight in how this looks to administrator.
