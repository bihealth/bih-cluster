# Slurm Overview

The "Slurm Scheduler" section assumes that you are already familiar with the HPC clusters in general, the BIH HPC in particular.
For example, you have read the sections

- Overview, starting with [For the Impatient](../overview/for-the-impatient), and
- First Steps/Tutorial, starting with [Episode 0](../first-steps/episode-0).

Thus, the intended audience has already used HPC clusters (or the BIH HPC) a bit and is familiar with the basics of Slurm.
The subsections will walk you through the different aspects of using Slurm.
As we will not be able to capture all aspects, we will link out to other resources on the web.

## Annotated Contents

- [Quickstart](quickstart.md) -- the most important Slurm commands, explained, with examples
- [Cheat Sheet](cheat-sheet.md) -- for quick reference
- [Job Scripts](job-scripts.md) -- how to setup job scripts with Slurm
- [and Snakemake](snakemake.md) -- how to use Snakemake with Slurm
- [X11 Forwarding](x11.md) -- X11 forwarding in Slurm (simple; short)
- [Rosetta Stone](rosetta-stone.md) -- lookup table for SGE <-> Slurm
- [Migrating from SGE](migrating.md) -- hints for migrating from SGE to Slurm

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