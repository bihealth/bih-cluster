# Software Installation with Conda

## Conda

For the management of the bioinformatics software on the BIH cluster we are using conda.
Conda is a package management system that is based on channels, and one of those
channels provides a huge selection of bioinformatics software.

Conda is written in Python and is based on recipes, such that everybody can
write recipes for missing software (if there is any). In general the packages
are pre-compiled and conda just downloads the binaries from the conda servers.

You are in charge of managing your own software stack, but conda makes it easy
to do so. We will provide you with a description on how to install conda and how
to use it. Of course there are many online resources that you can also use.
Please find a list at the end of the document.

Also note that some system-level software is managed through environment modules.
See [System-near Software Provided by HPC Administration](#system-near-software-provided-by-hpc-administration) below.

## Premise

When you logged into the cluster, please make sure that you also executed
`qrsh` to log into a computation node and perform the software installation
there.

## Installing conda

```terminal
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/work/miniconda
```

This will install conda to `$HOME/work/miniconda`.
This path can be changed to your liking.
Please note that the `$HOME` folder has limited space (an exception is the subfolder `$HOME/work` which has no space limit).

NB: `$HOME/scratch` is not appropriate as files placed there will be removed automatically after 4 weeks.

To make it available upon login, extend and export the `$PATH` variable with the
installation path + `/bin` and add it to your `$HOME/.bashrc`:

```bash
case "${HOSTNAME}" in
    med-login*)
        ;;
    *)
        export PATH=$HOME/work/miniconda/bin:$PATH
        ;;
esac
```

The above code makes sure that you don't have conda available on the login nodes,
where you are not allowed to start any computations.

To make bioinformatics software available, we have to add the `bioconda` and
some other channels to the conda configuration:

```terminal
$ conda config --add channels bioconda
$ conda config --add channels default
$ conda config --add channels conda-forge
```

You can also add channels to your liking.

## Installing software with conda

Installing packages with conda is straight forward:

```terminal
$ conda install <package>
```

This will install a package into the conda root environment. We will explain
environments in detail in the next section.

To search for a package, e.g. to find the correct name in conda or if it exists
at all, issue the command:

```terminal
$ conda search <string>
```

To choose a specific version (conda will install the latest version that is
compatible with the current installed Python version), you can provide the
version as follows:

```terminal
$ conda install <package>=<version>
```

## Creating an environment

Conda lets you create environments, such that you can test things in a different
environment or group your software. Another common use case is to have different
environments for the different Python versions. Since conda is Python-based,
conflicting packages will mostly struggle with the Python version.

By default, conda will install packages into its root environment. Please note
that software that does not depend on Python and is installed in the root
environment, is is available in all other environments.

To create a Python 2.7 environment and activate it, issue the following commands:

```terminal
$ conda create -n py27 python=2.7
$ source activate py27
(py27) $
```

From now on, conda will install packages into the `py27` environment when you issue
the `install` command. To switch back to the root environment, simply deactivate the
`py27` environment:

```terminal
(py27) $ source deactivate py27
$
```

## Recommended packages

| Program   | Version  | package           |
|---        |---       |---                |
| Python    | 3.5      | `python=3.5`      |
| snakemake | latest   | `snakemake drmaa` |
| samtools  | latest   | `samtools`        |

To get you going, issue the following command. This will install Python 3.5
into your root environment alongside with Snakemake and samtools. By default,
conda starts with Python 3.6, but most packages are not adapted to that Python
version yet.

```terminal
$ conda install python=3.5 snakemake drmaa samtools
```

Please also read [this document](../slurm/snakemake.md#snakemake-and-slurm) on how to use Snakemake with DRMAA.

