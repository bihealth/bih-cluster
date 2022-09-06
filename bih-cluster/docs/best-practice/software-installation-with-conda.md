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

When you logged into the cluster, please make sure that you also executed `srun` to log into a computation node and perform the software installation there.

## Installing conda

```bash
res-login-1:~$ srun --mem=5G --pty bash -i
med0127:~$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
med0127:~$ bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/work/miniconda
```

This will install conda to `$HOME/work/miniconda`.
This path can be changed to your liking.
Please note that the `$HOME` folder has limited space (an exception is the subfolder `$HOME/work` which has no space limit).

NB: `$HOME/scratch` is not appropriate as files placed there will be removed automatically after 2 weeks.

To make it available upon login, extend and export the `$PATH` variable with the
installation path + `/bin` and add it to your `$HOME/.bashrc`:

```bash
case "${SLURMD_NODENAME-${HOSTNAME}}" in
    login-*)
        ;;
    *)
        export PATH=$HOME/work/miniconda/condabin:$PATH
        ;;
esac
```

The above code makes sure that you don't have conda available on the login nodes,
where you are not allowed to start any computations.

To make bioinformatics software available, we have to add the `bioconda` and
some other channels to the conda configuration:

```bash
med0127:~$ conda config --add channels bioconda
med0127:~$ conda config --add channels default
med0127:~$ conda config --add channels conda-forge
```

You can also add channels to your liking.

## Installing software with conda

Installing packages with conda is straight forward:

```bash
med0127:~$ conda install <package>
```

This will install a package into the conda root environment. We will explain
environments in detail in the next section.

To search for a package, e.g. to find the correct name in conda or if it exists
at all, issue the command:

```bash
med0127:~$ conda search <string>
```

To choose a specific version (conda will install the latest version that is
compatible with the current installed Python version), you can provide the
version as follows:

```bash
med0127:~$ conda install <package>=<version>
```

Please note that new conda installs may ship with a recently update Python version and not all packages might have been adapted.
E.g., if you find out that some packages don't work after starting out/upgrading to Python 3.8, simply try to downgrade Python to 3.7 with `conda install python=3.7`.

!!! hint
    As resolving the dependency tree of an installation candidate can take a lot of
    time in Conda, especially when you are installing software from an `environment.yaml`
    file, an alternative resolver has been presented that you can use to install
    software into your Conda environment. The time savings are immense and an
    installation that took more than an hour can be resolved in seconds.

    Simply run

    ```bash
    med0127:~$ conda install mamba
    ```

    With that, you can install software into your environment using the same syntax as for Conda:

    ```bash
    med0127:~$ mamba install <package>
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

```bash
med0127:~$ conda create -n py27 python=2.7
med0127:~$ source activate py27
(py27) med0127:~$
```

From now on, conda will install packages into the `py27` environment when you issue
the `install` command. To switch back to the root environment, simply deactivate the
`py27` environment:

```bash
(py27) med0127:~$ source deactivate py27
med0127:~$
```

But of course, as Python 2.7 is not supported any more by the Python Software Foundation, you should switch over to Python 3 already!
