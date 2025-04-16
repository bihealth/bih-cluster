# Software Installation with Conda
## Conda
Users do not have the rights to install system packages on the BIH HPC cluster.
For the management of bioinformatics software we therefore recommend using the conda package manager.
Conda provides software in different “channels” and one of those channels contains a huge selection of bioinformatics software (bioconda).
Generally packages are pre-compiled and conda just downloads the binaries from the conda servers.

You are in charge of managing your own software stack, but conda makes it easy
to do so. We will provide you with a description on how to install conda and how
to use it. Of course there are many online resources that you can also use.
Please find a list at the end of the document.

!!! Warning
    Following a change in their terms of service Anaconda Inc. has started to demand payment from research institutions
    for using both Anaconda, Miniconda, and the defaults channel. As a consequence, 
    usage of this software is prohibited and we're recommending the alternative free "miniforge" distribution instead.

## Premise
When you logged into the cluster, please make sure that you also executed `srun` to log into a computation node and perform the software installation there.

## Installing conda

```bash
hpc-login-1:~$ srun --mem=5G --pty bash -i
hpc-cpu-123:~$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
hpc-cpu-123:~$ bash Miniforge3-Linux-x86_64.sh -b -f -p $HOME/work/miniforge
hpc-cpu-123:~$ eval "$(/$HOME/work/miniforge/bin/conda shell.bash hook)"
hpc-cpu-123:~$ conda init
hpc-cpu-123:~$ conda config --set auto_activate_base false
```

This will install conda to `$HOME/work/miniforge`.
You can change the path to your liking, but please note that your `$HOME` folder has limited space.
The `work` subfolder however has a bigger quota. More about this [here](../storage/home-quota.md).

miniforge uses `conda-forge` as primary channel.

To make bioinformatics software available, we have to add the `bioconda` channel to the conda configuration.

Note, that adding a channel put's it on top of the channel list.
So we suggest to re-order/-prioritize channels by adding `conda-forge` again and set a strict channel priority
(as suggested by [bioconda](https://bioconda.github.io/#usage)).

```bash
hpc-cpu-123:~$ conda config --add channels bioconda
hpc-cpu-123:~$ conda config --add channels conda-forge
hpc-cpu-123:~$ conda config --set channel_priority strict
```

We could also edit `~/.condarc` manually.

Either way, the recommended initial `~/.condarc` should look like this (in that order!):

```
hpc-cpu-123:~$ cat ~/.condarc
channels:
  - conda-forge
  - bioconda
channel_priority: strict
```

!!! Warning
    If you are just looking to re-configure a previous miniconda installation,
    make sure to remove the Anaconda Inc. channels `defaults` and `r` from your `~/.condarc`
    to prevent SSL errors due the licencing issues mentiond above.

## Installing software with conda
Installing packages with conda is straight forward:

```bash
hpc-cpu-123:~$ conda install <package>
```

This will install a package into the conda base environment. 
We will explain environments in detail in the next section.
To search for a package, e.g. to find the correct name in conda or if it exists
at all, issue the command:

```bash
hpc-cpu-123:~$ conda search <string>
```

To choose a specific version (conda will install the latest version that is
compatible with the current installed Python version), you can provide the
version as follows:

```bash
hpc-cpu-123:~$ conda install <package>=<version>
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
    hpc-cpu-123:~$ conda install mamba
    ```

    With that, you can install software into your environment using the same syntax as for Conda:

    ```bash
    hpc-cpu-123:~$ mamba install <package>
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
hpc-cpu-123:~$ conda create -n py27 python=2.7
hpc-cpu-123:~$ source activate py27
(py27) hpc-cpu-123:~$
```

From now on, conda will install packages into the `py27` environment when you issue
the `install` command. To switch back to the root environment, simply deactivate the
`py27` environment:

```bash
(py27) hpc-cpu-123:~$ source deactivate py27
hpc-cpu-123:~$
```

But of course, as Python 2.7 is not supported any more by the Python Software Foundation, you should switch over to Python 3 already!
