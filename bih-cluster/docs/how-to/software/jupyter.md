# How-To: Run Jupyter

!!! todo "Not yet updated to Slurm"

    TODO: This still needs to be updated to Slurm.

## What is Jupyter

[Project Jupyter](http://jupyter.org/) is a networking protocol for interactive computing that allows the user to write and execute code for a high number of different programming languages. The most used client is Jupyter Notebook that can be encountered in various form all over the web. Its basic principle is a document consisting of different cells, each of which contains either code (executed in place) or documentation (written in markdown). This allows one to handily describe the processed workflow.

## Setup and running Jupyter on the cluster

Install Jupyter on the cluster (via conda, by creating a custom environment)

```bash
med0xxx:~$ conda create -n jupyter jupyter
med0xxx:~$ source activate jupyter
```

(Install Jupyter language kernel as needed, see the kernel list below)

Now you can start the Jupyter server session (you may want to do this in a ```screen``` & ```qrsh``` session)
```bash
med0xxx:~$ jupyter notebook --no-browser
```

By default, Jupyter will create an access token (a link) to protect your notebook against unauthorized access which you have to save and enter in the accessing browser. You can change this to password base authorization via `jupyter notebook password`. If you are running multiple server on one or more nodes, one can separate them by changing the port number by adding `--port=$PORT`

## Connecting to the Running Session

This is slightly trickier as we have to create a SSH connection/tunnel via multiple hubs in between. The easiest is probably to configure your `.ssh/config` to automatically route your connection via the login-node (and MDC-jail if connecting from outside of the MDC network).

You have to add (replace with your `$CLUSTER_USER` and `$MDC_USER`)

```
Host med0*
  user $CLUSTER_USER
  ProxyCommand ssh $CLUSTER_USER@med-login2.bihealth.org -W %h:%p
```

(and if you are outside of the MDC network):

```
Host med-login2
  ProxyCommand ssh $MDC_USER@ssh1.mdc-berlin.de -W %h:%p
```

to your ```~/.ssh/config```

(If you have a newer version (7.2+) of SSH installed, you can also use `ProxyJump $CLUSTER_USER@med-login2.bihealth.org` instead of `ProxyCommand ...`)

See whether this works via i.e. `ssh med0110`

Now you setup a tunnel

```bash
workstation:~$ ssh -N -f -L 127.0.0.1:8888:localhost:${PORT} med0${NODE}
```

with the port of your server (usually `8888`) and the cluster node `qrsh` has send you to.

You should now be able to connect to the Jupyter server via `localhost:8888` in your webbrowser (see the note about token and password above) and start using Jupyter.

## Ending a Session

There are two independent steps in ending a session:

**Canceling the SSH tunnel**

- Identify the running SSH process

```bash
med0xxx:~$ ps aux | grep "$PORT"
```

This will give you something like this:

```
user        54  0.0  0.0  43104   784 ?        Ss   15:06   0:00 ssh -N -f -L 127.0.0.1:8888:localhost:8888 med0213
user        58  0.0  0.0  41116  1024 tty1     S    15:42   0:00 grep --color=auto 8888
```

from which you need the process ID (here `54`)

 - Terminate it the process

```bash
med0213:~$ kill -9 $PID
```

**Shutdown the Jupyter server**

Open the Jupyter session, cancel the process with {Ctrl} + {C} and confirm {y}. Make sure you saved your notebooks beforehand (though auto-save catches most things).

## Advanced

- [List of available Jupyter Kernels](https://github.com/jupyter/jupyter/wiki/Jupyter-kernels) for different programming languages
- [Jupyterlab](https://github.com/jupyterlab/jupyterlab) is a further development in the Jupyter ecosystem that creates a display similar to RStudio with panels for the current file system and different notebooks in different tabs.
- One can install Jupyter kernels or python packages while running a server or notebook without restrictions

If anyone has figured out, the following might also be interesting (please add):

- create a Jupyter-Hub
- multi-user support
