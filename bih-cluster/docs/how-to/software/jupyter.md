# How-To: Run Jupyter

!!! info "SSH Tunnels Considered Harmful"

    **Please use our [Open OnDemand Portal](../../ondemand/overview.md) for running Jupyter notebooks!**

    The information below is still accurate.
    However, many users find it tricky to get SSH tunnels working correctly.
    A considerable number of parts is involved and you have to get each step 100% correct.
    Helpdesk cannot support you in problems with SSH tunnels that are caused by incorrect usage.

## What is Jupyter

[Project Jupyter](http://jupyter.org/) is a networking protocol for interactive computing that allows the user to write and execute code for a high number of different programming languages. The most used client is Jupyter Notebook that can be encountered in various form all over the web. Its basic principle is a document consisting of different cells, each of which contains either code (executed in place) or documentation (written in markdown). This allows one to handily describe the processed workflow.

## Setup and running Jupyter on the cluster

Install Jupyter on the cluster (via conda, by creating a custom environment)

```bash
med0xxx:~$ conda create -n jupyter jupyter
med0xxx:~$ conda activate jupyter
```

(If you want to work in a language other than python, you can install more Jupyter language kernel, see the [kernel list](https://github.com/jupyter/jupyter/wiki/Jupyter-kernels))

Now you can start the Jupyter server session (you may want to do this in a ```screen``` & ```srun --pty bash -i``` session as jupyter keeps running while you are doing computations)
```bash
med0xxx:~$ jupyter notebook --no-browser
```

Check the port number (usually `8888`) in the on output and remember it for later:
```bash
[I 23:39:40.860 NotebookApp] The Jupyter Notebook is running at:
[I 23:39:40.860 NotebookApp] http://localhost:8888/
```

By default, Jupyter will create an access token (a link stated in the output) to protect your notebook against unauthorized access which you have to save and enter in the accessing browser. You can change this to password base authorization via `jupyter notebook password`.
If you are running multiple server on one or more nodes, one can separate them by changing the port number by adding `--port=$PORT`.

## Connecting to the Running Session

This is slightly trickier as we have to create a SSH connection/tunnel via multiple hubs in between. The easiest is probably to configure your `.ssh/config` to automatically route your connection via the login-node (and MDC-jail if connecting from outside of the MDC network).

You have to add (replace with your `$CLUSTER_USER` and `$MDC_USER`)

```
Host med0*
  user $CLUSTER_USER
  ProxyCommand ssh $CLUSTER_USER@login-2.research.hpc.bihealth.org -W %h:%p
```

(and if you are outside of the MDC network):

```
Host login-2.research.hpc.bihealth.org
  ProxyCommand ssh $MDC_USER@ssh1.mdc-berlin.de -W %h:%p
```

to your ```~/.ssh/config```

(If you have a newer version (7.2+) of SSH installed, you can also use `ProxyJump $CLUSTER_USER@login-2.research.hpc.bihealth.org` instead of `ProxyCommand ...`)

See whether this works via i.e. `ssh med0110`

!!! info "SSH Key Forwarding"

    Please note that only the login nodes query the central user databases ("Active Directory Servers") for keys.
    The compute nodes themselves are not reachable directly from the outside and use the usual `~/.ssh/authorized_keys` file.
    This means that you will have to append the content of ~/.ssh/id_rsa.pub` from the source host (your workstation/laptop) to `~/.ssh/authorized_keys` on the cluster.
    
    If you fail to do this then you will be able to login to login-1/login-2 but you will get a password prompt for the compute nodes (i.e., med0XXX).

Now you setup a tunnel

```bash
workstation:~$ ssh -N -f -L 127.0.0.1:8888:localhost:${PORT} med0${NODE}
```

with the port of your server (usually `8888`) and the cluster node `srun` has send you to.


You should now be able to connect to the Jupyter server via `localhost:8888` in your webbrowser (see the note about token and password above) and start using Jupyter.

## Losing connection

It can and will happen that will lose connection, either due to network problems or due to shut-down of your computer.
This is not a problem at all and you will not lose data, just reconnect to your session.
If your notebooks are also losing connection (you will see a colorful remark in the top right corner), reconnect and click the colorful button.
If this does not work, your work is still not lost as all cells that have been executed are automatically saved anyways.
Copy all unexecuted cells (those are only saved periodically) and reload the browser page (after reconnecting) with `F5`.
(you can also open a copy of the notebook in another tab, just be aware that there may be synchronisation problems)

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
