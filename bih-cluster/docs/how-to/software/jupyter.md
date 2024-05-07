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

This is slightly trickier as we have to create a SSH connection/tunnel with potentially multiple hops in between. The easiest way is probably to configure your `.ssh/config` to automatically route your connection via the login node (and possibly MDC jail). This is described in our [Advanced SSH config documentation](/connecting/advanced-ssh/linux.md#configure-ssh-client).

In short,add these lines to `~/.ssh/config` (replace curly parts):

```
Host bihcluster
  user {USERNAME}
  HostName hpc-login-2.cubi.bihealth.org

Host hpc-cpu*
  user {USERNAME}
  ProxyJump bihcluster
```

For MDC users outside the MDC network:
```
Host mdcjail
    HostName ssh1.mdc-berlin.de
    User {MDC_USER_NAME}

Host bihcluster
  user {USERNAME}
  HostName hpc-login-2.cubi.bihealth.org

Host hpc-cpu*
  user {USERNAME}
  ProxyJump bihcluster
```

Check that this config is working by connecting like this: `ssh hpc-cpu-1`. Please note that you cannot use any resources on this node without a valid Slurm session.

Now you setup a tunnel for your running Jupyter session:

```bash
workstation:~$ ssh -N -f -L 127.0.0.1:8888:localhost:{PORT} hpc-cpu-{x}
```
The port of your Jupyter server is usually `8888`. The cluster node `srun` has sent you to determines the last argument.

You should now be able to connect to your Jupyter server by typing `localhost:8888` in your webbrowser (see the note about token and password above).

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
