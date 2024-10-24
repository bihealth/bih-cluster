# Connecting via SSH with Unix-like Operating Systems
## Activating your Key in the SSH Key Agent
!!! note
    The big Linux distributions automatically manage ssh-agent for you and unlock your keys at login time.
    If this doesn't work for you, read on.

`ssh-agent` caches your SSH keys so that you do not need to type your passphrase every time it is used.
Activate it by making sure `ssh-agent` runs in the background and add your key:

```bash
$ eval "$(ssh-agent -s)"
$ ssh-add
```

or if you chose a custom key name, specify the file like so:

```bash
$ ssh-add ~/.ssh/mdc_id_rsa
```

### MacOS 
If you run into problems that your key is not accepted when connecting from MacOS,
please use:

```bash
$ ssh-add --apple-use-keychain
```
## Configure SSH Client
You can define a personal SSH configuration file to make connecting to the cluster more comfortable by reducing the typing necessary by a lot.
Add the following lines to the file `~/.ssh/config` file.
Replace `USER_NAME` with your cluster user name.
You can also adapt the Host naming as you like.

```
Host bihcluster
    HostName hpc-login-1.cubi.bihealth.org
    User USER_NAME

Host bihcluster2
    HostName hpc-login-1.cubi.bihealth.org
    User USER_NAME
```

Now, you can do type the following (and you don't have to remember the host name of the login node any more).

```bash
$ ssh bihcluster
```

This configuration works if you are inside Charité, the Charité VPN, or MDC.

### MDC users: Jail node
If you have an MDC user account and want to connect from the outside, you can use the following `~/.ssh/config` lines to set up a ProxyJump via the MDC SSH jail.

```
Host mdcjail
    HostName ssh1.mdc-berlin.de
    User MDC_USER_NAME
```

Now you can run

```bash
$ ssh -J mdcjail bihcluster1
```

If you are always connecting from outside the internal network, you can also add a permanent ProxyJump to the SSH configuration like so:

```
Host bihcluster
    HostName hpc-login-1.cubi.bihealth.org
    User USER_NAME
    ProxyJump mdcjail
```

## Connecting with another computer/laptop
If you need to connect to the cluster from another computer than the one
that contains the SSH keys that you submitted for the cluster login, you
have two possibilities.

1. Generate another SSH key pair and submit the public part as described
   beforehand.
2. Copy your _private_ part of the SSH key (`~/.ssh/id_rsa`) to the second
   computer into the same location.

!!! danger

    Do not leave the key on any USB stick. Delete it after file transfer.
    This is a sensible part of data. Make sure that the files are only readable for you.

```
$ cd ~/.ssh
$ chmod g-rwx id_rsa*
$ ssh-add  id_rsa
```

## File System mount via sshfs
```
$ sshfs <USERNAME>@hpc-transfer-1.cubi.bihealth.org:/ <MOUNTPOINT>
```

* `hpc-transfer-1:` follows the structure `<host>:<directory>` starting in the user home.
* `<MOUNTPOINT>` must be an empty but existing and readable directory on your local computer

### MacOS
Make sure you have both OSXFUSE and SSHFS installed. You can get both from here: https://osxfuse.github.io/ or the most recent version via Homebrew:
```
$ brew cask install osxfuse; brew install sshfs; brew link --overwrite sshfs
```
The last command is optional and unlinks any pre-existing links to older versions of sshfs.
Now you can run
```
$ sshfs -o follow_symlinks <USERNAME>@hpc-transfer-1<X>.cubi.bihealth.org:<directory_relative_to_Cluster_root> <MOUNTPOINT> -o volname=<BIH-FOLDER> -o allow_other,noapplexattr,noappledouble
```

## X11

!!! question "Do you really need to run a graphical application on the cluster?"

    Please note that running more complex Java applications, such as IGV may be not very efficient because of the connection speed. In most cases you can run them on your local workstation by [mounting them via SSHFS](#file-system-mount-via-sshfs).

Connect to one of the login nodes using X11 forwarding:

```bash
$ ssh -X -C -t <USERNAME>@hpc-login-1.bihealth.org
```

Once you get a login prompt, you can use the `srun` command with the `--x11` parameter to open a X11 session to a cluster node:

```bash
$ srun --pty --x11 bash
```

And finally you can start your X11 application, e.g.:

```bash
$ xterm
```

After a while Visual Terminal should start:

![](figures/xterm_linux.png)
