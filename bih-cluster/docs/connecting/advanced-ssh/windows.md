## Mounting the FS from within the MDC/Charite Network

!!! danger

    Mounting ssh on Windows is currently discouraged since relevant software is outdated
    (see also [hpc-talk](https://hpc-talk.cubi.bihealth.org/t/connecting-via-sshfs-win-manager/931)).
    Also, in most cases it is not really necessary to have a constant mount.
    For normal data transfer please use [WinSCP](../connecting-windows.md/#software-for-transfering-data-fromto-windows) instead.

Once WinSshFS is started, an icon will be added to your taskbar:

![](figures/winsshfs1.png)

Left-clicking that icon will bring up a window.
If not, right click the taskbar icon, select `Show Manager` and click `Add` in the menu.

![](figures/winsshfs2.png)

Fill out the marked fields:

- **Drive Name:**
Name that will show up in the windows explorer
- **Host:**
`hpc-transfer-1.cubi.bihealth.org`
- **Username:**
Your cluster username
- **Authentication method:**
`PrivateKey`. Select the `id_rsa` private key, not the `.ppk` format that is provided by PuTTY.
Enter the password that you used to secure your key with.
- **Directory:**
Cluster directory that will be mounted, you can choose any directory you have access to on the cluster.

Then click `Save` and then `Mount`.

Open the explorer. A new drive with the name you gave should show up:

![](figures/winsshfs3.png)

Finished!


## Connecting via MDC Jail Node

* This requires an active MDC account!

1. Additional to the steps above, click on the tab `Network settings`.
2. Check **Connect through SSH gateway (jump host)** and in the text field **Gateway SSH server** enter `ssh1.mdc-berlin.de` and in the field **User** your MDC username.
3. Check **Use private key** and select the SSH key that you uploaded to the MDC persdb (this might differ from your cluster key!).
4. Click **OK**

![](figures/mobaxterm_connect5.png)

## X11

!!! question "Do you really need to run a graphical application on the cluster?"

    Please note that running more complex Java applications, such as IGV may be not very efficient because of the connection speed. In most cases you can run them on your local workstation by [mounting them via SSHFS](#file-system-mount-via-sshfs).

Start MobaXterm, it should automatically fetch your saved Putty sessions as you can see on screen below:

![](figures/mobaxterm_main.png)

Connect to one of the login nodes, by double-click on saved profile, and then use `srun --pty --x11 bash` command to start X11 session to one of the nodes:

![](figures/mobaxterm_login.png)

Finally, start X11 application (below example of starting Visual Terminal):

![](figures/mobaxterm_xterm.png)
