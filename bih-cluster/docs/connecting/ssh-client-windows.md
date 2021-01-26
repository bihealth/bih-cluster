# Installing SSH Client for Windows

We recommend to use the program MobaXterm on Windows.
MobaXterm is a software that allows you to connect to an SSH server, much like PuTTy, but also maintains your SSH key.

!!! hint "Alternative SSH Clients for Windows"
    - Another popular option is [PuTTy](https://www.putty.org/) but many users have problems configuring it correctly with SSH keys.
    - On Windows 10, you can also install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10), e.g., together with [WSL Terminal](https://github.com/mskyaxl/wsl-terminal).
      This is not for the faint of heart (but great if you're a Unix head).

* Navigate to https://mobaxterm.mobatek.net/download-home-edition.html
* Download either the
    * **Portable edition** (blue button lefthand-side, **if you have no admin rights**, e.g. on a Charite or MDC workstation), or
    * **Installer edition** (green button righthand-side, **requires admin rights on your computer**).
* Install or unpack MobaXterm and start the software. As a Charite user, please cancel any firewall warnings that pop up.

## Software for Mounting the Cluster FS

In case you also want to mount the cluster file system in Windows, we recommend WinSshFS.

* Navigate to https://github.com/feo-cz/win-sshfs/releases/tag/1.5.12.8
* Download the **Release1.5.12.8.zip** file.
* Unpack and execute.
