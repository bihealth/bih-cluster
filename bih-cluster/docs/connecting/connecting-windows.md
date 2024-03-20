# Connecting via SSH on Windows
## Install SSH Client for Windows
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

### Software for transfering data from/to Windows
For transfering data from/to Windows, we recommand using WinSCP.
Install the latest version from here: https://winscp.net/eng/download.php

On the `Login` screen of WinSCP create a new login by selecting `New Site`.

Fill in the following parameters:

* `File protocol`: `SFTP`
* `Host name`: `hpc-transfer-1.cubi.bihealth.org` or `hpc-transfer-2.cubi.bihealth.org`
* `User name`: your user name

Go to `Advanced` > `SSH` > `Authentication` > `Authentication parameters` > `Private key file`
and select your private ssh key file (in `.ppk` format).
 
Press `Ok` then `Save`.

Press `Login` to connect.
It will ask for your private key passphrase, if you set one up.

If you need to convert your private ssh key file the `.ppk` format,
on the WinSCP login screen go to `Tools` > `PuTTYgen` and follow the steps here:
https://docs.acquia.com/cloud-platform/manage/ssh/sftp-key/

## Connecting from within MDC/Charite Network

Click on `Session`.

![](figures/mobaxterm_connect1.png)

Click on `SSH`.

![](figures/mobaxterm_connect2.png)

In **Basic SSH settings**, enter a hostname (`hpc-login-X.cubi.bihealth.org`, where `X` is 1 or 2), check **Specify username** and enter your username in the textfield.
Select the tab **Advanced SSH settings**, check **Use private key** and select your private SSH key file (possible choices described with the next to figures).

![](figures/mobaxterm_connect3.png)

Select the `id_rsa` file generated in Linux OR

![](figures/mobaxterm_connect3a.png)

select the `id_rsa.ppk` file generated in Windows with MobaXterm.

![](figures/mobaxterm_connect3b.png)

Afterwards hit the **OK** button and MobaXterm will connect.

![](figures/mobaxterm_connect4.png)

The session will be stored automatically and you can establish new connections later on, or also multiple ones at the same time, if you like.