# SSH Basics
## What is SSH?
SSH stands for **S** ecure **Sh** ell. It is a software that allows to establish a user-connection to a remote UNIX/Linux machine over the network and remote-control it from your local work-station.

Let's say you have an HPC cluster with hundreds of machines somewhere in a remote data-center and you want to connect to those machines to issue commands and run jobs. Then you would use SSH.

## Getting Started
### Installation
Simply install your distributions `openssh-client` package. You should be able to find plenty of good tutorials online.
On Windows you can consider using [MobaXterm (recommended)](../connecting/connecting-windows.md#install-ssh-client-for-windows) or [Putty](https://www.putty.org/).

### Connecting
Let's call your local machine the client and the remote machine you want to connect to the server.

You will usually have some kind of connection information, like a hostname, IP address and perhaps a port number. Additionally, you should also have received your user-account information stating your user-name, your password, etc.

Follow the instructions below to establish a remote terminal-session.

----

**If your are on Linux**

Open a terminal and issue the following command while replacing all the `<...>` fields with the actual data:

```
# default port
ssh <username>@<hostname-or-ip-address>

# non-default port
ssh <username>@<hostname-or-ip-address> -p <port-number>
```

**If you are on windows**

Start `putty.exe`, go into the `Session` category and fill out the form, then click the `Connect` button.
Putty also allows to save the connection information in different profiles so you don't have to memorize and retype all fields every time you want to connect.

----

## SSH-Keys
When you connect to a remote machine via SSH, you will be prompted for your password.
This will happen every single time you connect and can feel a bit repetitive at times, especially if you feel that your password is hard to memorize.
For those who don't want to type in their password every single time they connect, SSH keys are an alternative way of authentication.

Instead if being prompted for a password, SSH will simply use the key to authenticate.
As this key file should be device specific, this also increases security of the login process.

You can generate a new key by issuing:

```bash
client:~$ ssh-keygen -t ed25519

# 1. Choose file in which to save the key *(leave blank for default)*
# 2. Choose a passphrase of at least five characters
```

### How do SSH-Keys work?
An SSH key consists of two files, one private and one public key.
The public key is installed on remote machines and can only be validated with the matching private key, which is stored on client computers.
During the login process this is achieved via [public-key cryptography](https://en.wikipedia.org/wiki/Public-key_cryptography).

Traditionally the algorithm used for this was [RSA](https://en.wikipedia.org/wiki/RSA_(cryptosystem)).
Recently [elliptic curve cryptography](https://en.wikipedia.org/wiki/EdDSA) has been developed as a more secure and more performant alternative.
We recommend the `ed25519` type of SSH key.

### Passphrase
The security problem with SSH keys is that anyone with access to the private key has full access to all machines that have the public key installed.
Loosing the key or getting it compromised in another way imposes a serious security threat.
Therefore, it is best to secure the private key with a passphrase.
This passphrase is needed to unlock and use the private key.

Once you have your key-pair generated, you can easily change the passphrase of that key by issuing:

```bash
client:~$ ssh-keygen -p
```

### SSH-Agent
In order to avoid having to type the passphrase of the key every time we want to use it, the key can be loaded into an SSH-Agent.

For instance, if you have connected to a login-node via Putty and want to unlock your private key in order to be able to access cluster nodes, you cant configure the SSH-Agent.

```bash
client:~$ source <(ssh-agent)
```

*(The above command will load the required environment variables of the SSH-Agent into your shell environment, effectively making the agent available for your consumption.)*

Next, you can load your private key:

```bash
client:~$ ssh-add
```

*(You will be prompted for the passphrase of the key)*

You can verify that the agent is running and your key is loaded by issuing:

```bash
client:~$ ssh-add -l
# 'l' as in list-all-loaded-keys
```

*(The command should print at least one key, showing the key-size, the hash of the key-fingerprint and the location of the file in the file-system.)*

Since all home-directories are shared across the entire cluster and you created your key-pair inside your home-directory, you public-key (which is also in your home-directory) is automatically installed on all other cluster nodes, immediately.
Try connecting to any cluster node.
It should not prompt your for a password.

There is nothing you have to do to "unload" or "lock" the key-file.
Simply disconnect.