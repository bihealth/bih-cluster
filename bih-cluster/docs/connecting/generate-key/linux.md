# Generating an SSH Key in Linux

- You might already have one, check whether a file `~/.ssh/id_xxx.pub` is present.
- Otherwise, create key using the following command (marking your key with your email address will make it easier to reidentify your key later on):
  ```shell
  $ ssh-keygen -t ed25519 -C "your_email@example.com"
  ```
- Use the default location for your key
- Enter a passphrase twice to encrypt your key

!!! important "What is a key passphrase?"

    You should set a passphrase when generating your key pair.
    It is used for encrypting your private key in case it is stolen or lost.
    When using the key for login, you will have to enter the passphrase.
    Many desktop environments offer ways to automatically unlock your key on login.

    Read [SSH Basics](../ssh-basics.md) for more information.

The whole session should look something like this:

```shell
host:~$ ssh-keygen -t ed25519 -C "your_email@example.com"
Generating public/private ed25519 key pair.
Enter file in which to save the key (/home/USER/.ssh/id_ed25519): 
Created directory '/home/USER/.ssh'.
Enter passphrase (empty for no passphrase):
Enter same passphrase again: 
Your identification has been saved in /home/USER/.ssh/id_ed25519.
Your public key has been saved in /home/USER/.ssh/id_ed25519.pub.
The key fingerprint is:
SHA256:Z6InW1OYt3loU7z14Kmgy87iIuYNr1gJAN1tG71D7Jc your_email@example.com
The key's randomart image is:
+--[ED25519 256]--+
|.. . . o         |
|. . . + +        |
|.    . = . .     |
|.     . +oE.     |
|.       So= o o  |
| . .   . * = + + |
|  +   o + B o o .|
| oo+. .B + + .   |
|.ooooooo*.  .    |
+----[SHA256]-----+
```

The file content of `~/.ssh/id_ed25519.pub` should look something like this):

```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIFzuiaSVD2j5y6RlFxOfREB/Vbd+47ABlxF7du5160ZH your_email@example.com
```

## Submit Your Key

As a next step you need to submit the SSH key use these links as:


- [:hospital: Charite user](../submit-key/charite.md)
- [:microscope: MDC user](../submit-key/mdc.md)
