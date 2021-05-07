# Generating an SSH Key in Linux

- You might already have one, check whether the file `~/.ssh/id_rsa.pub` is present.
- Otherwise, create key using the following command (marking your key with your email address will make it easier to reidentify your key later on):
  ```shell
  host:~$ ssh-keygen -t rsa -C "your_email@example.com"
  ```
- Use the default location for your key
- Enter the passphrase twice to encrypt your key

!!! important "What is your key's passphrase?"

    You should set a passphrase when generating your private key.
    This passphrase is used for encrypting you private key to protect it against the private key file theft/being lost.
    When using the key for login, you will have to enter it (or the first time you load it into the SSH key agent).
    Note that when being asked for the **passphrase** this does not occur on the cluster (and is thus unrelated to it) but on your local computer.

    Also see [SSH Basics](/misc/ssh-basics/) for more information.

The whole session should look something like this:

```shell
host:~$ ssh-keygen -t rsa -C "your_email@example.com"
Generating public/private rsa key pair.
Enter file in which to save the key (/home/USER/.ssh/id_rsa): 
Created directory '/home/USER/.ssh'.
Enter passphrase (empty for no passphrase):
Enter same passphrase again: 
Your identification has been saved in /home/USER/.ssh/id_rsa.
Your public key has been saved in /home/USER/.ssh/id_rsa.pub.
The key fingerprint is:
55:dd:8f:88:84:1b:b6:f0:19:d3:fb:19:8e:7a:9e:7d your_email@example.com
The key's randomart image is:
+--[ RSA 2048]----+
|         o  .. . |
|      . * o.  . .|
|       + O.o . ..|
|        =.o o . .|
|        S  + o   |
|          . +    |
|         .       |
|        . .o  E  |
|         oo ..   |
```

The file content of `~/.ssh/id_rsa.pub` should look something like this):

```
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC/Rdd5rvf4BT38jsBlRrXpd1KDvjE1iZZlEmkB6809QK7hV6RCG13VcyPTIHSQePycfcUv5q1Jdy28MpacL/nv1UR/o35xPBn2HkgB4OqnKtt86soCGMd9/YzQP5lY7V60kPBJbrXDApeqf+H1GALsFNQM6MCwicdE6zTqE1mzWVdhGymZR28hGJbV9H4snMDDc0tW4i3FHGrDdmb7wHM9THMx6OcCrnNyA9Sh2OyBH4MwItKfuqEg2rc56D7WAQ2JcmPQZTlBAYeFL/dYYKcXmbffEpXTbYh+7O0o9RAJ7T3uOUj/2IbSnsgg6fyw0Kotcg8iHAPvb61bZGPOEWZb your_email@example.com
```

## Submit Your Key

As a next step you need to submit the SSH key use these links as:

- [:hospital: Charite user](../../submit-key/charite)
- [:microscope: MDC user](../submit-key/mdc.md)
