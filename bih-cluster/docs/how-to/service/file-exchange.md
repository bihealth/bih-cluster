# How-To: Use File Exchange

BIH HPC IT provides a file exchange server to be used by the BIH Omics facilities and their users.
The server is located in the BIH DMZ in Buch.
Clients from BIH and MDC network have direct access and clients from the outside can be granted temporary access.

## System Access

:exclamation:
There are two caveats for accessing the server.
:exclamation:

1. Currently, **access from the Charite network is not possible**.
   For transfering data to the Charite network, first copy it to the cluster and then copy it to your workstation via `scp` or `rsync` via `med-transfer1` or `med-transfer2`.
2. The BIH cluster is behind a firewall.
   You will only be able to open connections via FTPS to the machine from `med-transfer1.bihealth.org` and `med-transfer2`. **Access from the login and compute nodes will not work**.

## Account Creation

Accounts are created for each file exchange and are temporary, active for two weeks by default.
Users from the omics core units can request accounts by sending an email to hpc-admin@bihealth.de.
The process is as follows:

1. Send an email to hpc-admin@bihealth.de.
2. The request should have the following format:
    - Subject: File Exchange Account Creation
    - Body: Answer the following questions:
        - What is the (1) name, (2) email address, and (3) organization of the recipient.
        - Is the destination on the internal (BIH/Charite/MDC) networks or on the internet?
        - How much data is roughly to be transferred?
3. BIH HPC IT will create a temporary account with random password and send it to the requestor.
4. You upload the files (e.g., via `lftp`, see below).
5. For internal users, you send the user name and password to the data recipient who downloads the data (e.g., via `lftp`, see below).
6. For external users, also IP address of the user should be provided, so it can be allowed on internal firewall.
## Data Upload

The server can be uploaded to via FTPS.
E.g., using `lftp` on the `med-transfer1` node of the BIH cluster:

```terminal
~ # cd data-dir
data-dir # lftp
lftp: ~> open -u USER,PASSWORD file-exchange.bihealth.org
lftp USER@file-exchange.bihealth.org:~> ls
[connecting]
lftp USER@file-exchange.bihealth.org:/> mirror -R $"YOURFLOWCELL"
`FILE NAME' at 3676560 (0%) 751.8K/s eta:17h [Sending data/TLS]
[...]
lftp USER@file-exchange.bihealth.org:/> quit
data-dir #
```

## Data Download

E.g., you can download the data using `lftp` to the `med-transfer1` node of the BIH cluster via FTPS.

```
~ # data data-dir
data-dir # lftp
lftp: ~> open -u USER,PASSWORD file-exchange.bihealth.org
lftp USER@file-exchange.bihealth.org:~> mget *
[...]
lftp USER@file-exchange.bihealth.org:/> quit
data-dir #
```
