# Debugging Connection Problems

When you encounter problems with the login to the cluster although we indicated
that you should have access, depending on the issue, here is a list of how to
solve the problem:

### I'm getting a "connection refused"

The full error message looks as follows:

```
ssh: connect to host hpc-login-1.cubi.bihealth.org port 22: Connection refused
```

This means that your computer could not open a network connection to the server.

- HPC 4 Clinic is currently only available from the Charite (cabled) network.
- HPC 4 Research can be connected to from:
    - Charite (cabled) network
    - Charite VPN :point_right: **but only with [Zusatzantrag B](/connecting/from-external/#zusatzantrag-b-recommended)**. :point_left:
    - MDC (cabled) network
    - MDC VPN
    - BIH (cabled) network
- If you think that there is no problem with any of this then please include the output of the following command in your ticket (use the server that you want to read instead of `<DEST>`):
    - Linux/Mac
        ```
        ifconfig
        traceroute <DEST>
        ```
    - Windows
        ```
        ipconfig
        tracepath <DEST>
        ```

## I can connect, but it seems that my account has no access yet

```
You're logging into BIH HPC cluster! (login-1)

 ***Your account has not been granted cluster access yet.***

 If you think that you should have access, please contact
 hpc-helpdesk@bihealth.de for assistance.

 For applying for cluster access, contact hpc-gatekeeper@bihealth.de.

user@login-1's password:
```

!!! hint

    **This is the most common error**, and the main cause for this is a wrong username. Please take a couple of minutes to read the [article about how usernames are constructed](prerequisites.md#what-is-my-username)!

If you encounter this message **although we told you that you have access and you checked the username as mentioned above**,
please write to [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de),
always indicating the message you get and a detailed description of what you
did.

## I'm getting a passPHRASE prompt

```
You're logging into BIH HPC cluster! (login-1)

 *** It looks like your account has access. ***

 Login is based on **SSH keys only**, if you are getting a password prompt
 then please contact hpc-helpdesk@bihealth.de for assistance.

Enter passphrase for key '/home/USER/.ssh/id_rsa':
```

Here you have to enter the **passphrase that was used for encrypting your private key**.
Read [SSH Basics](/misc/ssh-basics/) for further information of what is going on here.

## I can connect, but I get a passWORD prompt

```
You're logging into BIH HPC cluster! (login-1)

 *** It looks like your account has access. ***

 Login is based on **SSH keys only**, if you are getting a password prompt
 then please contact hpc-helpdesk@bihealth.de for assistance.

user@login-1's password:
```

!!! important "This is diffeerent from passPHRASE prompt"

    Please see [I'm getting a passPHRASE prompt](#im-getting-a-passphrase-prompt) for more information.

When you encounter this message during a login attempt, there is an issue with
your SSH key. In this case, please connect with increased verbosity to the
cluster (`ssh -vvv ...`) and mail the output and a detailed description to
[hpc-helpdesk@bihealth.de](mailto:hpc-helpdesk@bihealth.de).

