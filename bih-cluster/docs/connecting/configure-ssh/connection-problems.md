# Debugging Connection Problems

When you encounter problems with the login to the cluster although we indicated
that you should have access, depending on the issue, here is a list of how to
solve the problem:

## I'm unable to connect to the cluster

When you can't reach the cluster, please make sure that you are in the right
network. Connections are only possible from within **Charité or MDC network**.
Connections from Eduroam or outside do not work unless you access via the MDC
jailnode (only for users with active MDC account).

**Solution**: Connect from within Charite or MDC.

## I can connect, but it seems that my account has no access yet

```
You're logging into BIH HPC cluster! (med-login1)

 ***Your account has not been granted cluster access yet.***

 If you think that you should have access, please contact
 hpc-helpdesk@bihealth.de for assistance.

 For applying for cluster access, contact hpc-gatekeeper@bihealth.de.

user@med-login1's password:
```

!!! hint

    **This is the most common error**, and the main cause for this is a wrong username. Please take a couple of minutes to read the [article about how usernames are constructed](prerequisites.md#what-is-my-username)!

If you encounter this message **although we told you that you have access and you checked the username as mentioned above**,
please write to [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de),
always indicating the message you get and a detailed description of what you
did.

## I can connect, but I get a password prompt

```
You're logging into BIH HPC cluster! (med-login1)

 *** It looks like your account has access. ***

 Login is based on **SSH keys only**, if you are getting a password prompt
 then please contact hpc-helpdesk@bihealth.de for assistance.

user@med-login1's password:
```

When you encounter this message during a login attempt, there is an issue with
your SSH key. In this case, please connect with increased verbosity to the
cluster (`ssh -vvv ...`) and mail the output and a detailed description to
[hpc-helpdesk@bihealth.de](mailto:hpc-helpdesk@bihealth.de).

!!! note "If you are asked to write to helpdesk"

    If you have an existent conversation with hpc-gatekeeper@bihealth.de and you are asked to open a ticket with hpc-helpdesk@bihealth.de: Please don't continue your conversation with hpc-helpdesk@bihealth.de, using the same ticket (i.e. same subject). The ticket system will automatically assign your reply to the existing ticket in hpc-gatekeeper. To avoid this, simply start with a fresh email.
