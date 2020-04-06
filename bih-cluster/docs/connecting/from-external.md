# Connecting from External Networks

This page describes how to connect to the BIH HPC from external networks (e.g., another university or from your home).
The options differ depending on your home organization and are described in detail below.

- MDC :microscope: users can use
    - the MDC SSH gateway/hop node, or
    - MDC VPN.
- Charite :hospital: users can use
    - the Charite VPN with "VPN Zusatzantrag B".

!!! faq "Getting Help with VPN and Gateway Nodes"

    Please note that the VPNs and gateway nodes are maintained by the central IT departments of Charite/MDC.
    BIH HPC IT cannot assist you in problems with these serves.
    Authorative information and documentation is provided by the central IT departments as well.

!!! tip "SSH Key Gotchas"

    You should use separate SSH key pairs for your workstation, laptop, home computer etc.
    As a reminder, you will have to [register the SSH keys](/connect/submit-key) with your home IT organization.
    When using gateway nodes, please make sure to use SSH key agents and agent forwarding (`ssh` flag "`-A`").

## For MDC Users

The general prerequisite is to [register the SSH keys](/connect/submit-key) with MDC IT via "persdb".

### MDC Gateway Node

The host name of the MDC gateway node is `ssh1.mdc-berlin.de`.
You will connect to this node with your **plain MDC user name**, e.g., `doej` **without** the "`_m`" suffix.
Do not forget the "`-A`" flag for SSH agent key forwarding.
Once you are on the gateway node, connect as if you were on your workstation:

```bash
# SSH key agent must be active at this point!
host:~$ ssh -A -l user ssh1.mdc-berlin.de
# If the SSH key agent does not run on your client host then the following
# will not work and the SSH key will not be available!
user@sl-it-p-ssh1:~$ ssh -l user_m med-login1.bihealth.org
```

### MDC VPN

You can find the instructions for getting MDC VPN access [here in the MDC intranet](https://www.mdc-berlin.info/anleitungen) below the "VPN" heading.
Please contact helpdesk@mdc-berlin.de for getting VPN access.

Install the VPN client and then start it.
Once VPN has been activated you can SSH to the HPC just as from your workstation.

```bash
host:~$ ssh -l user_m med-login1.bihealth.org
```

## For Charite Users

The general prerequisite is to [register the SSH keys](/connect/submit-key) with Charite IT via [zugang.charite.de](https://zugang.charite.de).

You will then have to apply for (1) general VPN access and (2) extended VPN access to BIH HPC.
Finally, you will be able to connect to BIH HPC from VPN.

### General Charite VPN Access

You need to apply for Charite VPN access, if you haven't done so already.
The form can be found in the [Charite Intranet](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Antrag_Mitarb_Stud.pdf) and contains further instructions.

### Zusatzantrag B (Recommended)

You can find [Zusatzantrag B](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Zusatzantrag_B.pdf) in the Charite intranet.
Fill it out and ship it in addition to the general VPN access form from above.
[Charite Helpdesk](helpdesk@charite.de) can help you with any questions.

Once you have been granted VPN access, start the client and connect to VPN.
You will then be able to connect from your client in the VPN just as you do from your workstation.

```
host:~$ ssh -l jdoe_c med-login1.bihealth.org
```

### Charite VDI

Alternative to using Zusatzantrag B, you can also get access to the Charite VDI (Virtual Desktop Infrastructure).
Here, you connect to a virtual desktop computer which is in the Charite network.
From there, you can connect to the BIH HPC system.

You need to apply for extended VPN access to be able to access the BIH VDI.
The form, partly pre-filled, can be found [here](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Zusatzantrag_O.pdf).
It is important to tick **Dienst(e)**, enter **HTTPS** and as target `view.bihealth.org`.
Please write to [bih-it-team@bihealth.de](bih-it-team@bihealth.de) with the request to access the [BIH VDI](https://view.bihealth.org).

When the access has been set up, follow the instructions on [client configuration](/connect/configure-ssh) for Windows, after logging in to the [BIH VDI](https://view.bihealth.org).
