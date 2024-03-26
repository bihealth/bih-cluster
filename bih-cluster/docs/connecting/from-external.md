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
    As a reminder, you will have to register the SSH keys with your home IT organization ([MDC](submit-key/mdc.md) or [Charite](submit-key/charite.md)).
    When using gateway nodes, please make sure to use SSH key agents and agent forwarding (`ssh` flag "`-A`").

## MDC Users
### Via Gateway Node
Use the following command to perform a proxy jump via the MDC SSH gateway (`ssh1` aka `jail1`) when connecting to a login node.
Note that for logging into the jail, the `<MDC_USER>` is required.

```bash
$ ssh -J <MDC_USER>@ssh1.mdc-berlin.de <HPC_USER>@hpc-login-1.cubi.bihealth.org
```

!!! note
    Please Note that the cluster login is independent of access to the MDC jail node ssh1.mdc-berlin.de.

    - Access to the cluster is granted by BIH HPC IT through hpc-helpdesk@bih-charite.de.
    - Access to the MDC jail node is managed by MDC IT.

### Via MDC VPN
You can find the instructions for getting MDC VPN access [here in the MDC intranet](https://www.mdc-berlin.info/anleitungen) below the "VPN" heading.
Please contact [helpdesk@mdc-berlin.de](mailto:helpdesk@mdc-berlin.de) for getting VPN access.

Install the VPN client and then start it.
Once VPN has been activated you can SSH to the HPC just as from your workstation.

```bash
$ ssh user_m@hpc-login-1.cubi.bihealth.org
```

## Charité Users
Access to BIH HPC from external networks (including Eduroam) requires a Charité VPN connection with special access permissions.

### General Charité VPN Access
You need to apply for general Charité VPN access if you haven't done so already.
The form can be found in the [Charite Intranet](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Antrag_Mitarb_Stud.pdf) and contains further instructions.
[Charité IT Helpdesk](mailto:helpdesk@charite.de) can help you with any questions.

### Zusatzantrag B
Special permissions form B is also required for HPC access.
You can find [Zusatzantrag B](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Zusatzantrag_B.pdf) in the Charité intranet.
Fill it out and send it to the same address as the general VPN access form above.

Once you have been granted VPN access, start the client and connect to VPN.
You will then be able to connect from your client in the VPN just as you do from your workstation.

```bash
$ ssh jdoe_c@hpc-login-1.cubi.bihealth.org
```

### Charité VDI (Not recommended)
Alternative to using Zusatzantrag B, you can also get access to the Charité VDI (Virtual Desktop Infrastructure).
Here, you connect to a virtual desktop computer which is in the Charité network.
From there, you can connect to the BIH HPC system.

You need to apply for extended VPN access to be able to access the BIH VDI.
The form can be found [here](https://intranet.charite.de/fileadmin/user_upload/portal/service/service_06_geschaeftsbereiche/service_06_14_it/VPN-Zusatzantrag_O.pdf).
It is important to tick **Dienst(e)**, enter **HTTPS** and as target `view.bihealth.org`.
Please write to [helpdesk@charite.de](mailto:helpdesk@charite.de) with the request to access the [BIH VDI](https://view.bihealth.org).

When the access has been set up, follow the instructions on [client configuration](advanced-ssh/windows.md) for Windows, after logging in to the [BIH VDI](https://view.bihealth.org).
