# Connecting to HPC 4 Research
HPC 4 Research is only available via the Charité, MDC, and BIH internal networks.
VPN access requires additional measures which are described in [Connecting from External Networks](./from-external.md).

There are two primary methods for interacting with BIH HPC:

1. Through the “Ondemand” web portal.
2. Via SSH and Slurm.

This part of the documentation only described direct console access via SSH.
For information regarding the web portal, please read [OnDemand Portal](../ondemand/overview.md).
In case you're not familiar with SSH, you should probably start via the web portal or (if you are determined to learn) read through our [SSH basics](ssh-basics.md) page.

## In brief
Follow these steps to connect to BIH HPC via the command line:

1. [Register an account](../admin/getting-access.md) via your PI.
2. Check your local PC for a SSH key pair or generate a new one. We have guides for [Linux](./generate-key/linux.md) and [Windows](./generate-key/windows.md) to help you out.
3. Submit your public key via [Charité](./submit-key/charite.md) or [MDC](./submit-key/mdc.md) web portal.
4. Connect to one of the two login nodes.
```bash
# Charite Users
$ ssh user_c@hpc-login-1.cubi.bihealth.org
$ ssh user_c@hpc-login-2.cubi.bihealth.org
# MDC Users
$ ssh user_m@hpc-login-1.cubi.bihealth.org
$ ssh user_m@hpc-login-2.cubi.bihealth.org
```

!!! hint
    There are two login nodes, `hpc-login-1` and `hpc-login-2`. There are two for
    redundancy reasons. Please do not perform big file transfers or an `sshfs`
    mount via the login nodes. For this purpose, we have `hpc-transfer-1` and
    `hpc-transfer-2`.

Please also read [Advanced SSH](./advanced-ssh/overview.md) for more custom scenarios how to connect to BIH HPC.
If you are using a Windows PC to access BIH HPC, please read [Connecting via SSH on Windows](./connecting-windows.md)

## What is my username?
Your username for accessing the cluster are composed of your username at your primary organization (Charité/MDC) and a suffix:

- Charite user: `<Charite username>_c -> doej_c`
- MDC user: `<MDC username>_m -> jdoe_m`

## How can I connect from the outside?
Please read [Connecting from External Networks](./from-external.md)

## I have problems connecting
Please read [Debugging Connection Problems](./connection-problems.md)