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

1. [Register an account](../admin/getting-access.md) via your PI. :memo: 
2. [Generate a SSH key pair :key: in Linux](generate-key/unix.md) or [Windows](generate-key/windows.md)
3. [Submit your public key :arrow_up: to Charite](submit-key/charite.md) or [to MDC](submit-key/mdc.md).
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

5. Allocate resources on a computation node using [Slurm](../slurm/overview.md). Do not compute on the login node!

    ```bash
    # Start interactive shell on computation node
    $ srun --pty bash -i
    ```

6. Bonus: [Configure your SSH client :wrench: on Linux and Mac](advanced-ssh/unix.md) or [Windows](advanced-ssh/windows.md).
7. Bonus: [Connect from external networks :flying_saucer:](./from-external.md).

!!! tip "tl;dr"

    - Web Access: https://hpc-portal.cubi.bihealth.org
    - SSH-Based Access:

        ```bash
        # Interactive login (choose one)
        ssh username@hpc-login-1.cubi.bihealth.org
        ssh username@hpc-login-2.cubi.bihealth.org
        srun --pty bash -i

        # File Transfer (choose one)
        sftp local/file username@hpc-transfer-1.cubi.bihealth.org:remote/file
        sftp username@hpc-transfer-2.cubi.bihealth.org:remote/file local/file

        # Interactive login into the transfer nodes (choose one)
        ssh username@hpc-transfer-1.cubi.bihealth.org
        ssh username@hpc-transfer-2.cubi.bihealth.org
        ```

## What is my username?
Your username for accessing the cluster are composed of your username at your primary organization (Charité/MDC) and a suffix:

- Charite user: `<Charite username>_c -> doej_c`
- MDC user: `<MDC username>_m -> jdoe_m`

## How can I connect from the outside?
Please read [Connecting from External Networks](./from-external.md)

## I have problems connecting
Please read [Debugging Connection Problems](./connection-problems.md)
