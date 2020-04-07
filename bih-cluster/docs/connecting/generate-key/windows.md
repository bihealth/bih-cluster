# Generating an SSH Key in Windows

!!! hint

    Please install an [SSH client for Windows](../ssh-client-windows.md).

## Generate the Key

Click on `Tools` and `MobaKeyGen (SSH key generator)`

![](figures/mobaxterm_keygen1.png)

In the section **Parameters** make sure to set the following properties:

* **Type of key to generate:** `RSA` (this is the `SSH-2` protocol)
* **Number of bits in a generated key:** `4096`

If all is set, hit the **Generate** button.

![](figures/mobaxterm_keygen2.png)

During generation, move the mouse cursor around in the blank area.

![](figures/mobaxterm_keygen3.png)

When finished, make sure to protect your generated key with a passphrase.
Save the private and public key. The default name under Linux for the public
key is `id_rsa.pub` and `id_rsa` for the private key, but you can name them
however you want (the `.pub` is NOT automatically added). Note that in the
whole cluster wiki we will use this file naming convention. Also note that the
private key will be stored in Putty format (`.ppk`, this extension is added
automatically).
**Important** The gibberish in the textbox is your public key in the
format how it has to be submitted to the MDC and Charite (links for this step
below). Thus, copy this text and paste it to the SSH-key-submission-web-service
of your institution.

![](figures/mobaxterm_keygen4.png)

**Important** Store the private key additionally in the OpenSSH format. To do
so, click `Conversions` and select `Export OpenSSH key`. To be consistent, give
the file the same name as your `.ppk` private key file above (just without the
`.ppk`).

![](figures/mobaxterm_keygen5.png)

## Summary

To summarize, you should end up with three files:

1. `id_rsa.pub`
The public key file, it is not required if you copy and submit the SSH public
key as described above and in the links below.
2. `id_rsa.ppk`
This file is only needed if you plan to use Putty.
3. `id_rsa`
This is your private key and the one and only most important file to access the
cluster. It will be added to the sessions in MobaXterm and WinSSHFS
(if required).

## Submit Your Key

As a next step you need to submit the SSH key use these links as:

- [:hospital: Charite user](../submit-key/charite)
- [:microscope: MDC user](../submit-key/mdc)
