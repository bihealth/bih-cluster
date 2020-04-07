1. You must have a valid Charite or MDC username.
2. You have applied for and been granted access to the cluster by the [gatekeeper](mailto:hpc-gatekeeper@bihealth.de).
3. You have an SSH key [generated](../generate-key/linux.md) and [submitted](../submit-key/charite.md).

## What is my username?

The username for accessing the cluster is composed of your username at your primary organization (Charite/MDC) and a suffix:

- Charite user: `<Charite username>_c`, e.g. `doej_c`
- MDC user: `<MDC username>_m`, e.g. `jdoe_m`

!!! hint
    There are two login nodes, `med-login1` and `med-login2`. There are two for
    redundancy reasons. Please do not perform big file transfers or an `sshfs`
    mount via the login nodes. For this purpose, we have `med-transfer1` and
    `med-transfer2`.

!!! note
    Please Note that the cluster login is independent of access to the MDC jail node ssh1.mdc-berlin.de.

    - Access to the cluster is granted by BIH HPC IT through hpc-gatekeeper@bihealth.de
    - Access to the MDC jail node is managed by MDC IT