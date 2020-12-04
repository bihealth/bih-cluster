This section offers help for common problems in connecting to the cluster.

### I'm getting a "connection refused"

The full error message looks as follows:

```
ssh: connect to host login-1.research.hpc.bihealth.org port 22: Connection refused
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
