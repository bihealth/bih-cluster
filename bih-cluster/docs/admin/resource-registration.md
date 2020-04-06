# Resource Registration

This page describes the necessary registration steps for accessing special resources on the cluster, such as:

- GPU nodes,
- high-memory nodes,
- Matlab licenses,
- access to the `critical` partition.

Overall, there are no technical restrictions to accessing these resources.
Rather, we expect users to informally register with hpc-gatekeeper@bihealth.de and only then use the resources.
Registrations have to be **refreshed every 6 months**.
We trust that our users are grown-ups and are able to arrange for fair usage without administration intervention.

!!! warning "House Rules"
    Using the GPU node without prior arrangement with hpc-gatekeeper, not being reachable via email or phone within one working day while blocking the node, or other uncooperative behaviour can lead to HPC administration killing your jobs.

## GPU Nodes

> **Make sure to read the FAQ entry "[I have problems connecting to the GPU node! What's wrong?](Manual-Useful-Tips-Frequently-Asked-Questions#i-have-problems-connecting-to-the-gpu-node-whats-wrong)".**

1. If you want to use the GPU node, please send an email to hpc-gatekeeper@bihealth.de with the following information:
    - For how long do you want to have access?
      Please limit yourself to up to 6 months per request, you can extend the request afterwards with the same process.
    - How frequently do you need the GPU node, for how long at a time?
    - Do you plan to do more interactive work (e.g., using `ipython`) or using batch jobs?
2. At the moment, all requests will be granted by hpc-gatekeeper.
3. Use the GPU node by following the instructions [How To: Connect to GPU Nodes](How-To-Connect-to-GPU-Nodes)
4. **Be nice and cooperative with other users, e.g., ask arrange sharing of the node via email and phone.**
   Type `getent passwd USER_NAME` on the cluster to see user's contact details.

## High-Memory Nodes

Note that the purpose of the high memory nodes is to run jobs that don't run on the remainder of the cluster.
As the normal cluster nodes have 126-189GB of RAM each, we expect many jobs to fit on the (plenty) cluster nodes and these don't have to run on the (few and sparse) high memory nodes.

1. If you want to use the high memory nodes, please send an email to hpc-gatekeeper@bihealth.de with the following information:
    - For how long do you want to have access?
      Please limit yourself to up to 6 months per request, you can extend the request afterwards with the same process.
    - How frequently do you need the nodes, for how long at a time?
    - Do you plan to do more interactive work (e.g., using `ipython`) or using batch jobs?
    - What kind of work do you plan to do (e.g., assembly, running R scripts that don't run on the usual nodes)
2. At the moment, all requests will be granted by hpc-gatekeeper.
3. Use the high-memory node by following the instructions [How-To: Connect to High-Memory Node](How-To-Connect-to-High-Memory-Nodes)
4. **Be nice and cooperative with other users, e.g., ask arrange sharing of the node via email and phone.**
   Type `getent passwd USER_NAME` on the cluster to see user's contact details.

## Matlab Licenses

!!! note "GNU Octave as Matlab alternative"
    Note that [GNU Octave](https://www.gnu.org/software/octave/) is an Open Source alternative to Matlab.
    While both packages are not 100% compatible, Octave is an alternative that does not require any license management.
    Further, you can [easily install it yourself using Conda](Manual-Software-Management-Software-Installation-with-Conda).


!!! question "Want to use the Matlab GUI?"
    Make sure you understand X forwarding as outline [in this FAQ entry](Manual-Useful-Tips-Frequently-Asked-Questions#how-can-i-access-graphical-user-interfaces-such-as-for-matlab-on-the-cluster).

1. If you want to use the Matlab nodes, please send an email to hpc-gatekeeper@bihealth.de with the following information:
    - For how long do you want to have access?
      Please limit yourself to up to 6 months per request, you can extend the request afterwards with the same process.
    - How frequently do you need the nodes, for how long at a time?
    - Do you plan to do more interactive work or using batch jobs?
2. At the moment, all requests will be granted by hpc-gatekeeper.
3. Use the high-memory node by following the instructions [How-To: Use Matlab](How-To-Use-Matlab)
4. **Be nice and cooperative with other users, e.g., ask arrange sharing of the node via email and phone.**
   Type `getent passwd USER_NAME` on the cluster to see user's contact details.

!!! info "Requesting Licenses"
    Before using matlab, you have to request a license by passing `-l license_matlab_r2016b=1` to your `qrsh`, `qlogin` or `qsub` job.
    Failure to do so can lead to other user's jobs crashing because license are not available.
    Violations of this rule can lead to HPC administration killing your jobs.

## Critical Partition

The cluster provides a `critical` partition for jobs with deadlines (e.g., on paper submission or for tasks supporting clinical applications).
Access to this partition has to be explicitely granted as with the other resources on this page.

1. If you want to use the `critical` partition, please send an email to hpc-gatekeeper@bihealth.de with the following information:
    - For how long do you want to have access?
      Please limit yourself to up to 6 months per request, you can extend the request afterwards with the same process.
    - How frequently do you need the nodes, for how long at a time?
    - Do you plan to do more interactive work or using batch jobs?
2. At the moment, all requests will be granted by hpc-gatekeeper.
3. Use the high-memory node by following the instructions [How-To: Use the Critical Partition](How-To-Use-Critical-Partition)
4. **Be nice and cooperative with other users, e.g., ask arrange sharing of the node via email and phone.**
   Type `getent passwd USER_NAME` on the cluster to see user's contact details.
