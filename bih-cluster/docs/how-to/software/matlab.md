# How-To: Use MATLAB

!!! todo "Matlab is not fully integrated yet"
    We also need to finalize installation and then document the compiler version

!!! note "GNU Octave as Matlab alternative"
    Note that [GNU Octave](https://www.gnu.org/software/octave/) is an Open Source alternative to Matlab.
    While both packages are not 100% compatible, Octave is an alternative that does not require any license management.
    Further, you can [easily install it yourself using Conda](../../best-practice/software-installation-with-conda).

!!! question "Want to use the Matlab GUI?"
    Make sure you understand X forwarding as outline [in this FAQ entry](../../help/faq.md#how-can-i-access-graphical-user-interfaces-such-as-for-matlab-on-the-cluster).

## Pre-requisites

You have to register with [hpc-gatekeeper@bihealth.de](mailto:hpc-gatekeeper@bihealth.de) for requesting access to the MATLAB licenses.
Register on [User Self Organisation MATLAB](../../admin/resource-registration.md#matlab-licenses) after registering with hpc-gatekeeper.

Afterwards, you can connect to the High-Memory using the `license_matlab_r2016b` resource (see below).

## How-To Use

BIH has a license of MATLAB R2016b for **16 seats** and various licensed packages (see below).
It is installed on all of the compute nodes:

```console
# The following is VITAL so the scheduler allocates a license to your session.
$ qrsh -l license_matlab_r2016b=1
$ module avail
----------------- /usr/share/Modules/modulefiles -----------------
dot         module-info null
module-git  modules     use.own

----------------------- /opt/local/modules -----------------------
cmake/3.11.0-0  llvm/6.0.0-0    openmpi/3.1.0-0
gcc/7.2.0-0     matlab/r2016b-0
$ module load matlab/r2016b-0 
Start matlab without GUI: matlab -nosplash -nodisplay -nojvm
    Start matlab with GUI (requires X forwarding (ssh -X)): matlab
$ matlab -nosplash -nodisplay -nojvm
                                               < M A T L A B (R) >
                                     Copyright 1984-2016 The MathWorks, Inc.
                                     R2016b (9.1.0.441655) 64-bit (glnxa64)
                                                September 7, 2016

 
For online documentation, see http://www.mathworks.com/support
For product information, visit www.mathworks.com.
 

	Non-Degree Granting Education License -- for use at non-degree granting, nonprofit,
	educational organizations only.  Not for government, commercial, or other organizational use.

>> ver
--------------------------------------------------------------------------------------------
MATLAB Version: 9.1.0.441655 (R2016b)
MATLAB License Number: 1108905
Operating System: Linux 3.10.0-862.3.2.el7.x86_64 #1 SMP Mon May 21 23:36:36 UTC 2018 x86_64
Java Version: Java is not enabled
--------------------------------------------------------------------------------------------
MATLAB                                                Version 9.1         (R2016b)
Bioinformatics Toolbox                                Version 4.7         (R2016b)
Global Optimization Toolbox                           Version 3.4.1       (R2016b)
Image Processing Toolbox                              Version 9.5         (R2016b)
Optimization Toolbox                                  Version 7.5         (R2016b)
Parallel Computing Toolbox                            Version 6.9         (R2016b)
Partial Differential Equation Toolbox                 Version 2.3         (R2016b)
Signal Processing Toolbox                             Version 7.3         (R2016b)
SimBiology                                            Version 5.5         (R2016b)
Statistics and Machine Learning Toolbox               Version 11.0        (R2016b)
Wavelet Toolbox                                       Version 4.17        (R2016b)
>> exit
```

## Running MATLAB UI

For starting the Matlab with GUI, make sure that your client is running a X11 server and you connect with X11 forwarding enabled (e.g., `ssh -X med-login1` from the Linux command line).
Then, make sure `qlogin` for connecting to a node with X11 forwarding enabled.

```bash
your-client $ ssh -X med-login1
[...]
med-login1 $ qlogin -l license_matlab_r2016b=1
[...]
med0203 $ module load matlab/r2016b-0
Start matlab without GUI: matlab -nosplash -nodisplay -nojvm
    Start matlab with GUI (requires X forwarding (ssh -X)): matlab
med0203 $ matlab
[UI will start]
```

For forcing starting in text mode can be done (as said after `module load`): `matlab -nosplash -nodisplay -nojvm`.

Also see [this FAQ entry](../../help/faq.md#how-can-i-access-graphical-user-interfaces-such-as-for-matlab-on-the-cluster).

## See Available Matlab Licenses

You can use `qstat -F license_matlab_r2016b` to see the currently available MATLAB license.
E.g., here I am running an interactive shell in which I have requested 1 of the 16 MATLAB licenses, so 15 more remain.

```
$ qstat -F license_matlab_r2016b  | head
queuename                      qtype resv/used/tot. load_avg arch          states
---------------------------------------------------------------------------------
all.q@med0104                  BP    0/9/32         15.93    lx-amd64      
	gc:license_matlab_r2016b=15
---------------------------------------------------------------------------------
all.q@med0105                  BP    0/24/32        18.28    lx-amd64      
	gc:license_matlab_r2016b=15
---------------------------------------------------------------------------------
all.q@med0106                  BP    0/8/32         24.56    lx-amd64      
	gc:license_matlab_r2016b=15
```

## A Working Example

Get a checkout of our MATLAB example.
Then, look around at the contents of this repository.

```console
$ git clone https://gitlab.com/bihealth/matlab_example.git
$ cat job_script.sh
[...]
$ cat example.m
[...]
```

For submitting the script, you can do the following

```console
$ qsub job_script.sh
```

This will submit a job with one Matlab license requested.
If you were to submit 17 of these jobs, then at least one of them would have to wait until all licenses are free.


!!! warning "Matlab License Server"
    Note that there is a Matlab license server running on the server that will check whether 16 or less Matlab sessions are currently running.
    If a Matlab session is running but this is not made known to the scheduler via `-l license_matlab_r2016b=1` then this can lead to scripts crashing as not enough licenses are available.
    If this happens to you, double-check that you have specified the license requirements correctly and notify hpc-helpdesk@bihealth.de in case of any problems.
    We will try to sort out the situation then.
