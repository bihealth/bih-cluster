# Slurm Command Format Strings

In the sections [Slurm Quickstart](quickstart.md) and [Slurm Cheat Sheet](cheat-sheet.md), we have seen that `sinfo` and `squeue` allow for the compact display *partitions/nodes* and *node* information.
In contrast, `scontrol show job <id>` and `scontrol show partition <id>` and `scontrol show node <id>` show comprehensive information that quickly gets hard to comprehend for multiple entries.

Now you might ask: *is there anything in between?*
And: *yes, there is*.

You can tune the output of `sinfo` and `squeue` using parameters, in particular by providing **format strings**.
All of this is described in the man pages of the commands that you can display with `man sinfo` and `man squeue` on the cluster.

## Tuning `sinfo` Output

Notable arguments of `sinfo` are:

- `-N, --Node` -- uncompress the usual lines and display one line per node and partition.
- `-s, --summarize` -- compress the node state, more compact display.
- `-R, --list-reasons` -- for nodes that are not *up*, display reason string provided by admin.
- `-o <fmt>, --format=<fmt>` -- use format string for display.

The most interesting argument is `-o/--format`.
The man page lists the following values that are used when using other arguments.
In other words, many of the display modifications could also be applied with `-o/--format`.

```
default        "%#P %.5a %.10l %.6D %.6t %N"
--summarize    "%#P %.5a %.10l %.16F  %N"
--long         "%#P %.5a %.10l %.10s %.4r %.8h %.10g %.6D %.11T %N"
--Node         "%#N %.6D %#P %6t"
--long --Node  "%#N %.6D %#P %.11T %.4c %.8z %.6m %.8d %.6w %.8f %20E"
--list-reasons "%20E %9u %19H %N"
--long --list-reasons
                "%20E %12U %19H %6t %N"
```

The best way to learn more about this is to play around with `sinfo -o`, starting out with one of the format strings above.
Details about the format strings are described in `man sinfo`.
Some remarks here:

- `%<num><char>` displays the value represented by `<char>` padded with spaces to the right such that a width of `<num>` is reached,
- `%.<num><char>` displays the value represented by `<char>` padded with spaces to the **left** such that a width of `<num>` is reached, and
- `%#<char>` displays the value represented by `<char>` padded with spaces to the max length of the value represented by `<char>` (this is a "virtual" value, used internally only, you cannot use this and you will have to place an integer here).

For example, to create a grouped display with reasons for being down use:

```bash
hpc-login-1:~$ sinfo -o "%10P %.5a %.10l %.16F  %40N %E"
PARTITION  AVAIL  TIMELIMIT   NODES(A/I/O/T)  NODELIST                                 REASON
debug*        up    8:00:00        0/0/16/16  med[0703-0710,0740-0742,0744-0745,0749,0 bogus node
debug*        up    8:00:00      18/98/0/116  med[0104-0124,0127,0133-0148,0151-0164,0 none
medium        up 7-00:00:00        0/0/16/16  med[0703-0710,0740-0742,0744-0745,0749,0 bogus node
medium        up 7-00:00:00      18/98/0/116  med[0104-0124,0127,0133-0148,0151-0164,0 none
long          up 28-00:00:0        0/0/16/16  med[0703-0710,0740-0742,0744-0745,0749,0 bogus node
long          up 28-00:00:0      18/98/0/116  med[0104-0124,0127,0133-0148,0151-0164,0 none
critical      up 7-00:00:00        0/0/16/16  med[0703-0710,0740-0742,0744-0745,0749,0 bogus node
critical      up 7-00:00:00      18/98/0/116  med[0104-0124,0127,0133-0148,0151-0164,0 none
highmem       up 14-00:00:0          0/4/0/4  med[0401-0404]                           none
gpu           up 14-00:00:0          3/1/0/4  med[0301-0304]                           none
```

## Tuning `squeue` Output

The standard squeue output might yield the following

```bash
hpc-login-1:~$ squeue | head
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
              3149    medium variant_ holtgrem PD       0:00      1 (Dependency)
              1177    medium     bash jweiner_  R 6-03:32:41      1 med0127
              1192    medium     bash jweiner_  R 5-12:48:57      1 med0127
              1210       gpu     bash hilberta  R 2-16:10:51      1 med0304
              1213      long     bash schubacm  R 2-15:22:44      1 med0127
              2401       gpu     bash ramkem_c  R 2-10:55:10      1 med0303
              3063      long     bash schubacm  R 1-09:52:54      1 med0127
              3066      long     bash schubacm  R 1-09:52:04      1 med0127
              3147    medium ngs_mapp holtgrem  R 1-03:13:42      1 med0148
```

Looking at `man squeue`, we learn that the default format strings are:

```
default        "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
-l, --long     "%.18i %.9P %.8j %.8u %.8T %.10M %.9l %.6D %R"
-s, --steps    "%.15i %.8j %.9P %.8u %.9M %N"
```

This looks a bit wasteful.
Let's cut down on the padding of the job ID and expand on the job name and remove some right paddings.

```bash
hpc-login-1:~$ squeue -o "%.6i %9P %30j %.10u %.2t %.10M %.6D %R %b" | head
 JOBID PARTITION NAME                                 USER ST       TIME  NODES NODELIST(REASON)
  3149 medium    variant_calling                holtgrem_c PD       0:00      1 (Dependency)
  1177 medium    bash                            jweiner_m  R 6-03:35:55      1 med0127
  1192 medium    bash                            jweiner_m  R 5-12:52:11      1 med0127
  1210 gpu       bash                           hilberta_c  R 2-16:14:05      1 med0304
  1213 long      bash                           schubacm_c  R 2-15:25:58      1 med0127
  2401 gpu       bash                             ramkem_c  R 2-10:58:24      1 med0303
  3063 long      bash                           schubacm_c  R 1-09:56:08      1 med0127
  3066 long      bash                           schubacm_c  R 1-09:55:18      1 med0127
  3147 medium    ngs_mapping                    holtgrem_c  R 1-03:16:56      1 med0148
```

### Displaying Resources

Now display how many of our internal projects still exist.

```bash
hpc-login-1:~$ squeue -o "%.6i %9P %30j %.10u %.2t %.10M %.6D %10R %s" | head
```

The next steps are (TODO):

- setup of certificate for containers
- opening firewall apropriately
- integrate with openmpi documentation