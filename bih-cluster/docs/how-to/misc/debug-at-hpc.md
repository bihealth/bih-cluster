# How-To: Debug Software on HPC Systems

!!! info "Please Contribute!"

    This guide is far from complete.
    Please feel free to contribute, e.g., refer to [How-To: Contribute to this Document](contribute.md).

Please make sure that you have read [How-To: Debug Software](debug-software.md) as a general primer.

As debugging is hard enough already, it makes one wonder how to do this on the HPC system in batch mode.
Here is a list of pointers.

## Attempt 1: Run it interactively!

First of all, you can of course get an interactive session using `srun --pty bash -i` and then run your program interactively.
Make sure to allocate appropriate memory and cores for your purpose.
You might also want to first start a `screen` or `tmux` session on the login node such that network interruptions to the login node don't harm your hard debugging work!

Does the program work correctly if you do this?
If yes, and it only fails when run in batch mode, consider the following behaviour of the scheduler.

The scheduler takes your resource requirements and tries to find a free slot.
Once it has found a free slot, it will attempt to run the program.
This mainly differs in running it interactively in standard input, output, and error streams.

- By default, stdin is connected to `/dev/null` such that no input is read.
  You can change this with the `--input=` flag to specify a file.
- By default stdout and stderr are joint and written to the file specified as `--output=`.
  You can use certain wildcards to make the output (but also the input files) depend on things like the job ID or job name.
- Please note that the **directory name to the output file** but exist before the job is launched.
  It is **not** sufficient to `mkdir` it in the job script itself.

Please refer to the [sbatch documentation](https://slurm.schedmd.com/sbatch.html) for details.

If your program fails without leaving any log file or any other trace, make sure that the path to the output file exists.
To the best of the author's knowledge, there is no way to tell apart a crash because this does not exist and a program failure (except maybe for the running time of 0 seconds and memory usage of 0 bytes).

## Attempt 2: Inspect the logs

Do you see any exception in your log files?
If not, continue.

If your job is canceled by `scancel` or stopped because it exhausted it maximal running time or allocated resources then you will find a note in the last line of your error output log (usually folded into the standard output).
Please note that if the previous output line did not include a line ending, the message might be at the very end of the last line.

The message will look similar to:

```
slurmstepd: error: *** JOB <your job id> ON med0xxx CANCELLED AT 2020-09-02T21:01:12 DUE TO TIME LIMIT ***
```

## Attempt 3: Increase logging/printing

Ideally, you can add one or more `--verbose`/`-v` flags to your program to increase verbosity.
See how far your program gets, see where it fails.
This attempt will be greatly helped by reproducible running on a minimal working example.

## Attempt 4: Use `sattach`

You can use [`sattach`](../../slurm/commands-sattach.md) for attaching your terminal to your running job.
This way, you can perform an interactive inspection of the commands.

You can combine this with one of the next attempst of using debuggers to e.g., get an `pdb` debugger at an important position of your program.
However, please note that `pdb` and `ipdb` will stop the program's execution if the standard input stream is at end of file (which `/dev/null` is and this is used by default in `sbatch` jobs).

## Attempt 5: Inspect Program Activity

Log into the node that your program runs on either using `srun --pty --nodelist=NODE` or using `ssh`.
Please note that you should never perform computational intensive things when logging into the node directly.
You can then use all activity inspection tips from [How-To: Debug Software](debug-software.md).

## Attempt 6: Use Debuggers

After having logged into the node running your program, you can of course also attach to the program with `gdb -p PID` or `cgdb -p PID`.

## Don't Despair

Here are some final remarks:

- Don't despair!
- The longer you search for the problem, the more fundamental it is.
  Chances are that you are just overlooking something obvious which is actually easy to fix.
- Keep old log files!
- Really, really, make sure that your program runs deterministically.
  You will save yourself a world of pain.
