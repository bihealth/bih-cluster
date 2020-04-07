# Migrating to Slurm

This page describes (an updated) list of steps that are commonly needed to migrate from SGE to Slurm.

## Summary

- Use `#SBATCH` instead of `#$` in your job scripts for resource definitions.
- Update your resource requirements parameters from SGE to Slurm using the [Slurm Rosetta Stone](rosetta-stone.md).
- Unset environment variable `DRMAA_LIBRARY_PATH` and remove from your `~/.bashrc`, job scripts, etc.

These steps are explained in detail below.

## Using `#SBATCH`

Slurm uses `#SBATCH` as the line marker at the head of your script file.

So instead of

```bash
#$ -N my_job
```

You now use

```bash
#SBATCH --job-name my_job
```

## Use Slurm Parameters

Of course, the Slurm scheduler parameters differ from the SGE ones.
So you also have to use Slurm syntax instead of SGE syntax (e.g., swich from `-N` to `-J/--job-name` as show in the previous section).
The page [Slurm Rosetta Stone](rosetta-stone.md) provides you with a translation table.

## Getting Rid of `DRMAA_LIBRARY_PATH`

For SGE, it was necessary to define the environment variable `DRMAA_LIBRARY_PATH` and the old document told you to put this into your `~/.bashrc` file.
You should remove the setting from your file.

To check the results, you can use `echo $DRMAA_LIBRARY_PATH` in your shell.

You can use `unset DRMAA_LIBRARY_PATH` to remove the environment variable from your current shell session.
