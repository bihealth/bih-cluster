# Data Migration Tips and tricks
Please use `hpc-transfer-1` and `hpc-transfer-2` for moving large amounts of files.
This not only leaves the compute notes available for actual computation, but also has no risk of your jobs being killed by Slurm.
You should also use `tmux` to not risk connection loss during long running transfers.

## Useful commands

1. Define source and target location and copy contents.
```sh
$ SOURCE=/data/gpfs-1/work/projects/my_project/
$ TARGET=/data/cephfs-2/unmirrored/projects/my-project/
$ rsync -ah --stats --progress --dry-run $SOURCE $TARGET
```

2. Remove the `--dry-run` flag to start the actual copying process.
3. Perform a second `rsync` to check if all files were successfully transferred.
   Paranoid users might want to add the `--checksums` flag to `rsync` or use `hashdeep`.
   Please note the flag `--remove-source-files` which will do exactly as the name suggests,
   but leaves empty directories behind.
```sh
$ rsync -ah --stats --remove-source-files --dry-run $SOURCE $TARGET
```
4. Again, remove the `--dry-run` flag to start the actual deletion.
5. Check if all files are gone from the SOURCE folder and remove the empty directories:
```sh
$ find $SOURCE -type f | wc -l
$ rm -r $SOURCE
```

!!! Warning 
    When defining your source location, do not use the `*` wildcard character.
    It will not match hidden (dot) files and leave them behind.

## Conda environments
Conda environment tend to not react well when the folder they are stored in is moved from its original location.
There are numerous ways to move the state of your environments, which are described [here](https://www.anaconda.com/blog/moving-conda-environments).

A simple way we can recommend is this:

1. Export all environments prior to the move.
```sh
#!/bin/bash
for env in $(ls .miniforge/envs/)
do
    conda env export -n $env -f $env.yml
done
```

2. Re-create them after the move:
```sh
$ conda env create -f environment.yml
```

!!! Note
    If you already moved your home folder, you can still activate your old environments like this:

    ```sh
    $ conda activate /fast/home/users/your-user/path/to/conda/envs/env-name-here
    ```
