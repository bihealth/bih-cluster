# Data Migration Tips and tricks
Please use `hpc-transfer-1` and `hpc-transfer-2` for moving large amounts of files.
This not only leaves the compute notes available for actual computation, but also has no risk of your jobs being killed by Slurm.
You should also use `tmux` to not risk connection loss during long running transfers.

## Moving a project folder
1. Define source and target location and copy contents.
   Please replace the parts in curly brackets with your actual folder names.
   It is important to end paths with a trailing slash (`/`) as this is interpreted by `sync` as “all files in this folder”.
```sh
$ SOURCE=/data/gpfs-1/work/projects/{my_project}/
$ TARGET=/data/cephfs-2/unmirrored/projects/{my-project}/
$ rsync -ahP --stats --dry-run $SOURCE $TARGET
```

2. Remove the `--dry-run` flag to start the actual copying process.

    !!! warning "Important"
        File ownership information will be lost during this process.
        This is due to non-root users
        [not being allowed](https://serverfault.com/questions/755753/preserve-ownership-with-rsync-without-root)
        to change ownership of arbitrary files.
        If this is a problem for you, please contact our admins again after completing this step.

3. Perform a second `rsync` to check if all files were successfully transferred.
   Paranoid users might want to add the `--checksum` flag to `rsync` or use `hashdeep`.
   Please note the flag `--remove-source-files` which will do exactly as the name suggests,
   but leaves empty directories behind.
```sh
$ rsync -ahX --stats --remove-source-files --dry-run $SOURCE $TARGET
```
4. Again, remove the `--dry-run` flag to start the actual deletion.
5. Check if all files are gone from the SOURCE folder and remove the empty directories:
```sh
$ find $SOURCE -type f | wc -l
0
$ rm -r $SOURCE
```

!!! Warning 
    When defining your SOURCE location, do not use the `*` wildcard character.
    It will not match hidden (dot) files and leave them behind.
    Its better to use a trailing slash which matches “All files in this folder”.

## Moving user work folders
1. All files within your own work directory can be transferred as follows.
   Please replace parts in curly braces with your cluster user name.
```sh
$ SOURCE=/data/gpfs-1/work/users/{username}/
$ TARGET=/data/cephfs-1/home/users/{username}/work/
$ rsync -ahP --stats --dry-run $SOURCE $TARGET
```

    !!! Note
        The `--dry-run` flag lets you check that rsync is working as expected without copying any files.
        Remove it to start the actual transfer.
    
2. Perform a second `rsync` to check if all files were successfully transferred.
   Paranoid users might want to add the `--checksums` flag or use `hashdeep`.
   Please note the flag `--remove-source-files` which will do exactly as the name suggests,
   but leaves empty directories behind.
```sh
$ rsync -ahP --stats --remove-source-files --dry-run $SOURCE $TARGET
```
4. Check if all files are gone from the SOURCE folder:
```sh
$ find $SOURCE -type f | wc -l
0
```

### Conda environments
Conda installations tend not to react well to moving their main folder from its original location.
There are numerous ways around this problem which are described [here](https://www.anaconda.com/blog/moving-conda-environments).

A simple solution we can recommend is this:

1. Export all environments prior to the move with this bash script:
```sh
#!/bin/bash
for env in $(ls $(conda info --base)/envs/)
do
    conda env export -n $env -f $env.yml
done
```
   If you run into errors it might be better to use `conda env export -n $env --no-builds -f $env.yaml`.

2. Install a fresh version of conda or mamba in your new work folder.
   Don't forget to turn off automatic base environment activation for less delay during login and reduced strain on the login nodes.
```sh
$ conda init
$ conda config --set auto_activate_base false
```

3. Re-create your old environments from the yaml files:
```sh
$ conda env create -f {environment.yml}
```

!!! Note "Tip"
    If we already moved your home folder, you can still activate your old environments like this:

    ```sh
    $ conda activate /fast/home/users/your-user/path/to/conda/envs/{env-name-here}
    ```
