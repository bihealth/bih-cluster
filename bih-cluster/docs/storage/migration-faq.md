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
$ rsync -ahPX --stats --dry-run $SOURCE $TARGET
```

    !!! warning "Important"
        Please note the importance of the -X flag to keep extended file attributes (ACLs) which
        we might have assigned to you if you are a delegate in charge of moving a project.

2. Remove the `--dry-run` flag to start the actual copying process.
3. Perform a second `rsync` to check if all files were successfully transferred.
   Paranoid users might want to add the `--checksums` flag to `rsync` or use `hashdeep`.
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

## Moving user home and work
1. First copy your home folder while skipping symbolic links.
   This is necessary because the locations of work and scratch changed and we don't want to drag along the outdated links.
   replace the parts in curly brackets with your actual user name and remove the `--dry-run` flag to perform the actual transfer.
   It is important to end paths with a trailing slash (`/`) as this is interpreted by `sync` as “all files in this folder”.
```sh
$ SOURCE=/data/gpfs-1/home/users/{username_c}/
$ TARGET=/data/cephfs-1/home/users/{username_c}/
$ rsync -ahP --stats --no-links --dry-run $SOURCE $TARGET
```
2. Rsync will not follow the symbolic link to your work folder.
   We therefore need to copy contents of your work directory separately.
```sh
$ SOURCE=/data/gpfs-1/work/users/{username_c}/
$ TARGET=/data/cephfs-1/home/users/{username_c}/work/
$ rsync -ahP --stats --dry-run $SOURCE $TARGET
```
3. Perform a second `rsync` per location to check if all files were successfully transferred.
   Paranoid users might want to add the `--checksums` flag to `rsync` or use `hashdeep`.
   Please note the flag `--remove-source-files` which will do exactly as the name suggests,
   but leaves empty directories behind.
   
    !!! Warning
        Check thoroughly that files were actually copied as expected before removing the `--dry-run` flag.
        Use absolute paths to not be confused by symbolic links.
 
```sh
$ rsync -ahX --stats --remove-source-files --dry-run $SOURCE $TARGET
```
4. Check if all files are gone from the SOURCE folder and remove the empty directories:
```sh
$ find $SOURCE -type f | wc -l
0
$ rm -r $SOURCE
```

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

2. install a new version of conda/mamba in your home (or better in `/data/cephfs-1/work/groups/<group>/users/<user>`) and run `source activate /path/to/new/conda/bin/activate`

3. Re-create them after the move:
```sh
$ conda env create -f environment.yml
```

(if you run into errors it might be better to do `conda env export -n $env --no-builds -f $env.yaml`)

!!! Note
    If you already moved your home folder, you can still activate your old environments like this:

    ```sh
    $ conda activate /fast/home/users/your-user/path/to/conda/envs/env-name-here
    ```
## Conda environments - Second option (experimental)
There is a way of moving conda environments that does not involve re-installing. Follow these steps:

1. Copy (cp -r) or rsync the miniconda directory from the fast location (``` /fast/work/users/<username>/miniconda ```) to the cephfs-1 (``` /data/cephfs-1/work/groups/<groupname>/users/<username>/miniconda ```, shorter:  ``` ~/work/miniconda ```)
2. Run ``` conda init ```. It will give you a list of all files conda looks at to find the conda path when initializing. Obviously, since you did not yet move conda, the paths of the files will be still pointing to  /data/gpfs-1 or /fast. These files are however now present also in a copy on the cephfs-1. These files could be for example: ``` /data/cephfs-1/work/groups/cubi/users/<username>/miniconda/condabin/conda ``` or ``` /data/cephfs-1/work/groups/cubi/users/<username>/miniconda/etc/profile.d/conda.sh ``` and multiple others (depending on the output of conda init). On each of these files you now run the command:
    ```sh
    sed -i 's/\/data\/gpfs-1\/work\/users\/<username>\/miniconda/\/data\/cephfs-1\/work\/groups\/<groupname>\/users\/<username>\/miniconda/g' <file_path/file_name>
    ```

    !!! Warning
        Be careful here and check again in the above mentioned file with less command how the old path look like and with what you want to replace them. 
        sed -i 's/<old_path>/<new_path>/g' file (backslash shown above are escape characters)
        Also note that new paths are absolute paths and not relative ones (for example use  ``` /data/cephfs-1/home/users/<username>/work ``` instead of  ``` $HOME/work  ```)

3. Now change your ~/.bashrc file in your home on cephfs-1 ( ``` /data/cephfs-1/home/users/<username>  ```) such that all paths that were pointing to the work dir on fast are now pointing to the work dir on cephfs-1. Also, if not already done, deactivate the automatic activation of the base env. This is how .bashrc should look like:
```sh
case "${SLURMD_NODENAME-${HOSTNAME}}" in
    login-*)
        ;;
    *)
        export PATH=/data/cephfs-1/home/users/<username>/work/miniconda/condabin:$PATH
        ;;
esac


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/cephfs-1/work/groups/cubi/users/<username>/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/cephfs-1/work/groups/cubi/users/<username>/miniconda/etc/profile.d/conda.sh" ]; then
        . "/data/cephfs-1/work/groups/cubi/users/<username>/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/data/cephfs-1/work/groups/cubi/users/<username>/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
export CONDA_AUTO_ACTIVATE_BASE=false
```

4. To really make sure that this will work  remove the old path pointing to miniconda on /fast from your $PATH variable (altough this variable should get initialized everytime you connect to the cluster). This can be done by simply printing the $PATH ( ``` echo $PATH ```), copying the content, deleting the path pointing to fast and copying the modifies content to $PATH again by  ``` export PATH=<paste content>  ```
5. After changing all these files and saving them(!), exit the cluster connection and connect again. As sanity check you can either open the .bashrc and check if all paths are pointing to cephfs-1 or run "echo $PATH" and see if the conda path is now pointing to cephfs-1.
6. Now you can use conda activate env_name as you did before

