# Automated Cleanup of Scratch

The `scratch` space is automatically cleaned up nightly with the following mechanism.

1. Create daily snapshots of the `scratch` folder and retain for 3 days.
2. Remove files which were not modified for the last 14 days.
3. Remove folders which were not modified for the last 14 days and are empty.
4. Erroneously deleted files can be manually retrieved from the snapshots.

!!! Info
    If you want to prevent a folder from being removed, place an empty `.keepdir`
    file inside of it. This will protect only the folder, not its content!

    ```
    $ touch ~/scratch/tmp/.keepdir
    ```

    There are also exception rules in place for the `.cache` folder many users
    symlink to scratch.

!!! Warning
    We specifically use the `mtime` attribute to determine if files in scratch 
    should be cleaned up. Copying or downloading files to scratch while preserving 
    the original `mtime` might lead to unexpected results.
