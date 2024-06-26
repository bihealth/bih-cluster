# Automated Cleanup of Scratch

The `scratch` space is automatically cleaned up nightly with the following mechanism.

1. Daily snapshots of the `scratch` folder are created and retained for 3 days.
2. Files which were not modified for the last 14 days are removed.
3. Erroneously deleted files can be manually retrieved from the snapshots.

!!! Warning
    We specifically use the `mtime` attribute to determine if files in scratch 
    should be cleaned up. Copying or downloading files to scratch while preserving 
    the original `mtime` might lead to unexpected results.