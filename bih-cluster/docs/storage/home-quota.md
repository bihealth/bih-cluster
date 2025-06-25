# Keeping your home folder clean
We set quite restrictive quotas for user homes, but in exchange you get file system [snapshots and mirroring](./storage-locations.md#snapshots-and-mirroring).
Your home folder should therefore only be used for scripts, your user config, and other small files.
Everything else should be stored in the `work` or `scratch` subdirectories, which effectively link to your group's shared storage space.
This document describes some common pitfalls and how to circumvent them.

!!! hint
    The tilde character (`~`) is shorthand for your home directory.

## Code libraries and other big folders
Various programs are used to depositing large folders in a user's home and can quickly use up your allotted storage quota.
These include:

- Python: `~/.local/lib/python*`
- conda: Location chosen by the user, but also `~/.conda`
- R: `~/R/x86_64-pc-linux-gnu-library`
- [HPC portal](../ondemand/overview.md): `~/ondemand`

Please note that directories whose name is starting with a dot are not shown by the normal `ls` command, but require the `ls -a` flag. You can search your home folder for large directories like so:
```bash
$ du -shc ~/.* ~/* --exclude=.. --exclude=.
```

You should move these locations to your `work` folder and create symbolic links in their place.
Conda installations should be installed in `work` from the very beginning as they do not react well to being moved around.

Here is an example for the `.local` folder.

```bash
$ mv ~/.local ~/work/.local
$ ln -s ~/work/.local ~/.local
```

## Temporary Files
Another usual culprit is the hidden `.cache` directory which contains temporary files.
This folder can be moved to the `scratch` volume in a similar manner as described above.

```bash
$ mv ~/.cache ~/scratch/.cache
$ ln -s ~/scratch/.cache ~/.cache
```

!!! warning "Important"
    Files placed in your `scratch` directory will be [automatically removed](./scratch-cleanup.md) after 2 weeks.
    Do not place any valuable files in there.
