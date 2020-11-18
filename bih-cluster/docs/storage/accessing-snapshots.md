# Accessing Snapshots and Backups

!!! info "HPC 4 Research Only"

    The following information is only valid for HPC 4 Research.

By now you have probably read that your home directory has strict quotas in place.
You pay this price in usability by the fact that snapshots and backups exist.

## Snapshot

Every night, a snapshot is created of all user, group, and project home directories.
The snapshots are placed in the following locations.

- `/fast/home/.snapshots/$SNAPSHOT/users/$USER`
- `/fast/home/.snapshots/$SNAPSHOT/groups/$GROUP`
- `/fast/home/.snapshots/$SNAPSHOT/projects/$PROJECT`

The snapshot name `$SNAPSHOT` simply is the 3-letter abbreviation of the week day (i.e., Mon, Tue, Thu, Fri, Sat, Sun).

The snapshots contains the state of the home directories at the time that they were made.
This also includes the permissions.
That is you can simply retrieve the state of your home directory of last Tuesday (if the directory already existed back then) at:

- `/fast/home/.snapshots/Tue/users/$USER`

e.g.,

```bash
$ ls /fast/home/.snapshots/Tue/users/$USER
...
```

## Backups

There are very few cases where backups will help you more than snapshots.
As with snapshots, there is a 7-day rotation for backups.
The backups fully reflect the snapshots and thus everything that is in the backups is also in the snapshots.

The only time that you will need backups is when the GPFS and its snapshots are damaged.
This protects the files in your home directory against technical errors and other catastrophes such as fire or water damage of the data center.
