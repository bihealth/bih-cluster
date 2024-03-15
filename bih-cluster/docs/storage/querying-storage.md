# Querying Storage Quotas

!!! info "Outdated"

    This document is only valid for the old, third-generation file system and will be removed soon.
    Quotas of our new CephFS storage are communicated via the [HPC Access](https://hpc-access.cubi.bihealth.org/) web portal.

As described elsewhere, all data in your user, group, and project volumes is subject to quotas.
This page quickly shows how to query for the current usage of data volume and file counts for your user, group, and projects.

## Query for User Data and File Usage

The file `/etc/bashrc.gpfs-quota` contains some Bash functions that you can use for querying the quota usage.
This file is automatically sourced in all of your Bash sessions.

For querying your user's data and file usage, enter the following command:

```
# bih-gpfs-quota-user holtgrem_c
```

You will get a report as follows.
As soon as usage reaches 90%, the data/file usage will be highlighted in yellow.
If you pass 99%, the data/file usage will be highlighted in red.

```
=================================
Quota Report for: user holtgrem_c
=================================

                           DATA       quota       GR- FILES      quota       GR-
ENTITY  NAME       FSET    USED       SOFT  HARD  ACE USED       SOFT  HARD  ACE
------- ---------- ------- ----- ---- ----- ----- --- ----- ---- ----- ----- ---
users   holtgrem_c home     103M  10%  1.0G  1.5G   -  2.5k  25%   10k   12k   -
users   holtgrem_c work     639G  62%  1.0T  1.1T   -  1.0M  52%  2.0M  2.2M   -
users   holtgrem_c scratch   42G   0%  200T  220T   -  207k 0.1%  200M  220M   -
[...]
```

## Query for Group Data and File Usage

```
# bih-gpfs-report-quota group ag_someag
=================================
Quota Report for: group ag_someag
=================================

                           DATA       quota       GR- FILES      quota       GR-
ENTITY  NAME       FSET    USED       SOFT  HARD  ACE USED       SOFT  HARD  ACE
------- ---------- ------- ----- ---- ----- ----- --- ----- ---- ----- ----- ---
groups  ag_someag  home        0   0%  1.0G  1.5G   -     4   0%   10k   12k   -
groups  ag_someag  work     349G  34%  1.0T  1.5T   -   302   0%  2.0M  2.2M   -
groups  ag_someag  scratch     0   0%  200T  220T   -     1   0%  200M  220M   -

[...]
```

## Query for Project Data and File Usage

```
# bih-gpfs-report-quota project someproj
==================================
Quota Report for: project someproj
==================================

                           DATA       quota       GR- FILES      quota       GR-
ENTITY  NAME       FSET    USED       SOFT  HARD  ACE USED       SOFT  HARD  ACE
------- ---------- ------- ----- ---- ----- ----- --- ----- ---- ----- ----- ---
groups  someproj   home        0   0%  1.0G  1.5G   -     4   0%   10k   12k   -
groups  someproj   work     349G  34%  1.0T  1.5T   -   302   0%  2.0M  2.2M   -
groups  someproj   scratch     0   0%  200T  220T   -     1   0%  200M  220M   -

[...]
```
