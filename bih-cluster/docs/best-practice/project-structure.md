# Project File System Structure

This Wiki page dscribes best pratices for managing your Bioinformatics projects on the file system.

## General Aims

Mostly, you can separate the files in your projects/pipelines into one of the following categories:

1. scripts (and their documentation)
2. configuration
3. data

Ideally, scripts and documentation are independent of a given project and can be separated from the rest.
Configuration is project-dependent and small and mostly does not contain any sensitive information (such as genotypes that allows for reidentification of donors).
In most cases, data might be large and is either also stored elsewhere or together with scripts and configuration can be regenerated easily.

!!! warning "There is no backup of `work` and `scratch`"
    The cluster GPFS file system `/fast` is not appropriate for keeping around single "master" copies of data.
    You should have a backup and archival strategy for your valuable "master" copy data.

## Best Practices

### Scripts

- Your scripts should go into version control, e.g., a Git repository.
- Your scripts should be driven by command line parameters and/or configuration such that no paths etc. are hard-coded.
  If for a second data set, you need to make a copy of your scripts and adjust some variables, e.g., at the top, you're doing something in a suboptimal fashion.
  Rather, get these values from the command line or a configuration file and only store (sensible) defaults in your script where appropriate.
- Thus, ideally your scripts are not project-specific.

### Configuration

- Your configuration usually **is** project-specific.
- Your configuration should also go into version contro, e.g., a Git repository.

In addition, you might need project-specific "wrapper" scripts that just call your project-independent script with the correct paths for your project.
These scripts rather fall into the "configuration" category and should then live together with your configuration.

### Data

- Your data should go into a location separate from your scripts and configuration.
- Ideally, the raw input data is separated from the work and output files such that you can make these files and directories read-only and don't accidentally damage these files.
- **You really should keep temporary files in a temporary directory, set the environment variable `TMPDIR` appropriately and automatically clean them up (see [Useful Tips: Temporary Files](Manual-Useful-Tips-Temporary-Files))**

## Best Practices in Practice

But how can we put this into practice?
Below, we give some examples of how to do this.
Note that for simplicity's sake we put all scripts and configuration into one directory/repository contrary to the best practices above.
This is for educational purposes only and you should strive for reuseable scripts where it makes sense and separate scripts and configuration.

We will limit this to simple Bash scripts for education's purposes.
You should be able to easily adapt this to your use cases.

Thus, the aim is to separate the data from the non-data part of the project such that we can put the non-data part of the project into a separate location and under version control.
We call the location for non-data part of the project the **home** location of your project and the location for the data part of the project the **work** location of your project.

Overall, we have three options:

1. Your processes are run in the home location and the sub directories used for execution are links into the work location using symlinks.
2. Your processes are run in the work location and
    1. the scripts are linked into the work location using symlinks, OR
    2. the scripts are called from the home location, maybe through project-specific wrapper scripts.

## Example: Link config/scripts into work location (Option 1)

Creating the work directory and copy the input files into `work/input`.

```bash
$ mkdir -p project/work/input
$ cp /fast/projects/cubit/tutorial/input/* project/work/input
```

Creating the home space.
We initialize a Git repository, properly configure the `.gitignore` file and add a `README.md` file.

```bash
$ mkdir -p project/home
$ cd project/home
$ cat <<EOF >.gitignore
*~
.*.sw?
EOF
$ cat <<EOF >README.md
# Example Project

This is an example project with config/scripts linked into work location.
EOF
$ git init
$ git add .gitignore README.md
$ git commit -m 'Initial project#
```

We then create the a simple script for executing the mapping step and a configuration file that gives the path to the index and list of samples to process.

```bash
$ mkdir scripts
$ cat <<"EOF" >scripts/run-mapping.sh
#!/bin/bash

# Unofficial Bash script mode, see:
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail

# Get directory to bash file, see
# https://stackoverflow.com/a/4774063/84349
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

# Helper function to print help to stderr.
help()
{
  >&2 echo "Run Mapping Step"
  >&2 echo ""
  >&2 echo "run-mapping.sh [-c config.sh] [-h]"
}

# Parse command line arguments into bash variables.
CONFIG=
while getopts "hs:" arg; do
  case $arg in
    h)
      help()
      exit
      ;;
    s)
      CONFIG=$OPTARG
      ;;
  esac
done

# Print the executed commands.
set -x

# Load default configuration, then load configuration file if any was given.
source $SCRIPTPATH/../config/default-config.sh
if [[ -z "$CONFIG" ]]; then
    source $CONFIG
fi

# Create output directory.
mkdir -p output

# Actually perform the mapping.  This assumes that you have
# made the bwa and samtools commands available, e.g., using conda.
for sample in $SAMPLES; do
    bwa mem \
        $BWA_INDEX \
        input/${sample}_R1.fq.gz \
        input/${sample}_R2.fq.gz \
    | samtools sort \
        -o output/${sample}.bam \
        /dev/stdin
done

EOF
$ chmod +x scripts/run-mapping.sh
$ mkdir -p config
$ cat <<"EOF" >config/default-config.sh
BWA_INDEX=/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa
SAMPLES=
EOF
$ cat <<"EOF" >config/project-config.sh
$ BWA_INDEX comes from default configuration already
SAMPLES=test
EOF
```

This concludes the basic project setup.
Now, to the symlinks:

```bash
$ cd ../work
$ ln -s ../home/scripts ../home/config .
```

And, to the execution...

```bash
$ ./scripts/run-mapping -c config/project-config.sh
[...]
```

## Example: Link Data Into Home (Option 2.1).

We can reuse the project up to the statement "This concludes the basic project setup" in the example for option 1.

Then, we can do the following:

```
$ cd ../work
$ mkdir -p output

$ cd ../home
$ cat <<"EOF" >>.gitignore

# Ignore all data
input/
work/
output/
EOF
$ git add .gitignore
$ git commit -m 'Ignoring data file in .gitignore'
$ ln -s ../work ../output .
```

And we can execute everything in the home directory.

```bash
$ ./scripts/run-mapping -c config/project-config.sh
[...]
```

## Example: Wrapper Scripts in Home (Option 2.2)

Again, we can reuse the project up to the statement "This concludes the basic project setup" in the example for option 1.

Then, we do the following:

```
$ cd ../work
$ cat <<"EOF" >do-run-mapping.sh
#!/bin/bash

../home/scripts/run-mapping.sh \
    -c ../home/config/project-config.sh
EOF
$ chmod +x do-run-mapping.sh
```

Note that the the `do-run.sh` script could also go into the project-specific Git repository and be linked into the work directory.

Finally, we can run our pipeline:

```bash
$ cd ../work
$ ./do-run-mapping.sh
[...]
```
