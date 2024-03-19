# First Steps: Episode 3

|Episode|Topic|
|---|---|
| 0 | [How can I install the tools?](episode-0.md) |
| 1 | [How can I use the static data?](episode-1.md) |
| 2 | [How can I distribute my jobs on the cluster (Slurm)?](episode-2.md) |
| **3** | **How can I organize my jobs with Snakemake?** |
| 4 | [How can I combine Snakemake and Slurm?](episode-4.md) |

In this episode we will discuss how we can parallelize steps in a pipeline that are not dependent on each other.
In the last episode we saw a case (the variant calling) that could have been potentially parallelized.

We will take care of that today. Please note that we are not going to use the `sbatch` command we learned
earlier. Thus, this tutorial will run on the same node where you execute the script. We will introduce you to
Snakemake, a tool with which we can model dependencies and run things in parallel. In the next
tutorial we will learn how to submit the jobs with `sbatch` and Snakemake combined.

For those who know `make` already, Snakemake will be familiar. You can think of Snakemake being a bunch
of dedicated bash scripts that you can make dependent on each other. Snakemake will start the next
script when a previous one finishes, and potentially it will run things in parallel if the dependencies allow.

Snakemake can get confusing, especially if the project gets big. This tutorial will only cover the very
basics of this powerful tool. For more, we highly recommend digging into the Snakemake documentation:

* https://snakemake.readthedocs.io/en/stable/
* http://slides.com/johanneskoester/deck-1#/

Every Snakemake run requires a `Snakefile` file. Create a new folder inside your tutorial folder and
copy the skeleton:

```terminal
(first-steps) $ mkdir -p /fast/users/${USER}/work/tutorial/episode3
(first-steps) $ pushd /fast/users/${USER}/work/tutorial/episode3
(first-steps) $ cp /data/cephfs-1/work/projects/cubit/tutorial/skeletons/Snakefile .
(first-steps) $ chmod u+w Snakefile
```

Your `Snakefile` should look as follows:

```python
rule all:
    input:
        'snps/test.vcf',
        'structural_variants/test.vcf'

rule alignment:
    input:
        '/data/cephfs-1/work/projects/cubit/tutorial/input/test_R1.fq.gz',
        '/data/cephfs-1/work/projects/cubit/tutorial/input/test_R2.fq.gz',
    output:
        bam='alignment/test.bam',
        bai='alignment/test.bam.bai',
    shell:
        r"""
        export TMPDIR=/fast/users/${{USER}}/scratch/tmp
        mkdir -p ${{TMPDIR}}

        BWAREF=/data/cephfs-1/work/projects/cubit/current/static_data/precomputed/BWA/0.7.17/GRCh37/g1k_phase1/human_g1k_v37.fasta

        bwa mem -t 8 \
            -R "@RG\tID:FLOWCELL.LANE\tPL:ILLUMINA\tLB:test\tSM:PA01" \
            ${{BWAREF}} \
            {input} \
        | samtools view -b \
        | samtools sort -O BAM -T ${{TMPDIR}} -o {output.bam}

        samtools index {output.bam}
        """

rule structural_variants:
    input:
        'alignment/test.bam'
    output:
        'structural_variants/test.vcf'
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        delly call -o {output} -g ${{REF}} {input}
        """

rule snps:
    input:
        'alignment/test.bam'
    output:
        'snps/test.vcf'
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        gatk HaplotypeCaller \
            -R ${{REF}} \
            -I {input} \
            -ploidy 2 \
            -O {output}
        """
```

Let me explain. The content resembles the same steps we took in the previous tutorials.
Although every step has its own rule (alignment, snp calling, structural variant calling), we could
instead have written everything in one rule. It is up to you to design your rules! Note that
the rule names are arbitrary and not mentioned anywhere else in the file.

But there is one primary rule: the `rule all`. This is the kickoff rule that makes everything run.

As you might have noticed, every rule has three main parameters: `input`, `output` and `shell`.
`input` defines the files that are going into the rule, `output` those that are produced
when executing the rule, and `shell` is the bash script that processes `input` to produce
`output`.

Rule `all` does not have any `output` or `shell`, it uses `input` to start the chain of rules.
Note that the input files of this rule are the output files of rule `snps` and
`structural_variants`. The input of those rules is the output of rule `alignment`. This
is how Snakemake processes the rules: It looks for rule `all` (or a rule that just has `input`
files) and figures out how it can create the required input files with other rules by looking at their
`output` files (the `input` files of one rule must be the `output` files of another rule).
In our case it traces the workflow back to rule `snps` and `structural_variants` as they
have the matching output files. They depend in return on the alignment,
so the `alignment` rule must be executed, and this is the first thing that will be done by Snakemake.

There are also some peculiarities about Snakemake:

* You can name files in `input` or `output` as is done in rule `alignment` with the output files.
* You can access the `input` and `output` files in the script by writing `{input}` or `{output}`.
    - If they are not named, they will be concatenated, separated by white space
    - If they are named, access them with their name, e.g., `{output.bam}`
    - Curly braces must be escaped with curly braces, e.g., for bash variables: `${{VAR}}` instead of `${VAR}` but not Snakemake internal variables like `{input}` or `{output}`
* In the rule `structural_variants` we cheat a bit because delly does not produce output files if it can't find variants.
    - We do this by `touching` (i.e., creating) the required output file. Snakemake has a function for doing so (call `touch()` on the filename).
* Intermediate folders in the path to output files are always created if they don't exist.
* Because Snakemake is Python based, you can write your own functions for it to use, e.g. for creating file names automatically.

But Snakemake can do more. It is able to parse the paths of the output files and set wildcards if you want.
For this your input (and output) file names have to follow a parsable scheme. In our case they do!
Our FASTQ files, our only initial input files, start with `test`. The output of the alignment as well
as the variant calling is also prefixed `test`. We now can modify the Snakemake file accordingly, by exchanging every
occurrence of `test` in each `input` or `output` field with `{id}` (note that you could also give a different
name for your variable). Only the input rule should not be touched, otherwise Snakemake would not know
which value this variable should have. Your `Snakefile` should look now like this:

```python
rule all:
    input:
        'snps/test.vcf',
        'structural_variants/test.vcf'

rule alignment:
    input:
        '/data/cephfs-1/work/projects/cubit/tutorial/input/{id}_R1.fq.gz',
        '/data/cephfs-1/work/projects/cubit/tutorial/input/{id}_R2.fq.gz',
    output:
        bam='alignment/{id}.bam',
        bai='alignment/{id}.bam.bai',
    shell:
        r"""
        export TMPDIR=/fast/users/${{USER}}/scratch/tmp
        mkdir -p ${{TMPDIR}}

        BWAREF=/data/cephfs-1/work/projects/cubit/current/static_data/precomputed/BWA/0.7.17/GRCh37/g1k_phase1/human_g1k_v37.fasta

        bwa mem -t 8 \
            -R "@RG\tID:FLOWCELL.LANE\tPL:ILLUMINA\tLB:test\tSM:PA01" \
            ${{BWAREF}} \
            {input} \
        | samtools view -b \
        | samtools sort -O BAM -T ${{TMPDIR}} -o {output.bam}

        samtools index {output.bam}
        """

rule structural_variants:
    input:
        'alignment/{id}.bam'
    output:
        'structural_variants/{id}.vcf'
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        delly call -o {output} -g ${{REF}} {input}
        """

rule snps:
    input:
        'alignment/{id}.bam'
    output:
        'snps/{id}.vcf'
    shell:
        r"""
        REF=/data/cephfs-1/work/projects/cubit/current/static_data/reference/GRCh37/g1k_phase1/human_g1k_v37.fasta

        gatk HaplotypeCaller \
            -R ${{REF}} \
            -I {input} \
            -ploidy 2 \
            -O {output}
        """
```

Before we finally run this, we can make a dry run. Snakemake will show you what it would do:

```terminal
(first-steps) $ snakemake -n
```

If everything looks green, you can run it for real. We provide it two cores
to allow two single-threaded jobs to be run simultaneously:

```terminal
(first-steps) $ snakemake -j 2
```
