# General Software Craftmanship

>
Computer software, or simply software, is a generic term that refers to a collection of data or computer instructions that tell the computer how to work, in contrast to the physical hardware from which the system is built, that actually performs the work.<br>
-- [Wikipedia: Software](https://en.wikipedia.org/wiki/Software)
>

As you will most probably never have contact with the HPC system hardware, everything you interact with on the HPC is *software*.
All of your scripts, your configuration files, programs installed by you or administration, and all of your data.

This should also answer the question why you should care about software and why you should try to create and use software of a minimal quality.

>
Software craftsmanship is an approach to software development that emphasizes the coding skills of the software developers themselves.<br>
-- [Wikipedia: Software Craftmanship](https://en.wikipedia.org/wiki/Software_craftsmanship)
>

This Wiki page is not mean to give you an introduction of creating good software but rather collect a (growing) list of easy-to-use and high-impact points to improve software quality.
Also, it provides pointers to resources elsewhere on the internet.

## Use Version Control

Use a version control system for your configuration and your code.
Full stop.
Modern version control systems are Git and Subversion.

- [Official Git Documentation](https://git-scm.com/doc)
- [Github Help](https://help.github.com/)
- [Fix Common Git Problems](https://ohshitgit.com/)

## Do not Share Git/SVN Checkouts for Multiple Users

Every user should have their own Git/Subversion checkout.
Otherwise you are inviting a large number of problems.

## Document Your Code

This includes

- programmer-level documentation in your source code, both inline and per code unit (e.g., function/class)
- top-level documentation, e.g., in README files.

## Document Your Data

Document where you got things from, how to re-download, etc.
E.g., put a README file into each of your data top level directories.

## Use Checksums

Use MD5 or other checksums for your data.
For example, `md5sum` and `hashdeep` are useful utilities for computing and checking them:

- [`md5sum` How-To](https://www.tecmint.com/generate-verify-check-files-md5-checksum-linux/) (tools such as `sha256sum` work the same...)
- [`hashdeep` How-To](http://www.vpsinfo.com/tutorial/file-integrity-with-hashdeep/)

## Use a Workflow Management System

Use some system for managing your workflows.
These systems support you by

- Detect failures and don't continue working with broken data,
- continue where you left off when someting breaks,
- make things more reproducible,
- allow distribution of jobs on the cluster.

[Snakemake](https://snakemake.readthedocs.org) is a popular workflow management system widely used in Bioinformatics.
A minimal approach is using [Makefiles](https://bsmith89.github.io/make-bml/).

## Understand Bash and Shell Exit Codes

If you don't want to use a workflow management system, e.g., for one-step jobs, you should at least understand Bash job management and exit codes.
For example, you can use `if/then/fi` in Bash together with exit codes to:

- Only call a command if the previous command succeded.
- Remove incomplete output files in case of errors.

```bash
if [[ ! -e file.md5 ]]; then
    md5sum file >file.md5 \
    || rm -f file.md5
fi
```

Also, learn about the [inofficial Bash strict mode](http://redsymbol.net/articles/unofficial-bash-strict-mode/).
