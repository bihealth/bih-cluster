# Miscellaneous Tips

This page collects some useful tips.

## Interpreting `core` Files

> A core file or core dump is a file that records the memory image of a running process and its process status (register values etc.). Its primary use is post-mortem debugging of a program that crashed while it ran outside a debugger. A program that crashes automatically produces a core file, unless this feature is disabled by the user.
>
> -- https://sourceware.org/gdb/onlinedocs/gdb/Core-File-Generation.html

**Quickstart**

This topic is too large to be discussed here (see the gdb and Wikipedia links).
If you just want to find out what caused the crash, here is how with `file`:

```bash
$ file core.111506
core.111506: ELF 64-bit LSB core file x86-64, version 1 (SYSV), SVR4-style, from '/fast/users/szyskam_c/work/miniconda/envs/fastq_search-env/bin/python3.7 -m sna'
```

**Also See**

- https://en.wikipedia.org/wiki/Core_dump