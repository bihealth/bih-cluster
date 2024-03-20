# `~/.bashrc` Guide

You can find the current default content of newly created user homes in `/etc/skel.bih`:

```bash
hpc-login-1:~$ head /etc/skel.bih/.bash*
==> /etc/skel.bih/.bash_logout <==
# ~/.bash_logout

==> /etc/skel.bih/.bash_profile <==
# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin

==> /etc/skel.bih/.bashrc <==
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=
```
