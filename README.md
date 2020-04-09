# BIH HPC Cluster Documentation

- You can find the built documentation here: https://bihealth.github.io/bih-cluster/

## Building the Documentation Locally

**Prerequisites**

```bash
host:~$ sudo pip install pipenv  # maybe pip3 install or python -m pip install pipenv
```

**Clone**

```bash
host:~$ git clone xxx
host:~$ cd bih-cluster
host:bih-cluster$ pipenv install
host:bih-cluster$ pipenv shell
(bih-cluster) host:bih-cluster$ cd bih-cluster
(bih-cluster) host:bih-cluster$ mkdocs serve
```
