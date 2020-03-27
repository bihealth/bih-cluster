# BIH HPC Documentation

This is the documentation of the BIH high-performance compute clusters, namely **hpc4research** and **hpc4clinic**
This documentation is maintained bih BIH HPC IT, BIH CUBI (Core Unit Bioinformatics), and the user community.

!!! note ""
    - Acknowledge usage of the cluster in your manuscript as *"Computation has been performed on the HPC for Research/Clinic cluster of the Berlin Institute of Health"*.
    - Please add your publications using the cluster to [this list](miscellaneous/publication-list).

??? note "Phasellus posuere in sem ut cursus"
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod
    nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor
    massa, nec semper lorem quam in massa.

- info, todo
- tip
- success, check, done
- question, help, faq
- abstract, summary, tldr
- warning, caution, attention
- failure fail missing
- danger
- bug
- example, snippet
- quote cite

``` bash tab="Bash"
#!/bin/bash

echo "Hello world!"
```

``` c tab="C"
#include <stdio.h>

int main(void) {
  printf("Hello world!\n");
}
```

``` c++ tab="C++"
#include <iostream>

int main() {
  std::cout << "Hello world!" << std::endl;
  return 0;
}
```

``` c# tab="C#"
using System;

class Program {
  static void Main(string[] args) {
    Console.WriteLine("Hello world!");
  }
}
```


asdf

``` python hl_lines="3 4"
""" Bubble sort """
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

Lorem ipsum[^1] dolor sit amet, consectetur adipiscing elit.[^2]

[^1]: Lorem ipsum dolor sit amet, consectetur adipiscing elit.

[^2]:
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod
    nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor
    massa, nec semper lorem quam in massa.

## Getting Started

- [For the Impatient](getting-started/for-the-impatient)


## Maintenance

[See here](admin/maintenance-windows) for status/progress and the next planned maintenance windows. We will notify you one or two weeks in advance of any maintenance (except for emergencies).
