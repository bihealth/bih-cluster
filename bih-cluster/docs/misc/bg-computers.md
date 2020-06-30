# Background In: Computers

This section gives a **very brief** overview on how modern computers work, the perspective of hardware.
Operating sytems are the topic of the next section [Background In: Systems](bg-systems.md).
Consider it an "abridge glossary" as the topic can fill lecture series and whole degrees.
The motivation for writing this is giving you "just enough knowledge to have a chance of *really* grasping what is going on a HPC cluster*.
The following two books are considered as standard literature; the authors neither consider them the best or definitive books, you'll probably be able to dig up more terse and better material through Wikipedia or your favourite search engine.
You might even find PDFs of the books or excerpts online.

- [John L. Hennessy & David A. Patterson -- Computer Architecture: A Quantitative Approach](https://www.amazon.com/dp/012383872X) (Hardware)
- [Andrew S. Tanenbaum -- Modern Operating Systems](https://www.amazon.com/dp/013359162X) (Software)
- [Michael Sipser -- Introduction to the Theory of Computation](https://www.amazon.com/dp/113318779X) (Models)

Consider this a *tour de force*, feel free to propose your fixes/additions through Github pull requests.

## The Random Access Machine (RAM)

The Random Access Machine ([Wikipedia](https://en.wikipedia.org/wiki/Random-access_machine)) is a useful model of computation.
Depending on your text book, there are many versions; the following is sufficient:

- A *program* in some kind of programming language (pointers (memory indirection), loops, and if/then/else are sufficient; variables, GOTO, and sub routines are useful but do not add power).
- An *instruction pointer* that points to the current instruction ("line").
- An infinite amount of linear *random access memory* that is addressable and has constant access time.

The RAM is a nice abstraction of modern computers.
The main difference is that real computers have *finite* memory and that access time is not linear.

We will not go into detail here but modern computers deviate in the following aspects (use the italic phrases as search term with "Wikipedia" to learn more).
Real computers have:

- *instruction pointer*,
- *processor registers*,
- *advanced vector extensions*,
- *multi-core processors*,
- *hardware threads*,
- *non-uniform memory access (NUMA)*,
- *translation lookaside buffers*,
- *instruction pieplining*, and much more.

Think of the RAM executing the program in your favourite programming language (as long as you do not use any threads or asynchronous features).
It does not matter whether it is Pyton, Java, C(++), or even Bash, the RAM is a nice abstraction of all of them.

## Central Processing Units (CPU)

Actual CPUs are "real-world" implementations of the RAM.
Ignoring parallelism (topic of the next section), they provide the following features:

- access to registers and main memory,
- CPU caches, organized in cache lines (a topic to its own),
- access to peripherials (e.g., hard drive disks, keyboards, the network, ...),
- executing a program (current position in the program is the *instruction pointer*),
- an actualy arithmetic logical unit (ALU) that allow for executing operations such as "plus", "minus", "if X == 0 then", etc.

The actual implementation might differ and important trade-offs will be taken from a CPU in a high-performance cluster vs. CPUs in your phone or a "smart" light bulp.
However, in their core they are all the same.

## Multi-Core and Multi-Thread CPUs

Until the 2000s, single-core computers would "double their speed" every two years because of something called *Moore's Law*.
The peaked with the Intel Pentium 4 which had an energy density that was somewhere close to that of an atomic fuel rod.
The days, the *number of transistors* (small element of memory and CPUs) per chip still doubles every two years but instead of CPUs becoming smaller and faster, we get multiple CPUs per chip.

Each of these is a **CPU core**.
For example, the BIH HPC contains CPUs that have 16 to 40 cores.
Each of these cores has its own arithmetic-logical unit (ALU) and can do its own independent computations.
(Actually this is only half of the truth, search the internet for *CPU/memory bus* and *von Neumann's bottleneck* to find out more).

The main characteristic of compute with multiple CPUs and CPU cores is that all CPUs and cores share the **same memory**.
That is, if they need to communicate they can just write/read data to the same location (following certain, "interesting", protocols) which is very fast.

Also since the 2000s, modern CPU architectures feature **(hardware) threads**.
Modern CPUs feature an overload of features and (mostly Intel) hardware architects realized that many of the hardware circuits were actually unused at any given time.
They just implemented a feature that would allow two programs to virtually run at the same time on the same CPU core (so-called "(hardware) threads").
For many workloads, they could double the amount of instructions that could be processed on a single core but the "Spectre" exploit (and many more until today) are the result of taking these shortcuts.

## Data Parallelism & General/Graphics Processing Units

The instruction pointer approach mentioned above implies that every program can only execute on instruction at a time.
However, in many applications one might want to do many of the same instruction at the same time.
Think about adding two vectors: we will want to add two series of millions of numbers that are all located next to each other in memory.
These are called "vector instructions", marketing called this "MMX", "SSE", and more recently "AVX".
The concept itself is called SIMD - single instruction, multiple data.

Instead of operating on *machine words* of 32 or 64 bits, one can make various operations on 512 bits at the same time.
Depending on the CPU architecture, one can operate on words of different length and also the type of operation differs (e.g., adding numbers was introduced relatively early while operations such as `popcnt` which counts the number of `1` bits were introduced relatively recently).

GPUs (general/graphics processing units) take the whole concept to the next level.
Originally designed for accelerating 3D graphics programs, they excel at performing the same operation on long vectors (*many* numbers in a row).
This all comes at a cost, of course, as their floating number resolution is limited.
Recently, the *graphics* has become *general* as the main vendor NVIDIA and third parties have come with clever libraries that accelerate non-classical applications such as deep learning (of course, there are also "classic" applications such as linear algebra which area also useful for machine learning).
