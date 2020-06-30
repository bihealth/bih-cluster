# Background In: Systems

This section gives a **very brief** overview on how modern computers work from the perspective of operating systems.
Hardware is the topic of the previous section [Background In: Computers](bg-computers.md).
Consider it an "abridge glossary" as the topic can fill lecture series and whole degrees.
The motivation for writing this is giving you "just enough knowledge to have a chance of *really* grasping what is going on a HPC cluster*.
The following two books are considered as standard literature; the authors neither consider them the best or definitive books, you'll probably be able to dig up more terse and better material through Wikipedia or your favourite search engine.
You might even find PDFs of the books or excerpts online.

- [John L. Hennessy & David A. Patterson -- Computer Architecture: A Quantitative Approach](https://www.amazon.com/dp/012383872X) (Hardware)
- [Andrew S. Tanenbaum -- Modern Operating Systems](https://www.amazon.com/dp/013359162X) (Software)
- [Michael Sipser -- Introduction to the Theory of Computation](https://www.amazon.com/dp/113318779X) (Models)

Consider this a *tour de force*, feel free to propose your fixes/additions through Github pull requests.

## The Three Big Abstractions

Computer history is **fast history** and only a few decades back, the way that computers were programmed differed **greatly** (the authors recommend getting an Arduino computer for a fun away to find out how it used to be).
We will not go into detail here, but these days we have the following "Big Three" of computing abstractions:

- the File.
- the Process.
- the Address Space.

What's the big deal?

### The File

Instead of writing data byte by byte (or maybe worse: block by block), you can just give your operating system a path to the file and just chose to read or write from it.
The worst that can happen is that you want to read from a non-existing file or write to a file in a non-existing directory.
Imagine a world in which there is no "name" but only numbers.
View a video that describes how "Legend of Zelda" for the "Nintendo Entertainment System" was programmed and you will know how lucky we are to live in the future.

### The Process

Your. Web. Browser. Has. Multiple. Tabs.

Each runs in its own "process" and is independent of the others.
Old-times computers were only able to exeucute one program at the same time.
Personal Computers (PCs) in the 1980s and 1990s were able to execute multiple programs by the programs cooperating and explicitely giving up control to the next.
These days, from your phone to the BIH HPC, the operating system ensures that no program just hogs the CPU and memory and the other programs cannot proceed.

### The Address Space

Similar to "The Process", the "Address Space" ensures that two programs cannot interfere.
Instead of isolating the instructions to be executed, the address space ensures that two programs cannot read each other's memory (your Youtube tab cannot read your online backing tab's memory).

### Bonus: Software Threads

Of course, one program can execute two "software threads" which can cooperate by sharing memory.
This way, your Youtube tab can play video and audio at the same time and use two of your hardware cores/threads to play the video.
Or maybe even the video decoder can use multiple threads at the same time.
All threads can communicate using shared memory of the current nodes.

Or: your NGS read mapping program can use multiple thread processing the millions of sequence.
Or: your image processing software can use multiple threads working on different parts of your image etc. etc.

### Bonus: Message Passing

How would two program instances running on two different computers (that do not share memory) work together?
The answer is: network communication / message passing!

The program instances need to know about each other (practically, this is taken care of by the cluster scheduler that sets up the environment).
Then, they they will send each other explicit messages over the network.
Clever network protocols use something called **remote direct memory access (RDMA)** which allows to transfer messages and (possibly large amounts of) data to the recipient with hardware support, allowing for low latencies and fast transfer rates.
