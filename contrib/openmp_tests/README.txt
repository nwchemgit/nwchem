OpenMP Tests
============

The purpose of this directory is collect a number of unit tests that explore
various issues that arise in programming a GA/MPI/OpenMP parallel program.

Issues to consider are:

1. Access to inter-node communication (either MPI or GA)
2. Thread safety of using:
   a. Common blocks
   b. Runtime database
   c. MA library
   d. I/O in general
3. Mechanisms for managing the work distribution
4. Load balancing
5. Thread management costs (startup time, wakeup time, etc.)

MPI/OpenMP
----------

MPI has different models for communication in relation to threads:
1. MPI_THREAD_MULTIPLE:   Fully multithreaded
2. MPI_THREAD_SERIALIZED: Only one thread at a time makes MPI calls
3. MPI_THREAD_FUNNELED:   Only the main thread makes MPI calls
4. MPI_THREAD_SINGLE:     Pure single threaded
Traditionally (i.e. with MPI-2) full multithreading support was rather limited
in practice in that the default MPI library typically did not support it.
Apparently with MPI-3 the matter should have improved a lot.

On the Knights Landing architecture one apparently does not get optimal 
performance with 1 MPI-rank per node. The main reason is that there are 
multiple communication channels available and every MPI rank can use only one.
So one needs multiple MPI-ranks per node.

The other aspect is how the work is organized within each process. Shan et al.
used a model in which every thread essentially operates as a single MPI-rank
apart from the fact that memory model behind it is a shared memory one. The
MPI thread model they used was MPI_THREAD_SERIALIZED. I don't like that idea
because of 2 reasons:
1. To keep load balancing overheads small this requires handing out small
   tasks which generates a lot of communication.
2. All communications need to occur in an OpenMP Single block which leads to
   a serialization of the code when a lot of threads are used.

Another way to do the calculation is to dedicate the master thread to doing
all the communication, memory allocations and other "tricky" things. It should
then break a macro-task up into micro-tasks and farm these out to the other
threads. I.e. use a master-slave model. This has as the advantage:
1. Large macro-tasks can be handed out reducing the total amount of
   communication while intra-node parallelization can keep the task elapse
   time low (avoiding increasing load balancing overheads). 
The main potential disadvantage is:
1. The master thread may become a bottleneck.

