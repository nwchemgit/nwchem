gprof2tree
==========

Author: Huub van Dam
Web: http://www.cse.scitech.ac.uk/ccg/software/gprof2tree/ [August 24, 2015]

Introduction
------------

In large scale code it takes a considerable effort to find your way around.
In addition tracing how the execution path change if particular directives
are added to the input is not easy. Yet in order to enhance a code knowing
what parts are involved in implementing a particular option is essential. 

This tool addresses this issue by taking the output from the profiler gprof
and prints the call tree of that part of the program that was executed. In
particular the subroutine invocations at every level of the call tree can be
sorted. This way a diff program can be used to compare the call trees of the
same program run with different input options, highlighting the particular
parts of the code that are involved in implementing a particular option.
Alternatively it provides a quick way to get an overview of where particular
routines are used in a code. 

Limitations
-----------

Important in examining call trees is the ability to prune the trees eliminating
all the shrubbery that you are not interested in. Currently you need to specify
every routine you want to suppress, there is no facility to use regular 
expressions. I.e. if you are not interested in the calls to the RTDB library
you currently cannot specify that you want to stop at "rtdb_*". Instead you
have to list every rtdb_ function separately. The only help in this case is that
you can list all the rtdb_ functions in one rtdb related file, and create 
similar files for each other module. You can then specify the list of files
related to all the modules that you are not interested in with the --stop-at
flag. However, this still is relatively painful compared to using regular 
expressions.

More information
----------------

See the Gprof2tree.pdf file.

