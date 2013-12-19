/* This file contains Doxygen documentation only, and no code at all.
   It provides general information about the task module.
*/

/**
\defgroup task Task interface
\brief Interface to control the main tasks

In NWChem the overall work of a job is broken up into tasks. 
In general, each task is defined by an energy expression and a quantity to
calculate with it. For example, `task scf optimize` calculates a minimum energy
molecular structure using the Hartree-Fock energy expression.

The routines in this module implement how tasks are actually accomplished.
I.e. the task in the above example is implemented in the function 
`task_optimize`.
*/
/* $Id$ */
