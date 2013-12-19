/* This file only contains Doxygen documentation about the RunTime DataBase.
   This file does not contain any source code at all.
*/

/**
\defgroup rtdb The RunTime DataBase (RTDB)
\brief Notes on the RunTime DataBase

The run time data base is the parameter and information repository for
the independent modules (e.g., SCF, RIMP2) comprising NWChem.
This approach is similar in spirit to the GAMESS
dumpfile or the Gaussian checkpoint file.  The only way modules can
share data is via the database or via files, the names of which are stored in
the database (and may have default values).  Information is stored
directly in the database as typed arrays, each of which is described by

- a name, which is a simple string of ASCII characters (e.g., 
  `reference energies`),

- the type of the data (real, integer, logical, or character), 

- the number of data items, and

- the actual data (an array of items of the specified type).

A database is simply a file and is opened by name. Usually there is
just one database per calculation, though multiple databases may be
open at any instant.  

By default, access to all open databases occur in parallel, meaning
that

- all processes must participate in any read/write of any database
  and any such operation has an implied synchronization

- writes to the database write the data associated with process
  zero but the correct status of the operation is returned to all
  processes

- reads from the database read the data named by process zero and
  broadcast the data to all processes, checking dimensions and types
  of provided arrays

Alternatively, database operations can occur sequentially.  This means
that only process zero can read/write the database, and this happens
with no communication or synchronization with other processes.
Any read/write operations by any process other than process zero is
an error.

Usually, all processes will want the same data at the same time from
the database, and all processes will want to know of the success or
failure of operations.  This is readily done in the default parallel
mode.  An exception to this is during the reading of input.

Usually, only process zero will read the input and needs to store the
data directly into the database without involving the other processes.
This is done using sequential mode.

The following subsections contain a detailed listing of the C and Fortran API.  
Programs using RTDB routines must include the appropriate header file; 
`rtdb.fh` for Fortran, or `rtdb.h` for C.   These files define the return
types for all RTDB functions.  In addition, `rtdb.fh` specifies the 
following parameters 

- `rtdb_max_key` --- an integer parameter that defines the maximum
  length of a character string key

- `rtdb_max_file` --- an integer parameter that defines the maximum
  length of a file name

The Fortran routines return logical values; .true.  on success, .false. 
on failure.  The C routines return integers; 1 on success, or 0 on failure.
All `rtdb_*` functions are also mirrored by routines `rtdb_par_*`
in which process 0 performs the operation and all other processes
are broadcast the result of a read and discard writes.

*/
/* $Id$ */
