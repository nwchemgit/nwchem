/* This file contains only Doxygen documentation on geometry objects.
   There is no actual code here.
*/

/**
\defgroup geom The Geometry Object
\brief Notes on geometry objects

The geometry object is used in NWChem to store and manipulate important
information describing the molecular system to be modeled, not all of
which is specifically connected with the geometry of the system.  The
geometry object serves four main purposes; 

- provides a definition of the coordinate system and positioning in space
  (including lattice vectors for periodic systems)

- defines an association of names/tags with coordinates in space

- specifies the external potential (nuclear multipole
  moments, external fields, effective core potentials, ...) that
  define the Hamiltonian for all electronic structure methods

- stores most Hamiltonian related information (but
  not wavefunction related information).

The tag associated with a geometric center serves a number of purposes in 
NWChem. It provides a convenient and unambiguous way to refer to

- a specific chemical element (which provides default values for information
  such as nuclear charge, mass, number of electrons, ...)

- the name of an `atomic' basis set

- a DFT grid

The tag can also serve as a test for symmetry equivalence, since lower symmetry 
can be forced by specifying different tags for otherwise symmetry equivalent
centers.

The data contained in the geometry object (or information that can be derived
from data in the object) include the following;

- A description of the coordinates of all types of centers (e.g.,
  atom, basis function)

- Charges (or optionally, ECPs, ...) associated with those centers

- Tags (names) of centers

- Masses associated with centers

- Variables for optimization (e.g., via constrained cartesians
  or Z-matrix variables)

- Symmetry information

- Any other simple scalar/vector attribute associated
  specifically with a center

Specific geometries are referenced through an integer handle.   
Multiple geometries can be defined such that any one of them 
may be accessible at any instant for a given problem.  However,
geometries can consume a large amount of memory, so it is usually
advisable to keep the number of simultaneously `open' geometries to a minimum.

Logical functions return .true. on sucess, .false. on failure.  Below the
various functions are described in more detail.

*/
/* $Id$ */
