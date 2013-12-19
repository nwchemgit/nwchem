/*
$Id$

This file contains only Doxygen documentation and no code at all.
The contents of this file is constructed from various bits of documentation 
found in this directory.
*/

/**
\defgroup bas Basis Set Objects
\brief Notes on Basis Set Objects

Based on notes by Rick A. Kendall and Robert J. Harrison (April/May 1994)

The basis set object/api was designed to get all information from the
basis set based on a unique handle.  The internal data structures
store only the unique tag information which is different from what a
symmetry unique atom might be.  

The implementation of the basis set objects is based on a number of
considerations:

- The basis set objects should be able to cope with both general and segmented
  contracted Gaussian basis sets. Other basis sets should be anticipated.

- Gaussian basis sets can be either entirely Cartesian or spherical harmonic
  but not a mixture of the two representations.

- The basis functions are associated with atomic tags, not with coordinates.
  The tags establish the link between the basis functions and the coordinates
  in a geometry.

- All shells associated with a tag will be numbered consecutively. Also all
  basis functions in a shell are numbered consecutively. This implies that
  all basis functions associated with a tag are number consecutively as well.

For reference Cartesian basis functions are defined as
\f{eqnarray*}{
   \chi_{klm}(\vec{R}_A;\vec{r}_1)
   &=&N(X_A-x_1)^k(Y_A-y_1)^l(Z_A-z_1)^m\sum_i c_i e^{-\alpha_i r_{A1}^2}
\f}
where \f$ \vec{R}_A \f$ refers to the position of atom \f$ A \f$, and 
\f$ \vec{r}_1 \f$ refers to the position of electron \f$ 1 \f$. The terms in the
summation are refered to as primitive Gaussians, the sum of them as a 
contraction. The coefficients \f$ c_i \f$ are the contraction coefficients and
the \f$ \alpha_i \f$ are referred to as the exponents. 
The polynomial factor in front of the contraction is the angular momentum part
of the basis function.
Finally, \f$ N \f$ is the normalization constant, usually defined by requiring 
that \f$ \langle \chi_{klm} | \chi_{klm} \rangle = 1 \f$.
The set of basis functions \f$ \chi_{klm} \f$ for which \f$ k+l+m \f$ is a 
given constant is called a shell. In some cases this concept has been
generalized a little. For example an SP-shell is a set of functions for which
\f$ k+l+m \f$ is either \f$ 0 \f$ or \f$ 1 \f$. The exponents are the same in
either case, but the contraction coefficients are different for the S functions
from the ones for the P functions. 

For spherical harmonic basis functions the contraction of the above definition
remains the same, just the angular momentum part is replaced. In practice
the angular momentum part of spherical harmonic basis functions is represented
by linear combinations of Cartesian basis functions. The linear combinations
are chosen such that only the basis functions with the maximum angular momentum
corresponding to \f$ k+l+m \f$ remain. The other functions of lower angular
momentum that appear in Cartesian shells are referred to as contaminants.

Internally the data associated with every basis set is stored in tables. This
was required to be able to store the data in the usual Fortran77 data types.
The advantage of this is that the data structure remains very flat which makes
it easy to write the data to an external data store and also to retrieve it
from there. 

However the way the data stored is slightly different on file than in memory. On
file only the definitions of the basis functions is stored (i.e. not mapped to
any atoms). In memory the mapping between basis functions and atoms is generated
to enable fast access to data.

*/
