// 
// File:          NWChem_Chemistry_QC_GaussianShell_Impl.cc
// Symbol:        NWChem.Chemistry_QC_GaussianShell-v0.4
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for NWChem.Chemistry_QC_GaussianShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "NWChem_Chemistry_QC_GaussianShell_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._includes)

// user-defined constructor.
void NWChem::Chemistry_QC_GaussianShell_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._ctor)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._ctor)
}

// user-defined destructor.
void NWChem::Chemistry_QC_GaussianShell_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._dtor)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._dtor)
}

// static class initializer.
void NWChem::Chemistry_QC_GaussianShell_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._load)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Get the number of contractions in the shell. 
 * @return Number of contractions. 
 */
int64_t
NWChem::Chemistry_QC_GaussianShell_impl::get_n_contraction ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_n_contraction)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_n_contraction} (get_n_contraction method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_n_contraction)
}

/**
 * Get the number of primitives in the shell.
 * @return Number of primitives. 
 */
int64_t
NWChem::Chemistry_QC_GaussianShell_impl::get_n_primitive ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_n_primitive)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_n_primitive} (get_n_primitive method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_n_primitive)
}

/**
 * Get the coefficient for an unnormalized primitive in a contraction.
 * @param connum Contraction number.
 * @param expnum Primitive number.
 * @return The contraction coefficient. 
 */
double
NWChem::Chemistry_QC_GaussianShell_impl::get_contraction_coef (
  /* in */ int64_t connum,
  /* in */ int64_t expnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_contraction_coef)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_contraction_coef} (get_contraction_coef method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_contraction_coef)
}

/**
 * Get the exponent for a primitive.
 * @param expnum The primitive number.
 * @return The exponent. 
 */
double
NWChem::Chemistry_QC_GaussianShell_impl::get_exponent (
  /* in */ int64_t expnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_exponent)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_exponent} (get_exponent method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_exponent)
}

/**
 * Get the angular momentum for a single contraction.
 * @param connum Contraction number.
 * @return Angular momentum value. 
 */
int64_t
NWChem::Chemistry_QC_GaussianShell_impl::get_angular_momentum (
  /* in */ int64_t connum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_angular_momentum)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_angular_momentum} (get_angular_momentum method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_angular_momentum)
}

/**
 * Get the max angular momentum of any contraction in the shell.
 * @return Maximum angular momentum value. 
 */
int64_t
NWChem::Chemistry_QC_GaussianShell_impl::get_max_angular_momentum ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum} (get_max_angular_momentum method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum)
}

/**
 * Get the angular type for a single contraction.
 * @param connum Contraction number.
 * @return enum AngularType {CARTESIAN,SPHERICAL,MIXED} 
 */
::Chemistry::QC::GaussianBasis::AngularType
NWChem::Chemistry_QC_GaussianShell_impl::get_contraction_angular_type (
  /* in */ int64_t connum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type} (get_contraction_angular_type method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type)
}

/**
 * Get the angular type for the shell.
 * @return enum AngularType {CARTESIAN,SPHERICAL,MIXED} 
 */
::Chemistry::QC::GaussianBasis::AngularType
NWChem::Chemistry_QC_GaussianShell_impl::get_angular_type ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_angular_type)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_angular_type} (get_angular_type method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_angular_type)
}

/**
 * Print the shell data. 
 */
void
NWChem::Chemistry_QC_GaussianShell_impl::print_shell ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.print_shell)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.print_shell} (print_shell method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.print_shell)
}


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._misc)

