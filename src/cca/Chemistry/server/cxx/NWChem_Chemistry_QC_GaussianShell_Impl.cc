// 
// File:          NWChem_Chemistry_QC_GaussianShell_Impl.cc
// Symbol:        NWChem.Chemistry_QC_GaussianShell-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Chemistry_QC_GaussianShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Chemistry_QC_GaussianShell-v0.4.xml
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
 * @return number of contractions 
 */
int32_t
NWChem::Chemistry_QC_GaussianShell_impl::get_n_contraction ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_n_contraction)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_n_contraction} (get_n_contraction method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_n_contraction)
}

/**
 * Get the number of primitives in the shell.
 * @return number of primitives 
 */
int32_t
NWChem::Chemistry_QC_GaussianShell_impl::get_n_primitive ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_n_primitive)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_n_primitive} (get_n_primitive method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_n_primitive)
}

/**
 * Get the coefficient for an unnormalized primitive 
 * in a contraction.
 * @param connum contraction number
 * @param expnum primitive number
 * @return contraction coefficient 
 */
double
NWChem::Chemistry_QC_GaussianShell_impl::get_contraction_coef (
  /* in */ int32_t connum,
  /* in */ int32_t expnum ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_contraction_coef)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_contraction_coef} (get_contraction_coef method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_contraction_coef)
}

/**
 * Get the exponent for a primitive.
 * @param expnum primitive id number
 * @return exponent 
 */
double
NWChem::Chemistry_QC_GaussianShell_impl::get_exponent (
  /* in */ int32_t expnum ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_exponent)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_exponent} (get_exponent method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_exponent)
}

/**
 * Get the angular momentum for a single contraction.
 * @param connum contraction id number
 * @return angular momentum value 
 */
int32_t
NWChem::Chemistry_QC_GaussianShell_impl::get_angular_momentum (
  /* in */ int32_t connum ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_angular_momentum)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_angular_momentum} (get_angular_momentum method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_angular_momentum)
}

/**
 * Get the max angular momentum, considering all contractions 
 * in the shell.
 * @return maximum angular momentum value 
 */
int32_t
NWChem::Chemistry_QC_GaussianShell_impl::get_max_angular_momentum ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum} (get_max_angular_momentum method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_max_angular_momentum)
}

/**
 * Get the angular type for a single contraction.
 * @param connum contraction number
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
NWChem::Chemistry_QC_GaussianShell_impl::get_contraction_angular_type (
  /* in */ int32_t connum ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type} (get_contraction_angular_type method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.get_contraction_angular_type)
}

/**
 * Get the angular type.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
NWChem::Chemistry_QC_GaussianShell_impl::get_angular_type ()
throw ( 
  ::sidl::BaseException
)
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
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell.print_shell)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell.print_shell} (print_shell method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell.print_shell)
}


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._misc)

