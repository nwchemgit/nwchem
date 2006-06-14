// 
// File:          NWChem_IntegralEvaluator1_Impl.cc
// Symbol:        NWChem.IntegralEvaluator1-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.IntegralEvaluator1
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.IntegralEvaluator1-v0.4.xml
// 
#include "NWChem_IntegralEvaluator1_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1._includes)
// Insert-Code-Here {NWChem.IntegralEvaluator1._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1._includes)

// user-defined constructor.
void NWChem::IntegralEvaluator1_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1._ctor)
  // Insert-Code-Here {NWChem.IntegralEvaluator1._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1._ctor)
}

// user-defined destructor.
void NWChem::IntegralEvaluator1_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1._dtor)
  // Insert-Code-Here {NWChem.IntegralEvaluator1._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1._dtor)
}

// static class initializer.
void NWChem::IntegralEvaluator1_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1._load)
  // Insert-Code-Here {NWChem.IntegralEvaluator1._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
NWChem::IntegralEvaluator1_impl::get_buffer (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.get_buffer)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.get_buffer} (get_buffer method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.get_buffer)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
NWChem::IntegralEvaluator1_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.get_descriptor)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.get_descriptor} (get_descriptor method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.get_descriptor)
}

/**
 * Compute a shell singlet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param deriv_level Derivative level. 
 */
void
NWChem::IntegralEvaluator1_impl::compute (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.compute)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.compute} (compute method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.compute)
}

/**
 * Compute a shell singlet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
NWChem::IntegralEvaluator1_impl::compute_array (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.compute_array)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.compute_array} (compute_array method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.compute_array)
}

/**
 * Compute integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
double
NWChem::IntegralEvaluator1_impl::compute_bounds (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.compute_bounds)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.compute_bounds} (compute_bounds method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.compute_bounds)
}

/**
 * Compute array of integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
::sidl::array<double>
NWChem::IntegralEvaluator1_impl::compute_bounds_array (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1.compute_bounds_array)
  // Insert-Code-Here {NWChem.IntegralEvaluator1.compute_bounds_array} (compute_bounds_array method)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1.compute_bounds_array)
}


// DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator1._misc)
// Insert-Code-Here {NWChem.IntegralEvaluator1._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator1._misc)

