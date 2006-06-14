// 
// File:          NWChem_IntegralEvaluator4_Impl.hh
// Symbol:        NWChem.IntegralEvaluator4-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.IntegralEvaluator4
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.IntegralEvaluator4-v0.4.xml
// 

#ifndef included_NWChem_IntegralEvaluator4_Impl_hh
#define included_NWChem_IntegralEvaluator4_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_IntegralEvaluator4_IOR_h
#include "NWChem_IntegralEvaluator4_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescr_hh
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralDescr_hh
#include "Chemistry_QC_GaussianBasis_IntegralDescr.hh"
#endif
#ifndef included_NWChem_IntegralEvaluator4_hh
#include "NWChem_IntegralEvaluator4.hh"
#endif
#ifndef included_sidl_BaseException_hh
#include "sidl_BaseException.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator4._includes)
// Insert-Code-Here {NWChem.IntegralEvaluator4._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator4._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.IntegralEvaluator4" (version 0.4)
   */
  class IntegralEvaluator4_impl
  // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator4._inherits)
  // Insert-Code-Here {NWChem.IntegralEvaluator4._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator4._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    IntegralEvaluator4 self;

    // DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator4._implementation)
    // Insert-Code-Here {NWChem.IntegralEvaluator4._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator4._implementation)

  private:
    // private default constructor (required)
    IntegralEvaluator4_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    IntegralEvaluator4_impl( struct NWChem_IntegralEvaluator4__object * s ) : 
      self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~IntegralEvaluator4_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Get buffer pointer for given type.
     * @return Buffer pointer. 
     */
    void*
    get_buffer (
      /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc
    )
    throw ( 
      ::sidl::BaseException
    );

    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
    get_descriptor() throw ( 
      ::sidl::BaseException
    );

    /**
     * Compute a shell quartet of integrals.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3.
     * @param shellnum4 Gaussian shell number 4.
     * @param deriv_level Derivative level. 
     */
    void
    compute (
      /* in */ int64_t shellnum1,
      /* in */ int64_t shellnum2,
      /* in */ int64_t shellnum3,
      /* in */ int64_t shellnum4
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Compute a shell quartet of integrals and return as a borrowed
     * sidl array.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3.
     * @param shellnum4 Gaussian shell number 4.
     * @return Borrowed sidl array buffer. 
     */
    ::sidl::array<double>
    compute_array (
      /* in */ int64_t shellnum1,
      /* in */ int64_t shellnum2,
      /* in */ int64_t shellnum3,
      /* in */ int64_t shellnum4
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Compute integral bounds.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3.
     * @param shellnum4 Gaussian shell number 4. 
     */
    double
    compute_bounds (
      /* in */ int64_t shellnum1,
      /* in */ int64_t shellnum2,
      /* in */ int64_t shellnum3,
      /* in */ int64_t shellnum4
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Compute array of integral bounds.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3.
     * @param shellnum4 Gaussian shell number 4. 
     */
    ::sidl::array<double>
    compute_bounds_array (
      /* in */ int64_t shellnum1,
      /* in */ int64_t shellnum2,
      /* in */ int64_t shellnum3,
      /* in */ int64_t shellnum4
    )
    throw ( 
      ::sidl::BaseException
    );

  };  // end class IntegralEvaluator4_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.IntegralEvaluator4._misc)
// Insert-Code-Here {NWChem.IntegralEvaluator4._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.IntegralEvaluator4._misc)

#endif
