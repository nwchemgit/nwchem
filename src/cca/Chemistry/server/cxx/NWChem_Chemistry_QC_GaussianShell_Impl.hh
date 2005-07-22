// 
// File:          NWChem_Chemistry_QC_GaussianShell_Impl.hh
// Symbol:        NWChem.Chemistry_QC_GaussianShell-v0.4
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for NWChem.Chemistry_QC_GaussianShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 

#ifndef included_NWChem_Chemistry_QC_GaussianShell_Impl_hh
#define included_NWChem_Chemistry_QC_GaussianShell_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_GaussianShell_IOR_h
#include "NWChem_Chemistry_QC_GaussianShell_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hh
#include "Chemistry_QC_GaussianBasis_AngularType.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_GaussianShell_hh
#include "NWChem_Chemistry_QC_GaussianShell.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.Chemistry_QC_GaussianShell" (version 0.4)
   */
  class Chemistry_QC_GaussianShell_impl
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._inherits)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_GaussianShell self;

    // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._implementation)
    // Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_GaussianShell_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_GaussianShell_impl( struct 
      NWChem_Chemistry_QC_GaussianShell__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_GaussianShell_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Get the number of contractions in the shell. 
     * @return Number of contractions. 
     */
    int64_t
    get_n_contraction() throw () 
    ;

    /**
     * Get the number of primitives in the shell.
     * @return Number of primitives. 
     */
    int64_t
    get_n_primitive() throw () 
    ;

    /**
     * Get the coefficient for an unnormalized primitive in a contraction.
     * @param connum Contraction number.
     * @param expnum Primitive number.
     * @return The contraction coefficient. 
     */
    double
    get_contraction_coef (
      /* in */ int64_t connum,
      /* in */ int64_t expnum
    )
    throw () 
    ;


    /**
     * Get the exponent for a primitive.
     * @param expnum The primitive number.
     * @return The exponent. 
     */
    double
    get_exponent (
      /* in */ int64_t expnum
    )
    throw () 
    ;


    /**
     * Get the angular momentum for a single contraction.
     * @param connum Contraction number.
     * @return Angular momentum value. 
     */
    int64_t
    get_angular_momentum (
      /* in */ int64_t connum
    )
    throw () 
    ;


    /**
     * Get the max angular momentum of any contraction in the shell.
     * @return Maximum angular momentum value. 
     */
    int64_t
    get_max_angular_momentum() throw () 
    ;

    /**
     * Get the angular type for a single contraction.
     * @param connum Contraction number.
     * @return enum AngularType {CARTESIAN,SPHERICAL,MIXED} 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_contraction_angular_type (
      /* in */ int64_t connum
    )
    throw () 
    ;


    /**
     * Get the angular type for the shell.
     * @return enum AngularType {CARTESIAN,SPHERICAL,MIXED} 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type() throw () 
    ;

    /**
     * Print the shell data. 
     */
    void
    print_shell() throw () 
    ;
  };  // end class Chemistry_QC_GaussianShell_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._misc)

#endif
