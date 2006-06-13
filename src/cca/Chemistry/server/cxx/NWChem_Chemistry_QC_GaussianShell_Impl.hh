// 
// File:          NWChem_Chemistry_QC_GaussianShell_Impl.hh
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
#ifndef included_sidl_BaseException_hh
#include "sidl_BaseException.hh"
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
     * @return number of contractions 
     */
    int32_t
    get_n_contraction() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the number of primitives in the shell.
     * @return number of primitives 
     */
    int32_t
    get_n_primitive() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the coefficient for an unnormalized primitive 
     * in a contraction.
     * @param connum contraction number
     * @param expnum primitive number
     * @return contraction coefficient 
     */
    double
    get_contraction_coef (
      /* in */ int32_t connum,
      /* in */ int32_t expnum
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get the exponent for a primitive.
     * @param expnum primitive id number
     * @return exponent 
     */
    double
    get_exponent (
      /* in */ int32_t expnum
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get the angular momentum for a single contraction.
     * @param connum contraction id number
     * @return angular momentum value 
     */
    int32_t
    get_angular_momentum (
      /* in */ int32_t connum
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get the max angular momentum, considering all contractions 
     * in the shell.
     * @return maximum angular momentum value 
     */
    int32_t
    get_max_angular_momentum() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the angular type for a single contraction.
     * @param connum contraction number
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_contraction_angular_type (
      /* in */ int32_t connum
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get the angular type.
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type() throw ( 
      ::sidl::BaseException
    );

    /**
     * Print the shell data. 
     */
    void
    print_shell() throw ( 
      ::sidl::BaseException
    );
  };  // end class Chemistry_QC_GaussianShell_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianShell._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianShell._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianShell._misc)

#endif
