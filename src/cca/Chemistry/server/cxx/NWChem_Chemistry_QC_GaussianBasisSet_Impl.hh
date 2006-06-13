// 
// File:          NWChem_Chemistry_QC_GaussianBasisSet_Impl.hh
// Symbol:        NWChem.Chemistry_QC_GaussianBasisSet-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Chemistry_QC_GaussianBasisSet
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Chemistry_QC_GaussianBasisSet-v0.4.xml
// 

#ifndef included_NWChem_Chemistry_QC_GaussianBasisSet_Impl_hh
#define included_NWChem_Chemistry_QC_GaussianBasisSet_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_GaussianBasisSet_IOR_h
#include "NWChem_Chemistry_QC_GaussianBasisSet_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_Molecule_hh
#include "Chemistry_Molecule.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hh
#include "Chemistry_QC_GaussianBasis_AngularType.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Atomic_hh
#include "Chemistry_QC_GaussianBasis_Atomic.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_GaussianBasisSet_hh
#include "NWChem_Chemistry_QC_GaussianBasisSet.hh"
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


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianBasisSet._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianBasisSet._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianBasisSet._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.Chemistry_QC_GaussianBasisSet" (version 0.4)
   */
  class Chemistry_QC_GaussianBasisSet_impl
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianBasisSet._inherits)
  // Insert-Code-Here {NWChem.Chemistry_QC_GaussianBasisSet._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianBasisSet._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_GaussianBasisSet self;

    // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianBasisSet._implementation)
    // Insert-Code-Here {NWChem.Chemistry_QC_GaussianBasisSet._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianBasisSet._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_GaussianBasisSet_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_GaussianBasisSet_impl( struct 
      NWChem_Chemistry_QC_GaussianBasisSet__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_GaussianBasisSet_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Get the user specified name.
     * @return name 
     */
    ::std::string
    get_label() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the number of basis functions.
     * @return number of functions 
     */
    int64_t
    get_n_basis() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the number of shells.
     * @return number of shells 
     */
    int64_t
    get_n_shell() throw ( 
      ::sidl::BaseException
    );

    /**
     * Get the max angular momentum for any contraction in the 
     * basis set.
     * @return max angular momentum 
     */
    int32_t
    get_max_angular_momentum() throw ( 
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
     * Get an atomic basis set.
     * @param atomnum atom number 
     * @return Atomic 
     */
    ::Chemistry::QC::GaussianBasis::Atomic
    get_atomic (
      /* in */ int64_t atomnum
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get the molecule.
     * @return Molecule 
     */
    ::Chemistry::Molecule
    get_molecule() throw ( 
      ::sidl::BaseException
    );

    /**
     * Print the molecular basis data. 
     */
    void
    print_molecular() throw ( 
      ::sidl::BaseException
    );
  };  // end class Chemistry_QC_GaussianBasisSet_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_GaussianBasisSet._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_GaussianBasisSet._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_GaussianBasisSet._misc)

#endif
