// 
// File:          NWChem_Physics_Units_Impl.hh
// Symbol:        NWChem.Physics_Units-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Physics_Units
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Physics_Units-v0.4.xml
// 

#ifndef included_NWChem_Physics_Units_Impl_hh
#define included_NWChem_Physics_Units_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_Physics_Units_IOR_h
#include "NWChem_Physics_Units_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_NWChem_Physics_Units_hh
#include "NWChem_Physics_Units.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._includes)
// Insert-Code-Here {NWChem.Physics_Units._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(NWChem.Physics_Units._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.Physics_Units" (version 0.4)
   */
  class Physics_Units_impl
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._inherits)
  // Insert-Code-Here {NWChem.Physics_Units._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Physics_Units self;

    // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._implementation)
    // Insert-Code-Here {NWChem.Physics_Units._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(NWChem.Physics_Units._implementation)

  private:
    // private default constructor (required)
    Physics_Units_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Physics_Units_impl( struct NWChem_Physics_Units__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Physics_Units_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Initializes the units as a human readable string
     * options are "angstroms" or "bohr" 
     */
    void
    initialize (
      /* in */ const ::std::string& unitname
    )
    throw () 
    ;


    /**
     * Returns the units as a human readable string. 
     */
    ::std::string
    get_unit_name() throw () 
    ;

    /**
     * Converts from self's units to the given unit name. 
     */
    double
    convert_to (
      /* in */ const ::std::string& unitname
    )
    throw () 
    ;


    /**
     * Converts to self's units from the given unit name. 
     */
    double
    convert_from (
      /* in */ const ::std::string& unitname
    )
    throw () 
    ;

  };  // end class Physics_Units_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._misc)
// Insert-Code-Here {NWChem.Physics_Units._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Physics_Units._misc)

#endif
