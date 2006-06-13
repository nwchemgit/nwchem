// 
// File:          NWChem_Chemistry_QC_Model_Impl.hh
// Symbol:        NWChem.Chemistry_QC_Model-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Chemistry_QC_Model
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Chemistry_QC_Model-v0.4.xml
// 

#ifndef included_NWChem_Chemistry_QC_Model_Impl_hh
#define included_NWChem_Chemistry_QC_Model_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_Model_IOR_h
#include "NWChem_Chemistry_QC_Model_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_Molecule_hh
#include "Chemistry_Molecule.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_Model_hh
#include "NWChem_Chemistry_QC_Model.hh"
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


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_Model._includes} (includes or arbitrary code)
#include <string>
#include "Chemistry_Chemistry_Molecule.hh"
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.Chemistry_QC_Model" (version 0.4)
   */
  class Chemistry_QC_Model_impl
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._inherits)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_Model self;

    // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._implementation)
    // Insert-Code-Here {NWChem.Chemistry_QC_Model._implementation} (additional details)
    Chemistry::Chemistry_Molecule molecule_;
    // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_Model_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_Model_impl( struct NWChem_Chemistry_QC_Model__object * s ) : 
      self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_Model_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    initialize (
      /* in */ const ::std::string& scratch_directory
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    change_theory (
      /* in */ const ::std::string& theory
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    change_basis (
      /* in */ const ::std::string& basis
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    setCoordinatesFromFile (
      /* in */ const ::std::string& molecule_filename
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    int32_t
    getNumCoordinates() throw () 
    ;
    /**
     * user defined non-static method.
     */
    ::sidl::array<double>
    get_coor() throw () 
    ;
    /**
     * user defined non-static method.
     */
    void
    set_coor (
      /* in */ ::sidl::array<double> x
    )
    throw () 
    ;


    /**
     * Set the molecule. @param molecule The new molecule. 
     */
    void
    set_molecule (
      /* in */ ::Chemistry::Molecule molecule
    )
    throw () 
    ;


    /**
     * Returns the molecule.  @return The Molecule object. 
     */
    ::Chemistry::Molecule
    get_molecule() throw () 
    ;
    /**
     * user defined non-static method.
     */
    double
    get_energy() throw ( 
      ::sidl::BaseException
    );

    /**
     * Sets the accuracy for subsequent energy calculations.
     * @param acc The new accuracy. 
     */
    void
    set_energy_accuracy (
      /* in */ double acc
    )
    throw () 
    ;


    /**
     * Returns the accuracy to which the energy is already computed.
     * The result is undefined if the energy has not already 
     * been computed.
     * @return The energy accuracy. 
     */
    double
    get_energy_accuracy() throw () 
    ;

    /**
     * This allows a programmer to request that if any result 
     * is computed,
     * then the energy is computed too.  This allows, say, for a request
     * for a gradient to cause the energy to be computed.  This computed
     * energy is cached and returned when the get_energy() member 
     * is called.
     * @param doit Whether or not to compute the energy.
     */
    void
    set_do_energy (
      /* in */ bool doit
    )
    throw () 
    ;


    /**
     * Returns the Cartesian gradient.  
     */
    ::sidl::array<double>
    get_gradient() throw ( 
      ::sidl::BaseException
    );

    /**
     * Sets the accuracy for subsequent gradient calculations
     * @param acc The new accuracy for gradients. 
     */
    void
    set_gradient_accuracy (
      /* in */ double acc
    )
    throw () 
    ;


    /**
     * Returns the accuracy to which the gradient is already computed.
     * The result is undefined if the gradient has not already 
     * been computed.
     * @return The current gradient accuracy. 
     */
    double
    get_gradient_accuracy() throw () 
    ;

    /**
     * Returns the Cartesian Hessian. @return The Hessian. 
     */
    ::sidl::array<double>
    get_hessian() throw ( 
      ::sidl::BaseException
    );

    /**
     * Sets the accuracy for subsequent Hessian calculations.
     * @param acc The new accuracy for Hessians. 
     */
    void
    set_hessian_accuracy (
      /* in */ double acc
    )
    throw () 
    ;


    /**
     * Returns the accuracy to which the Hessian is already computed.
     * The result is undefined if the Hessian has not already 
     * been computed. 
     */
    double
    get_hessian_accuracy() throw () 
    ;

    /**
     * Returns a Cartesian guess Hessian. 
     */
    ::sidl::array<double>
    get_guess_hessian() throw ( 
      ::sidl::BaseException
    );

    /**
     * Sets the accuracy for subsequent guess Hessian calculations.
     * @param acc The new accuracy for guess Hessians. 
     */
    void
    set_guess_hessian_accuracy (
      /* in */ double acc
    )
    throw () 
    ;


    /**
     * Returns the accuracy to which the guess Hessian is 
     * already computed.  The result is undefined if the guess Hessian 
     * has not already been computed.
     * @return The guess hessian accuracy.  
     */
    double
    get_guess_hessian_accuracy() throw () 
    ;

    /**
     * This can be called when this Model object is no longer needed.  
     * No other members may be called after finalize. 
     */
    int32_t
    finalize() throw () 
    ;
  };  // end class Chemistry_QC_Model_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_Model._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._misc)

#endif
