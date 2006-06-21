// 
// File:          NWChem_Chemistry_QC_ModelFactory_Impl.hh
// Symbol:        NWChem.Chemistry_QC_ModelFactory-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Chemistry_QC_ModelFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/vidhya/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Chemistry_QC_ModelFactory-v0.4.xml
// 

#ifndef included_NWChem_Chemistry_QC_ModelFactory_Impl_hh
#define included_NWChem_Chemistry_QC_ModelFactory_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_ModelFactory_IOR_h
#include "NWChem_Chemistry_QC_ModelFactory_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_Molecule_hh
#include "Chemistry_Molecule.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory_hh
#include "Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh"
#endif
#ifndef included_Chemistry_QC_Model_hh
#include "Chemistry_QC_Model.hh"
#endif
#ifndef included_NWChem_Chemistry_QC_ModelFactory_hh
#include "NWChem_Chemistry_QC_ModelFactory.hh"
#endif
#ifndef included_gov_cca_CCAException_hh
#include "gov_cca_CCAException.hh"
#endif
#ifndef included_gov_cca_Services_hh
#include "gov_cca_Services.hh"
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


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._includes} (includes or arbitrary code)
#include <string>
#include "cca.h"
#include <gov_cca_ports_ParameterPortFactory.hh>
#include <gov_cca_ports_ParameterPort.hh>
#include "dc/babel/babel-cca/server/ccaffeine_TypeMap.hh"
#include "util/IO.h"
//#include "jc++/jc++.h"
//#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "port/portInterfaces.h"
#include "port/supportInterfaces.h"
#include "NWChem_Chemistry_QC_Model_Impl.hh"

#include "Chemistry_MoleculeFactory.hh"
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._includes)

namespace NWChem { 

  /**
   * Symbol "NWChem.Chemistry_QC_ModelFactory" (version 0.4)
   */
  class Chemistry_QC_ModelFactory_impl
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._inherits)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_ModelFactory self;

    // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._implementation)
    // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._implementation} (additional details)
    std::string theory_;
    std::string basis_;
    std::string molecule_filename_;
    std::string nwchem_filename_;
    std::string config_filename_;
    std::string scratch_directory_;
    Chemistry::Molecule molecule_;
    Chemistry::MoleculeFactory molecule_factory_;

    gov::cca::Services services_;

    //
    // parameter stuff
    //
    gov::cca::Services myServices;
    gov::cca::ports::ParameterPortFactory ppf_;
    gov::cca::ports::ParameterPort pp_;
    
    StringParameter *scratchParameter;
    StringParameter *coordParameter;
    StringParameter *configParameter;
    StringParameter *basisSetParameter;
    StringParameter *theoryParameter;
    BoolParameter *utest;
    bool dynTestDone;  // dynamic parameter test

        //ConfigurableParameterPort*
        //setupParameters(ConfigurableParameterFactory *cfp);
        //extra public method
  public:

        //bool updateParameterPort(ConfigurableParameterPort *opp);

    // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_ModelFactory_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_ModelFactory_impl( struct 
      NWChem_Chemistry_QC_ModelFactory__object * s ) : self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_ModelFactory_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Set the theory name for Model's created with get_model.
     * @param theory A string giving the name of the theory, 
     * for example, B3LYP.
     */
    void
    set_theory (
      /* in */ const ::std::string& theory
    )
    throw () 
    ;


    /**
     * Set the basis set name for Model's created with get_model.
     * @param basis The basis set name to use, for example, aug-cc-pVDZ.
     */
    void
    set_basis (
      /* in */ const ::std::string& basis
    )
    throw () 
    ;


    /**
     * Set the Molecule to use for Model's created with get_model.
     * @param molecule An object of type Molecule.
     */
    void
    set_molecule (
      /* in */ ::Chemistry::Molecule molecule
    )
    throw () 
    ;


    /**
     * Set the object to use to compute integrals for Model's 
     * created with get_model.
     * @param intfact An object of type 
     * GaussianBasis.IntegralEvaluatorFactory.
     */
    void
    set_integral_factory (
      /* in */ ::Chemistry::QC::GaussianBasis::IntegralEvaluatorFactory intfact
    )
    throw () 
    ;


    /**
     * Returns a newly created Model.  Before get_model can be called, 
     * set_theory, set_basis, and set_molecule must be called.
     * @return The new Model instance.
     */
    ::Chemistry::QC::Model
    get_model() throw ( 
      ::sidl::BaseException
    );

    /**
     * This can be called when this Model object is no longer needed.  
     * No other members may be called after finalize. 
     */
    int32_t
    finalize() throw () 
    ;

    /**
     * Starts up a component presence in the calling framework.
     * @param services the component instance's handle on the framework world.
     * Contracts concerning Svc and setServices:
     * 
     * The component interaction with the CCA framework
     * and Ports begins on the call to setServices by the framework.
     * 
     * This function is called exactly once for each instance created
     * by the framework.
     * 
     * The argument Svc will never be nil/null.
     * 
     * Those uses ports which are automatically connected by the framework
     * (so-called service-ports) may be obtained via getPort during
     * setServices.
     */
    void
    setServices (
      /* in */ ::gov::cca::Services services
    )
    throw ( 
      ::gov::cca::CCAException
    );

  };  // end class Chemistry_QC_ModelFactory_impl

} // end namespace NWChem

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._misc)

#endif
