// 
// File:          NWChem_Chemistry_QC_ModelFactory_Impl.cc
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
#include "NWChem_Chemistry_QC_ModelFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._includes} (additional includes or code)
#include "sidl.hh"
#include "gov_cca_Port.hh"
#include "gov_cca_Services.hh"

// CCA stuff
#include "util/IO.h"
#include "dc/babel/babel-cca/server/ccaffeine_ports_PortTranslator.hh"
// Classic includes:
#include "cca.h"
#include <stdPorts.h>
#include "util/IO.h"
#include "parameters/parametersStar.h"
#include <unistd.h>
#include <fstream.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._includes)

// user-defined constructor.
void NWChem::Chemistry_QC_ModelFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._ctor)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._ctor} (constructor)
  #ifdef DEBUG
       cout << "NWChemFactory::Component created" << endl;
  #endif // DEBUG
  dynTestDone = false;
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._ctor)
}

// user-defined destructor.
void NWChem::Chemistry_QC_ModelFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._dtor)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._dtor} (destructor)
  #ifdef DEBUG
      cout << "NWChemFactory::Component destroyed" << endl;
  #endif // DEBUG
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._dtor)
}

// static class initializer.
void NWChem::Chemistry_QC_ModelFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._load)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Set the theory name for Model's created with get_model.
 * @param theory A string giving the name of the theory, 
 * for example, B3LYP.
 */
void
NWChem::Chemistry_QC_ModelFactory_impl::set_theory (
  /* in */ const ::std::string& theory ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.set_theory)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.set_theory} (set_theory method)
  theory_ = theory;
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.set_theory)
}

/**
 * Set the basis set name for Model's created with get_model.
 * @param basis The basis set name to use, for example, aug-cc-pVDZ.
 */
void
NWChem::Chemistry_QC_ModelFactory_impl::set_basis (
  /* in */ const ::std::string& basis ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.set_basis)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.set_basis} (set_basis method)
  basis_ = basis;
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.set_basis)
}

/**
 * Set the Molecule to use for Model's created with get_model.
 * @param molecule An object of type Molecule.
 */
void
NWChem::Chemistry_QC_ModelFactory_impl::set_molecule (
  /* in */ ::Chemistry::Molecule molecule ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.set_molecule)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.set_molecule} (set_molecule method)
  molecule_=molecule; 
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.set_molecule)
}

/**
 * Set the object to use to compute integrals for Model's 
 * created with get_model.
 * @param intfact An object of type 
 * GaussianBasis.IntegralEvaluatorFactory.
 */
void
NWChem::Chemistry_QC_ModelFactory_impl::set_integral_factory (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralEvaluatorFactory intfact ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.set_integral_factory)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.set_integral_factory} (set_integral_factory method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.set_integral_factory)
}

/**
 * Returns a newly created Model.  Before get_model can be called, 
 * set_theory, set_basis, and set_molecule must be called.
 * @return The new Model instance.
 */
::Chemistry::QC::Model
NWChem::Chemistry_QC_ModelFactory_impl::get_model ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.get_model)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.get_model} (get_model method)
  string symbols[]={
          "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne",
          "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca",
          "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
          "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr",
          "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
          "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
          "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
          "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg",
          "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
          "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
          "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"};
	  
   NWChem::Chemistry_QC_Model model = NWChem::Chemistry_QC_Model::_create();

   scratch_directory_="/scratch";
   scratch_directory_="/home/vidhya/tmp";
   theory_="scf";
   basis_="sto-3g";
   molecule_filename_="/home/vidhya/CCA/mcmd-paper/nwchem/src/cca/tests/water.cmp";

   // Getting values from the parameter ports
   gov::cca::TypeMap tm = pp_.readConfigurationMap();    
   theory_ = std::string( tm.getString("theory",
                                       "failed theory fetch") );
   basis_  = std::string( tm.getString("basis",
                                       "failed basis fetch") );
   molecule_filename_ =
      std::string( tm.getString("molecule_filename",
                                "failed molecule_filename fetch") ); 
   nwchem_filename_ =
      std::string( tm.getString("nwchem_filename",
                                "failed nwchem_filename fetch") ); 

// This is the old way of doing parameter ports     
#if 0
   if(theory_.size()==0)
      theory_ = std::string(theoryParameter->value);
   if(basis_.size()==0)
      basis_ = std::string(basisSetParameter->value);
   if(molecule_filename_.size()==0)
      molecule_filename_ = std::string(coordParameter->value);
   if(config_filename_.size()==0)
      config_filename_ = std::string(configParameter->value);
   if(scratch_directory_.size()==0)
     scratch_directory_ = std::string(scratchParameter->value);
#endif

   
#ifndef DEBUG
   cout << "Got from parameter port theory=" << theory_ << endl;
   cout << "Got from parameter port basis=" << basis_ << endl;
   cout << "Got from parameter port molecule_filename=" << molecule_filename_ << endl;
   cout << "Got from parameter port config_filename=" << config_filename_ << endl;
   cout << "Got from parameter port scratch_directory=" << scratch_directory_ << endl;
#endif // DEBUG

   // Set the actual values in the model
   model.initialize(scratch_directory_);
   model.change_theory(theory_);
   model.change_basis(basis_);

// Currently we need to read in the coordinates from an NWChem file so
// that the labels, masses, etc. get set up correctly
   model.setCoordinatesFromFile(nwchem_filename_); 

  /*
   Currently two options at this point:
     1) we are using python GUI, set_molecule() has already been called
        and !molecule evaluates to FALSE
     2) we are using caffeine proper, we execute the following block to get
        a molecule from the molecule factory
  */

// Setting up the molecule object in CCA
  if( !molecule_ ) {
  // pass coordinate file and get molecule from MoleculeFactory
      molecule_factory_ = services_.getPort("MoleculeFactory");
      molecule_factory_.set_molecule_filename(molecule_filename_);
      molecule_ = molecule_factory_.get_molecule();
  }

// Make sure that the model has the molecule object
  model.set_molecule(molecule_);
//  double conv = molecule_.get_units().convert_to("bohr");

  return model;
  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.get_model)
}

/**
 * This can be called when this Model object is no longer needed.  
 * No other members may be called after finalize. 
 */
int32_t
NWChem::Chemistry_QC_ModelFactory_impl::finalize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.finalize)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.finalize} (finalize method)
  if (molecule_factory_._not_nil())
      services_.releasePort("MoleculeFactory");  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.finalize)
}

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
NWChem::Chemistry_QC_ModelFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory.setServices)
  // Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory.setServices} (setServices method)
   molecule_=0;
   gov::cca::Port port = self;
   if (port._is_nil()) {
      cout << "unable to cast self as port NWChemFactory::....::setServices"
         << endl;
   }

   services_ = services;
   if (services_._is_nil())
   {
      cout << "services is nill in NWChemFactory::....::setServices"
         << endl;
      return;
   }

   try {
      services_.addProvidesPort(port, std::string("ModelFactory"),
            std::string("Chemistry.QC.ModelFactory"), 0);
   }
   catch (gov::cca::CCAException e) {
      std::cout << "Error using provide ModelFactory: " << e.getNote() << std::endl;
   }

   try {
      services_.registerUsesPort("MoleculeFactory",
                                 "Chemistry.MoleculeFactory", 0);
   }
   catch (gov::cca::CCAException e) {
      std::cout << "Error using use MoleculeFactory: " << e.getNote() << std::endl;
   }

   //
   //  Parameter stuff
   //
#if 1
   try {
      services_.registerUsesPort("ppf",
                                 "gov.cca.ports.ParameterPortFactory", 0);
      services_.registerUsesPort("BasisName",
                                 "Util.StringProvider", 0);
      services_.registerUsesPort("TheoryName",
                                 "Util.StringProvider", 0);
      services_.registerUsesPort("MoleculeFile",
                                 "Util.StringProvider", 0); 
   }
   catch (gov::cca::CCAException e) {
      std::cout << "Error using ParameterPort: " << e.getNote() << std::endl;
   }
   // setup parameters
   try {

      gov::cca::TypeMap tm_ = services_.createTypeMap(); // empty typemap
      if(tm_._is_nil()) {
         std::cerr << "TypeMap is nill\n";
         abort();
      }
      ppf_ = services_.getPort("ppf");
      ppf_.initParameterData(tm_, "CONFIG");
      ppf_.setBatchTitle(tm_,"NWChem ModelFactory Options");
      ppf_.setGroupName(tm_,"Job Specification");
      ppf_.addRequestString(tm_, "theory", "Theory name",
                            "Theory", "HF");
      ppf_.addRequestString(tm_, "basis", "AO basis name",
                            "Basis", "STO-3G");
      ppf_.addRequestString(tm_, "molecule_filename",
                            "Full path to molecule file",
                            "Molecule filename", "");
      ppf_.addRequestString(tm_, "nwchem_filename",
                            "Full path to nwchem file",
                            "NWChem filename", "");
      ppf_.addRequestString(tm_, "keyval_filename",
                            "Full path to keyval input file",
                            "Keyval filename", "");
      ppf_.addRequestString(tm_, "integral_buffer", "Integral buffer approach",
                            "Integral buffer", "opaque");

      ppf_.addParameterPort(tm_, services_);
      services_.releasePort("ppf");

      pp_ = services_.getPort("CONFIG");
      if (pp_._is_nil()) {
         std::cerr << "getport failed\n";
         abort();
      }

   }
   catch(std::exception& e) {
      std::cerr << "Error in parameter port setup: " << e.what() << std::endl;
   } 
#else
   try {
      this->services_.registerUsesPort("classicParam",
            "gov.cca.ParameterPortFactoryService", tm);
   } catch (gov::cca::CCAException e) {
      cout << "Error registering 'classicParam' usesPort in NWChemFactory: " <<
         e.getNote() << endl;
   }

   gov::cca::Port p = services_.getPort("classicParam");

   // Get the underlying classic parameter port
   ccaffeine::ports::PortTranslator portX = p;
   if (portX._not_nil()) {
      classic::gov::cca::Port *cp =
         static_cast<classic::gov::cca::Port *>(portX.getClassicPort());

      ConfigurableParameterFactory *cpf =
         dynamic_cast<ConfigurableParameterFactory *>(cp);

      pp = setupParameters(cpf);

      classic::gov::cca::Port *clscp = dynamic_cast<classic::gov::cca::Port*>(pp);
      if (clscp == NULL)
         cout << "Cannot cast ConfigurableParameterPort* to classic::gov::cca::Port*" << endl;

      // Pass the configurable parameter port to the Port Translator
      void * vp = static_cast<void*>(clscp);
      if (!vp) {
         cout << "Could not cast clscp to void*" << endl;
      }

      ccaffeine::ports::PortTranslator provideX =
         ccaffeine::ports::PortTranslator::createFromClassic(vp);

      // Add the configure port we want to provide
      services_.addProvidesPort(provideX, "configure", "ParameterPort", tm);

      // release classic parameter
      services_.releasePort("classicParam");
      services_.unregisterUsesPort("classicParam");
   }  
#endif
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_ModelFactory._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_ModelFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_ModelFactory._misc)

