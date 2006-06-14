// 
// File:          NWChem_Chemistry_QC_IntEvalFactory_Impl.cc
// Symbol:        NWChem.Chemistry_QC_IntEvalFactory-v0.4
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for NWChem.Chemistry_QC_IntEvalFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/windus/CCA/mcmd-paper/nwchem/src/cca/repo/NWChem.Chemistry_QC_IntEvalFactory-v0.4.xml
// 
#include "NWChem_Chemistry_QC_IntEvalFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory._includes)

// user-defined constructor.
void NWChem::Chemistry_QC_IntEvalFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory._ctor)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory._ctor)
}

// user-defined destructor.
void NWChem::Chemistry_QC_IntEvalFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory._dtor)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory._dtor)
}

// static class initializer.
void NWChem::Chemistry_QC_IntEvalFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory._load)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  get_name[]
 */
::std::string
NWChem::Chemistry_QC_IntEvalFactory_impl::get_name ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_name)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_name} (get_name method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_name)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
NWChem::Chemistry_QC_IntEvalFactory_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_descriptor)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_descriptor} (get_descriptor method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_descriptor)
}

/**
 * Method:  is_supported[]
 */
bool
NWChem::Chemistry_QC_IntEvalFactory_impl::is_supported (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.is_supported)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.is_supported} (is_supported method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.is_supported)
}

/**
 * Set available storage
 * @param storage Available storage in bytes 
 */
void
NWChem::Chemistry_QC_IntEvalFactory_impl::set_storage (
  /* in */ int64_t storage ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.set_storage)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.set_storage} (set_storage method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.set_storage)
}

/**
 * Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1
NWChem::Chemistry_QC_IntEvalFactory_impl::get_evaluator1 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator1)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_evaluator1} (get_evaluator1 method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator1)
}

/**
 * Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
NWChem::Chemistry_QC_IntEvalFactory_impl::get_evaluator2 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator2)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_evaluator2} (get_evaluator2 method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator2)
}

/**
 * Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3
NWChem::Chemistry_QC_IntEvalFactory_impl::get_evaluator3 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator3)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_evaluator3} (get_evaluator3 method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator3)
}

/**
 * Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4
NWChem::Chemistry_QC_IntEvalFactory_impl::get_evaluator4 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator4)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.get_evaluator4} (get_evaluator4 method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.get_evaluator4)
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
NWChem::Chemistry_QC_IntEvalFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory.setServices)
  // Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory.setServices} (setServices method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_IntEvalFactory._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_IntEvalFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_IntEvalFactory._misc)

