// 
// File:          NWChem_Chemistry_QC_Model_Impl.cc
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
#include "NWChem_Chemistry_QC_Model_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._includes)
// Insert-Code-Here {NWChem.Chemistry_QC_Model._includes} (additional includes or code)
#include "NWChemWrap.fh"
#include <string>
#include <cctype>

#include <iostream>
using namespace std;
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._includes)

// user-defined constructor.
void NWChem::Chemistry_QC_Model_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._ctor)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._ctor)
}

// user-defined destructor.
void NWChem::Chemistry_QC_Model_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._dtor)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._dtor)
}

// static class initializer.
void NWChem::Chemistry_QC_Model_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._load)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize[]
 */
void
NWChem::Chemistry_QC_Model_impl::initialize (
  /* in */ const ::std::string& scratch_directory ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.initialize)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.initialize} (initialize method)
  int len;
  len=scratch_directory.length();
  nwchem_nwchemstart_(scratch_directory.c_str(),len);
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.initialize)
}

/**
 * Method:  change_theory[]
 */
void
NWChem::Chemistry_QC_Model_impl::change_theory (
  /* in */ const ::std::string& theory ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.change_theory)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.change_theory} (change_theory method)
  nwchem_settheory_(theory.c_str());  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.change_theory)
}

/**
 * Method:  change_basis[]
 */
void
NWChem::Chemistry_QC_Model_impl::change_basis (
  /* in */ const ::std::string& basis ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.change_basis)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.change_basis} (change_basis method)
  int len;
  len=basis.length();
  nwchem_setbasisset_(basis.c_str(),len);  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.change_basis)
}

/**
 * Method:  setCoordinatesFromFile[]
 */
void
NWChem::Chemistry_QC_Model_impl::setCoordinatesFromFile (
  /* in */ const ::std::string& molecule_filename ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.setCoordinatesFromFile)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.setCoordinatesFromFile} (setCoordinatesFromFile method)
  std::cout << "\n\nNWCHEM COORDS FROM FILE: " << molecule_filename;
  nwchem_setcoordinatesfromfile_(molecule_filename.c_str()); 
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.setCoordinatesFromFile)
}

/**
 * Method:  getNumCoordinates[]
 */
int32_t
NWChem::Chemistry_QC_Model_impl::getNumCoordinates ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.getNumCoordinates)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.getNumCoordinates} (getNumCoordinates method)
  int num;
  nwchem_getnumcoordinates_(&num);
  return num;
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.getNumCoordinates)
}

/**
 * Method:  get_coor[]
 */
::sidl::array<double>
NWChem::Chemistry_QC_Model_impl::get_coor ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_coor)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_coor} (get_coor method)
  int num;
  nwchem_getnumcoordinates_(&num);

  sidl::array<double> x = sidl::array<double>::create1d(num);
  double* dataPtr = x.first();
  nwchem_getcoordinates_(dataPtr);
  return x;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_coor)
}

/**
 * Method:  set_coor[]
 */
void
NWChem::Chemistry_QC_Model_impl::set_coor (
  /* in */ ::sidl::array<double> x ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_coor)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_coor} (set_coor method)
  double* dataPtr = x.first();
  nwchem_setcoordinates_(dataPtr); 
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_coor)
}

/**
 * Set the molecule. @param molecule The new molecule. 
 */
void
NWChem::Chemistry_QC_Model_impl::set_molecule (
  /* in */ ::Chemistry::Molecule molecule ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_molecule)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_molecule} (set_molecule method)
  molecule_ = molecule;
  ::sidl::array<double> x;
  int coord_num;

  double conv = molecule_.get_units().convert_to("bohr");
  coord_num=molecule_.get_n_atom()*3;
  x=molecule_.get_coor();
  double* dataPtr = x.first();

  for(int i=0;i<coord_num;++i)
  {
    *(dataPtr+i)=(*(dataPtr+i))*conv;
  }

  nwchem_setcoordinates_(dataPtr);

  return;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_molecule)
}

/**
 * Returns the molecule.  @return The Molecule object. 
 */
::Chemistry::Molecule
NWChem::Chemistry_QC_Model_impl::get_molecule ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_molecule)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_molecule} (get_molecule method)
  return molecule_;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_molecule)
}

/**
 * Method:  get_energy[]
 */
double
NWChem::Chemistry_QC_Model_impl::get_energy ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_energy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_energy} (get_energy method)
  double f;
  nwchem_taskenergy_(&f);
  return f;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_energy)
}

/**
 * Sets the accuracy for subsequent energy calculations.
 * @param acc The new accuracy. 
 */
void
NWChem::Chemistry_QC_Model_impl::set_energy_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_energy_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_energy_accuracy} (set_energy_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_energy_accuracy)
}

/**
 * Returns the accuracy to which the energy is already computed.
 * The result is undefined if the energy has not already 
 * been computed.
 * @return The energy accuracy. 
 */
double
NWChem::Chemistry_QC_Model_impl::get_energy_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_energy_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_energy_accuracy} (get_energy_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_energy_accuracy)
}

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
NWChem::Chemistry_QC_Model_impl::set_do_energy (
  /* in */ bool doit ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_do_energy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_do_energy} (set_do_energy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_do_energy)
}

/**
 * Returns the Cartesian gradient.  
 */
::sidl::array<double>
NWChem::Chemistry_QC_Model_impl::get_gradient ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_gradient)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_gradient} (get_gradient method)
  int num;
  nwchem_getnumcoordinates_(&num);
  sidl::array<double> g = sidl::array<double>::create1d(num);
  double* gradPtr = g.first();
  nwchem_taskgradient_(gradPtr);
  return g;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_gradient)
}

/**
 * Sets the accuracy for subsequent gradient calculations
 * @param acc The new accuracy for gradients. 
 */
void
NWChem::Chemistry_QC_Model_impl::set_gradient_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_gradient_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_gradient_accuracy} (set_gradient_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_gradient_accuracy)
}

/**
 * Returns the accuracy to which the gradient is already computed.
 * The result is undefined if the gradient has not already 
 * been computed.
 * @return The current gradient accuracy. 
 */
double
NWChem::Chemistry_QC_Model_impl::get_gradient_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_gradient_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_gradient_accuracy} (get_gradient_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_gradient_accuracy)
}

/**
 * Returns the Cartesian Hessian. @return The Hessian. 
 */
::sidl::array<double>
NWChem::Chemistry_QC_Model_impl::get_hessian ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_hessian)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_hessian} (get_hessian method)
  int num;
  nwchem_getnumcoordinates_(&num);
  int lower[2]={1,1}, upper[2]={num,num};
  sidl::array<double> h = sidl::array<double>::createCol(2,lower,upper);
  double* hessPtr = h.first();
  nwchem_taskhessian_(hessPtr);
  return h;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_hessian)
}

/**
 * Sets the accuracy for subsequent Hessian calculations.
 * @param acc The new accuracy for Hessians. 
 */
void
NWChem::Chemistry_QC_Model_impl::set_hessian_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_hessian_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_hessian_accuracy} (set_hessian_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_hessian_accuracy)
}

/**
 * Returns the accuracy to which the Hessian is already computed.
 * The result is undefined if the Hessian has not already 
 * been computed. 
 */
double
NWChem::Chemistry_QC_Model_impl::get_hessian_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_hessian_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_hessian_accuracy} (get_hessian_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_hessian_accuracy)
}

/**
 * Returns a Cartesian guess Hessian. 
 */
::sidl::array<double>
NWChem::Chemistry_QC_Model_impl::get_guess_hessian ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_guess_hessian)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_guess_hessian} (get_guess_hessian method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_guess_hessian)
}

/**
 * Sets the accuracy for subsequent guess Hessian calculations.
 * @param acc The new accuracy for guess Hessians. 
 */
void
NWChem::Chemistry_QC_Model_impl::set_guess_hessian_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.set_guess_hessian_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.set_guess_hessian_accuracy} (set_guess_hessian_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.set_guess_hessian_accuracy)
}

/**
 * Returns the accuracy to which the guess Hessian is 
 * already computed.  The result is undefined if the guess Hessian 
 * has not already been computed.
 * @return The guess hessian accuracy.  
 */
double
NWChem::Chemistry_QC_Model_impl::get_guess_hessian_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.get_guess_hessian_accuracy)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.get_guess_hessian_accuracy} (get_guess_hessian_accuracy method)
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.get_guess_hessian_accuracy)
}

/**
 * This can be called when this Model object is no longer needed.  
 * No other members may be called after finalize. 
 */
int32_t
NWChem::Chemistry_QC_Model_impl::finalize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model.finalize)
  // Insert-Code-Here {NWChem.Chemistry_QC_Model.finalize} (finalize method)
  nwchem_nwchemend_();
  return 0;  
  // DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model.finalize)
}


// DO-NOT-DELETE splicer.begin(NWChem.Chemistry_QC_Model._misc)
// Insert-Code-Here {NWChem.Chemistry_QC_Model._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Chemistry_QC_Model._misc)

