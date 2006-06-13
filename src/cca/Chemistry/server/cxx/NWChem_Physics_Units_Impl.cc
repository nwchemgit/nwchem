// 
// File:          NWChem_Physics_Units_Impl.cc
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
#include "NWChem_Physics_Units_Impl.hh"

// DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._includes)
// Insert-Code-Here {NWChem.Physics_Units._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(NWChem.Physics_Units._includes)

// user-defined constructor.
void NWChem::Physics_Units_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._ctor)
  // Insert-Code-Here {NWChem.Physics_Units._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units._ctor)
}

// user-defined destructor.
void NWChem::Physics_Units_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._dtor)
  // Insert-Code-Here {NWChem.Physics_Units._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units._dtor)
}

// static class initializer.
void NWChem::Physics_Units_impl::_load() {
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._load)
  // Insert-Code-Here {NWChem.Physics_Units._load} (class initialization)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Initializes the units as a human readable string
 * options are "angstroms" or "bohr" 
 */
void
NWChem::Physics_Units_impl::initialize (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units.initialize)
  // Insert-Code-Here {NWChem.Physics_Units.initialize} (initialize method)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units.initialize)
}

/**
 * Returns the units as a human readable string. 
 */
::std::string
NWChem::Physics_Units_impl::get_unit_name ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units.get_unit_name)
  // Insert-Code-Here {NWChem.Physics_Units.get_unit_name} (get_unit_name method)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units.get_unit_name)
}

/**
 * Converts from self's units to the given unit name. 
 */
double
NWChem::Physics_Units_impl::convert_to (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units.convert_to)
  // Insert-Code-Here {NWChem.Physics_Units.convert_to} (convert_to method)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units.convert_to)
}

/**
 * Converts to self's units from the given unit name. 
 */
double
NWChem::Physics_Units_impl::convert_from (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(NWChem.Physics_Units.convert_from)
  // Insert-Code-Here {NWChem.Physics_Units.convert_from} (convert_from method)
  // DO-NOT-DELETE splicer.end(NWChem.Physics_Units.convert_from)
}


// DO-NOT-DELETE splicer.begin(NWChem.Physics_Units._misc)
// Insert-Code-Here {NWChem.Physics_Units._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(NWChem.Physics_Units._misc)

