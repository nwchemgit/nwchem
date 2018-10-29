Using NWChem with Microsoft Quantum Development Kit
===================================================

Instructions on building and running NWChem are available
[here](http://www.nwchem-sw.org/index.php/Main_Page).

Alternatively, you can run the NWChem docker image using the
instructions at:

<to be filled in>

The general format of NWChem input files is described
[here](https://github.com/nwchemgit/nwchem/wiki). You can also start
with any of the input files in `QA/chem_library_tests` in the NWChem
repository.

## NWChem Input File Format

To generate inputs to Microsoft QDK, the NWChem input file currently
needs to have the following components:

* `echo` directive to capture geometry and basis set information

* `set tce:print_integrals T` directive to obtain one and two electron integrals

* `set tce:qorb <#orbitals>` directive specifying number of orbitals

* `set tce:qela  <#alpha electrons>` directive specifying number of alpha electrons

* `set tce:qelb  <#beta electrons>` directive specifying number of alpha electrons

* `tce` block with `ccsd` directive is needed to obtain initial state
  suggestions and reference energies

mcscf can be used to extract representative FCI energies. In the
absence of mcscf, tce ccsd energies will be used instead.

## Generating QDK Input

The output file from running NWChem can be used to generate input to
Microsort QDK using the `export_chem_library_yaml.py` script in this
directory as:

`$python3 export_chem_library_yaml.py < <nwchem_output_file> > <yaml_output>`

