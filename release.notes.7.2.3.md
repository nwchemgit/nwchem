NWChem Version 7.2.3 Release Notes
==================================

NWChem is now available on Github at
https://github.com/nwchemgit/nwchem

User Manual available from the NWChem website
https://nwchemgit.github.io

NWChem 7.2.3 is released as open-source under the ECL 2.0 license.

NWChem 7.2.3 will be released with the latest Global Arrays Toolkit (v5.8.2).

The  7.2.3 release is a maintenance release containing fixes/enhancements for the NWChem 7.2.0 tree

The change log below is relative to the 7.2.2 code base.

NEW FUNCTIONALITY
-----

   N/A

BUG FIXES/ENHANCEMENTS
-----

   * added code to deal with elements up to Z=120
   * added IMOM (Initial Maximum Overlap Method)
   * reworked makefile structure to compile with GNU make 4.4
   * fixes for Python interfaces
   * avoid need of USE_MPI for make clean
   * compiler updates
   * bug fix for VEM in TDDFT
   * updates for IBO and Pipek-Mezey localization
   

GITHUB ISSUES ADDRESSED
----
   * [Cosmo generates NaN when lineq=1](https://github.com/nwchemgit/nwchem/issues/990)
   * [dft-3d URL update](https://github.com/nwchemgit/nwchem/issues/962)
   * [x2c incompatible with cd fitting ](https://github.com/nwchemgit/nwchem/issues/931)
   * [wrong memory setup with floating point input ](https://github.com/nwchemgit/nwchem/issues/930)
   * [modelpotential input not processing bq elements](https://github.com/nwchemgit/nwchem/issues/926)
   * [memory printout cannot handle large memory sizes](https://github.com/nwchemgit/nwchem/issues/838)
   * [O(10^3) global arrays are created from the property code](https://github.com/nwchemgit/nwchem/issues/831)
   
