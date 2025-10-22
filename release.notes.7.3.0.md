NWChem Version 7.3.0 Release Notes
==================================

NWChem is now available on Github at
https://github.com/nwchemgit/nwchem

User Manual available from the NWChem website
https://nwchemgit.github.io

NWChem 7.3.0 is released as open-source under the ECL 2.0 license.

NWChem 7.3.0 will be released with the latest Global Arrays Toolkit (v5.9.2).

The change log below is relative to the 7.3.0 code base.

NEW FUNCTIONALITY
-----

   * Solvation: new cavity construction approach based on the solvent-excluding surface (SES),
     using the GEPOL algorithm
   * Solvation: added .cosmo file generation to be used by COSMO-RS and COSMO-SAC 
   * Bethe-Salpeter Equation

BUG FIXES/ENHANCEMENTS
-----

   * port to GNU 15 compiler
   * port to LLVM Flang 21 compiler
   * port to Intel OneAPI 2025 compiler [PR#1116](https://github.com/nwchemgit/nwchem/pull/1116)
   * enhancement to the socket driver [PR#1145](https://github.com/nwchemgit/nwchem/pull/1145)
   * removed dependence from PeIGS library (that can only be installed when USE_PEIGS=Y)
   * IBO code improvements [PR#1065](https://github.com/nwchemgit/nwchem/pull/1065)
   
   
   

GITHUB ISSUES ADDRESSED
----
   * [RT-TDDFT treference keyword not working](https://github.com/nwchemgit/nwchem/issues/1160)
   * [BUG: NWChem socket driver exits in a single-shot](https://github.com/nwchemgit/nwchem/issues/1144)
   * [Bug in fitted Coulomb when geometry has dummy atoms](https://github.com/nwchemgit/nwchem/issues/1139)
   * [Bug in libxc interface while using a metaGGA functional and cgmin](https://github.com/nwchemgit/nwchem/issues/1137)
   * [Multiple task operations in single command line](https://github.com/nwchemgit/nwchem/issues/1119)
   * [Crash when geometry name with more than 63 characters](https://github.com/nwchemgit/nwchem/issues/1106)
   * [Solvation: egas variable not initialized when do_gasphase=false](https://github.com/nwchemgit/nwchem/issues/1101)
   * [rt_tddft: data corruption when more than 100 geometries are used](https://github.com/nwchemgit/nwchem/issues/1075)
   * [r2SCAN implementation error](https://github.com/nwchemgit/nwchem/issues/1067)
   * [Problem with downloading DFT-D3 library](https://github.com/nwchemgit/nwchem/issues/1053)
   * [task bsse not compatible with sodft](https://github.com/nwchemgit/nwchem/issues/1049)
   * [ECP atom label incorrectly applied](https://github.com/nwchemgit/nwchem/issues/1037)
   * [Dead URL in contrib/getfiles.nwchem for dftd3.tgz](https://github.com/nwchemgit/nwchem/issues/1011)
   * [pspw_md QA fails on arm64 MacOS with gfortran-12 and clang-15](https://github.com/nwchemgit/nwchem/issues/970)
