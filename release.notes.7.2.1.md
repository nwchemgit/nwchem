NWChem Version 7.2.1 Release Notes
==================================

NWChem is now available on Github at
https://github.com/nwchemgit/nwchem

User Manual available from the NWChem website
https://nwchemgit.github.io

NWChem 7.2.1 is released as open-source under the ECL 2.0 license.

NWChem 7.2.1 will be released with the latest Global Arrays Toolkit (v5.8.2).

The  7.2.1 release is a maintenance release containing fixes/enhancements for the NWChem 7.2.0 tree

NEW FUNCTIONALITY
-----

   N/A
   

BUG FIXES/ENHANCEMENTS
-----

* fix compilation for Cray compilers

* fixes for macOS Xcode 15

* fixes for Intel OneAPI 2023 releases

* memory fixes for ARMCI_NETWORK=MPI-PR https://github.com/GlobalArrays/ga/pull/310

GITHUB ISSUES ADDRESSED
----

   * [Shifter Image Parallelization Issues](https://github.com/nwchemgit/nwchem/issues/775)
   * [Custom dielectric constant is ignored in COSMO calculations](https://github.com/nwchemgit/nwchem/issues/776)
   * [NMR hyperfine coupling when wavefunction is closed-shell](https://github.com/nwchemgit/nwchem/issues/788)
   * [Triplet CDSpectrum calculations failing](https://github.com/nwchemgit/nwchem/issues/796)
   * [7.2.0 fails computing gradients of bare ECPs](https://github.com/nwchemgit/nwchem/issues/801)
   * [SegV when USE_INTERNALBLAS=1 for version 7.2.0](https://github.com/nwchemgit/nwchem/issues/804)
   * [Undefined symbol "ycnrm2"](https://github.com/nwchemgit/nwchem/issues/817)
   * [Wrong symmetry assignments in TDDFT](https://github.com/nwchemgit/nwchem/issues/828)
   * [GW might lead to incorrect results when symmetry is on](https://github.com/nwchemgit/nwchem/issues/829)
   * [ch5n_nbo QA test fails](https://github.com/nwchemgit/nwchem/issues/864)
   
