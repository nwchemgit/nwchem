/* This file contains only Doxygen documentation. It does not actually
   contain any source code.
*/

/**
\defgroup nwint Gaussian integrals

\brief Notes on the Gaussian integral API

The API documented here provides access to a number of Gaussian
integral codes. At the moment NWChem contains 4 different integral
packages based on particular formulations of the Gaussian integrals
[1-3]. These packages are the SP-rotated axis code, the Rys polynomial
code [1], the McMurchie-Davidson code [2], and the Obara-Saika code [3].
Each package has its own strengths. The API provides
a unified interface to these integral packages. This documentation 
explains how to drive the integral codes and what integrals can be
calculated using them.

[1] M. Dupuis, J. Rys, H. F. King,
    "Evaluation of molecular integrals over Gaussian basis functions",
    Journal of Chemical Physics (1976), <b>65</b>, pp 111-116,
    DOI: <a href="https://doi.org/10.1063/1.432807">10.1063/1.432807</a>
    (Note: Hondo integrals using Rys quadratures).

[2] L. E. McMurchie, E. R. Davidson,
    "One- and two-electron integrals over cartesian gaussian functions",
    Journal of Computational Physics (1978), <b>26</b>, pp 218-231,
    DOI: <a href="https://doi.org/10.1016/0021-9991(78)90092-X">10.1016/0021-9991(78)90092-X</a>
    (Note: NWChem integrals using McMurchie-Davidson).

[3] S. Obara, A. Saika, 
    "Efficient recursive computation of molecular integrals over Cartesian Gaussian functions",
    Journal of Chemical Physics (1986), <b>84</b>, pp 3963-3974,
    DOI: <a href="https://doi.org/10.1063/1.450106">10.1063/1.450106</a>
    (Note: Texas integrals using Obara-Saika).
*/
/* $Id$ */
