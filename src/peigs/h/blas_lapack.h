/*
 $Id$
  c defines */

#define NO_EVEC 0

/*
  CPU definitions and machine precision definitions
  */

#ifdef ALPHA

#define DLAMCHE        2.22044604925031308e-16
#define DLAMCHP        2.22044604925031308e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.22507385850720138e-308
#define DLAMCHU        2.22507385850720138e-308
#define SLAMCHE        2.22044604925031308e-16
#define SLAMCHP        2.22044604925031308e-16
#define SLAMCHB        2.e0
#define SLAMCHS        2.22507385850720138e-308
#define SLAMCHU        2.22507385850720138e-308


/*
 Single Precision results
 depsilon      5.9604644775390625e-08
 dbase      2.0000000000000000e+00
 dsafeulp      1.1754943508222875e-38
 */
#endif

#ifdef HPPA
/*
 Double Precision results
 depsilon   1.1102230246251565e-16
 dbase      2.0000000000000000e+00
 dsafeulp   2.2250738585072013e-308
 dlamch(u)  2.2250738585072013e-308
 */
#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072013e-308
#define DLAMCHU 2.2250738585072013e-308

#endif

#ifdef IEEE
/*
 Double Precision results
 depsilon   1.1102230246251565e-16
 dbase      2.0000000000000000e+00
 dsafeulp   2.2250738585072013-308
 dlamch(u)  2.2250738585072013-308
 */
#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072013e-308
#define DLAMCHU 2.2250738585072013e-308

#endif

#ifdef SPARC

/*
  sparc
  */
/*
 Double Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014e-308

 Single Precision results
   depsilon      5.9604644775390625e-08
      dbase      2.0000000000000000e+00
   dsafeulp      1.1754943508222875e-38
*/

#define DLAMCHE        1.1102230246251565e-16
#define DLAMCHP        1.1102230246251565e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.2250738585072014e-308
#define DLAMCHU        2.2250738585072014e-308

#endif

#ifdef SPARC64
/*./teslamch
 
  Double Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014-308
 dlamch(u)       2.2250738585072014-308
 
  Single Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014-308
 slamch(u)       2.2250738585072014-308
*/
#define DLAMCHE        2.22044604925031308e-16
#define DLAMCHP        2.22044604925031308e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.22507385850720138e-308
#define DLAMCHU        2.22507385850720138e-308
/* this values were OK for WS 5.0, break for WS 6.0
#define DLAMCHE        1.1102230246251565e-16
#define DLAMCHP        1.1102230246251565e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.2250738585072014e-308
#define DLAMCHU        2.2250738585072014e-308
*/

#endif
#ifdef __crayx1
#define DLAMCHE 2.2204460492503130e-16
#define DLAMCHP 2.2204460492503130e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308

#define dscal_  sscal_
#define ddot_   sdot_
#define daxpy_  saxpy_
#define dnrm2_  snrm2_
#define dasum_  sasum_
#define dcopy_  scopy_
#define idamax_ isamax_

/*
  lapack
  */

#define dlagtf_ slagtf_
#define dlagts_ slagts_
#define dlamch_ slamch_
#define dlarnv_ slarnv_


#endif
#ifdef PENTIUM
/* wild ass guess; same as sparc */
#define DLAMCHE 2.2204460492503131e-16
#define DLAMCHP 2.2204460492503131e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308

#endif
#ifdef RS6000

/* rs6000 */

#define DLAMCHE 0.111022302462515654e-15
#define DLAMCHP 0.111022302462515654e-15
#define DLAMCHB 2.e0
#define DLAMCHS 0.22250738585072013e-307
#define DLAMCHU 0.22250738585072013e-307


/*
  depsilon  0.111022302462515654e-15 
  dbase   2.00000000000000000     
  dsafeulp  0.22250738585072013e-307 
  depsilon  0.5960464478e-07 
  dbase   2.000000000     
  dsafeulp  0.1175494351e-37 
  */
#endif
#ifdef RS600064

/* rs6000 64 -bit*/

#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.e0
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308


/*
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014e-308
 dlamch(u)       2.2250738585072014e-308
  */
#endif
#ifdef i860

/*
  with -Knoieee
  */ 

#define DLAMCHE 1.1102230246251565e-016
#define DLAMCHB 2.e0
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHE 1.1102230246251565e-016
#define DLAMCHP 1.1102230246251565e-016
#define DLAMCHB  2.e0
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHU  2.2250738585072014e-308


#endif

/* DLAMCH guesses when they are not set */
#ifndef DLAMCHE
#define DLAMCHE 2.2204460492503131e-16
#endif
#ifndef DLAMCHP
#define DLAMCHP 2.2204460492503131e-16
#endif
#ifndef DLAMCHB
#define DLAMCHB 2.0000000000000000e+00
#endif
#ifndef DLAMCHS
#define DLAMCHS 2.2250738585072013e-308
#endif
#ifndef DLAMCHU
#define DLAMCHU 2.2250738585072013e-308
#endif
