/*
 $Id$
*/
#ifndef _DFT_H_
#define _DFT_H_
/* dft.h - 6/9/95
   author - Eric Bylaska

   This file contains a routine for finding the hartree
   potential and energy of a density rho, defined on
   a log grid.

*/
#include	"hartree.h"
#include	"dirac_exchange.h"
#include        "pbe_exchange.h"
#include        "revpbe_exchange.h"
#include	"perdew_zunger.h"
#include        "perdew_wang.h"
#include	"vosko.h"
#include        "pbe_correlation.h"
#include        "revpbe_correlation.h"
#include        "becke_exchange.h"
#include        "lyp_correlation.h"


extern	void	init_DFT(char*);
extern	double	R_Hartree_DFT();
extern	void	R_Exchange_DFT();
extern	void	R_Correlation_DFT();
extern  void    R_Screening_Cut();



/* Hartree type: Hartree_Type */
#define Hartree_On	-8101
#define Hartree_Off	-8102


/* Exchange type: Exchange_Type */
#define Exchange_Dirac          -8201
#define Exchange_PBE96          -8202
#define Exchange_Becke          -8204
#define Exchange_revPBE         -8205
#define Exchange_Off            -8206


/* Correlation type: Correlation_Type */
#define Correlation_Vosko               -8301
#define Correlation_Perdew_Zunger       -8302
#define Correlation_Perdew_Wang         -8303
#define Correlation_PBE96               -8304
#define Correlation_LYP                 -8305
#define Correlation_revPBE              -8306
#define Correlation_Off                 -8307


/* used for setting solver parameters */
extern void	set_Hartree_DFT();
extern void	set_Exchange_DFT();
extern void	set_Correlation_DFT();

extern	char	*hartree_Name_DFT();
extern	char	*exchange_Name_DFT();
extern	char	*correlation_Name_DFT();

#endif
