#ifndef _OUTPUT_
#define _OUTPUT_

#define E_OUT_OK   0	/* Successful completion */
#define E_OUT_DIM  1	/* Error in matrix dimensions */
#define E_OUT_FMT  3	/* Error in element format */
#define E_OUT_WID  5	/* Can't fit at least one element per line */
#define E_OUT_UPLO 7	/* Error specifying argument UPLO */

#define OUT_NUM     1   /* Label only with numbers */
#define OUT_TXT     2   /* Label only with text */
#define OUT_NUM_TXT 3   /* Label with both number and text */

/* #include "fp_precision.h"  /* Make sure working precision is defined */

#ifdef SNGLPR
#   define DBGEWR	SBGEWR
#   define DBSPWR       SBSPWR
#   define DBSYWR	SBSYWR
#   define DBVWR	SBVWR
#   define DOSYWR       SOSYWR
#   define DGEWR	SGEWR
#   define DGEWRL	SGEWRL
#   define DGEWR2	SGEWR2
#   define DSYWR	SSYWR
#   define DSPWR        SSPWR
#endif SNGLPR

#endif !_OUTPUT_
