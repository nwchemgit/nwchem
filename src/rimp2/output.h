#ifndef _OUTPUT_
#define _OUTPUT_

#define E_OUT_OK   0
#define E_OUT_DIM  1
#define E_OUT_FMT  3
#define E_OUT_WID  5
#define E_OUT_UPLO 7

#define OUT_NUM     1
#define OUT_TXT     2
#define OUT_NUM_TXT 3

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
#endif

#endif
/* $Id$ */
