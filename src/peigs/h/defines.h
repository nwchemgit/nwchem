*
* $Id$
*
c #define USESTDIN
#ifdef USESTDIN
#define INUNIT *
#else
#define INUNIT 5
#endif

#ifdef HERM
#ifdef TPUTER_NODES
#define DATATYPE complex
#define DATASIZE 8
#else
#define DATATYPE complex*16
#define DATASIZE 16
#endif
#define CONVERT dcmplx
#define SUMFUNC sumdc
#define INITI2 chiti2
#define INITV2 chitv2
#define SPLITV2 chspl2
#else
#define DATATYPE DoublePrecision precision
#define DATASIZE 8
#define CONVERT dble
#define SUMFUNC sumd
#define INITI2 rsiti2
#define INITV2 rsitv2
#define SPLITV2 rsspl2
#endif
#ifdef HERM
#define INTERACT1 chint1
#define INTERACT2 chint2
#define SCROLLFUNC chscr
#else
#define INTERACT1 rsint1
#define INTERACT2 rsint2
#define SCROLLFUNC rsscr
#endif
#ifdef HERM
#define KMLENGTH 4
#define DATAFUNC chdata
#define EVALFUNC cheval
#define EVECFUNC chevec
#define EVEC2FUNC chevc2
#define REDISTFUNC chrdst
#define INITFUNC chinit
#define INIT2FUNC chini2
#define JACFUNC chjac
#define JAC2FUNC chjac2
#define MEMFUNC chmem
#define MSWEEPFUNC chmswp
#define  INT1FUNC chint1
#define  INT2FUNC chint2
#else
#define KMLENGTH 3
#define DATAFUNC rsdata
#define EVALFUNC rseval
#define EVECFUNC rsevec
#define EVEC2FUNC rsevc2
#define REDISTFUNC rsrdst
#define INITFUNC rsinit
#define INIT2FUNC rsini2
#define JACFUNC rsjac
#define JAC2FUNC rsjac2
#define MEMFUNC rsmem
#define MSWEEPFUNC rsmswp
#define  INT1FUNC rsint1
#define  INT2FUNC rsint2
#endif
