*
* $Id: datasizes.h,v 1.2 1999-07-28 00:39:02 d3e129 Exp $
*
*     Number of bytes in an INTEGER variable.
      INTEGER    NBYTEI
#ifdef STD_INT
      PARAMETER (NBYTEI = 4)
#else
      PARAMETER (NBYTEI = 8)
#endif
