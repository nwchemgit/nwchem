*
* $Id$
*
*     Number of bytes in an INTEGER variable.
      INTEGER    NBYTEI
#ifdef STD_INT
      PARAMETER (NBYTEI = 4)
#else
      PARAMETER (NBYTEI = 8)
#endif
