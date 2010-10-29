* $Id$
      SUBROUTINE txs_SECOND(T)
      common /intgop/ ncache,maxprice,iprint,iblock
      double precision T
c     REAL*8 T
      double precision util_cpusec
      external util_cpusec
*
*     rjh ... no point timing if not printing
*
      t = 0.d0
      if (iprint .le. 0) return
*
C      DIMENSION ITIME(4) 
C     I1=TIMES(ITIME)
C     WRITE(*,*) ITIME
C     T=DBLE(ITIME(1)+ITIME(2)+ITIME(3)+ITIME(4))/100.0D0         
C      T=MCLOCK()/100.0D0
      T = util_cpusec()
C  
      RETURN
      END
C
*rak:      SUBROUTINE FDATE(X)
*rak:      CHARACTER*24 BLANK,X
*rak:      DATA BLANK /'                        '/
*rak:      X=BLANK
*rak:      END
