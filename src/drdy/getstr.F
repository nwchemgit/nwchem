C $Id$
      Subroutine drdy_getstr(xnum,cstr,nstr,nten,ndec)
C
C   get string of numbers for a floating point number xnum
C      ndec = digits after decimal point (must be >0)
C      nten = digits before decimal point
C   on return nstr should equal nten if ndec=0 and ndec+nten+1 for ndec>0
C 
      Implicit None
#include "errquit.fh"
C
      Integer nstr,nten,ndec
      Character*1 cstr(1)
      Double Precision xnum
C
      Integer ndim
      Parameter (ndim=20)
      Integer i,idig,idigx,int1,int2,ndig
      Integer ldig(ndim)
      Character*1 cnum(10)
      Double Precision r,r1,r2
C
      Data cnum/'0','1','2','3','4','5','6','7','8','9'/
      Integer fu6
      Data fu6/6/
C      
      if (nten+ndec.ge.ndim) then
         write(fu6,*) ' in getstr, nten,ndec=',nten,ndec
         write(fu6,*) ' but their sum must be <',ndim
         call errquit('drdy_getstr: problem',555, UNKNOWN_ERR)
      else if (nten.le.0) then
         write(fu6,*) ' in getstr, nten=',nten
         write(fu6,*) ' but they must be > 0'
         call errquit('drdy_getstr: problem',555, UNKNOWN_ERR)
      else if (ndec.lt.0) then
         write(fu6,*) ' in getstr, ndec=',ndec
         write(fu6,*) ' but they must be >= 0'
         call errquit('drdy_getstr: problem',555, UNKNOWN_ERR)
      endif
      ndig = 0
      r = xnum*(0.1d00**nten)
      int2 = int(r)
      r2 = float(int2)
      do i = 1,nten
         r = r*10.0d00
         int1 = int(r)
         r1 = float(int1)
         ndig = ndig + 1
         ldig(ndig) = int1 - int(10.0d00*r2)
         int2 = int1
         r2 = r1
      enddo
      if (ndec.gt.0) then
         r = xnum
         int2 = int(r)
         r2 = float(int2)
         do i = 1,ndec
            r = r*10.0d00
            int1 = int(r)
            r1 = float(int1)
            ndig = ndig + 1
            ldig(ndig) = int1 - int(10.0d00*r2)
            int2 = int1
            r2 = r1
         enddo
      endif
      r = r*10.0d00
      int1 = int(r)
      r1 = float(int1)
      idigx = int1 - int(10.0d00*r2)
      if (idigx.ge.5) then
         ldig(ndig) = ldig(ndig)+1
         idig = ndig+1
         do i = 1, ndig-1
            idig = idig - 1
            if (ldig(idig).gt.9) then
               ldig(idig) = 0
               ldig(idig-1) = ldig(idig-1) + 1
            endif
         enddo
         if(ldig(1).gt.9) ldig(1) = 0
      endif
      nstr = 0
      do idig = 1,nten
         nstr = nstr + 1
         cstr(nstr) = cnum(ldig(idig)+1)
      enddo
      if (ndec.gt.0) then
         nstr = nstr + 1
         cstr(nstr) = '.'
         do idig = 1,ndec
            nstr = nstr + 1
            cstr(nstr) = cnum(ldig(nten+idig)+1)
         enddo
      endif
      return
      end
