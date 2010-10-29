C $Id$
c ***************************************************************************
c     icint
c ***************************************************************************
c
c This routine is no longer used!  It is now part of io.F.
c
      integer function drdy_icint(svalue)
c
c   Function which takes a character string which is really a integer and
c   converts it to an integer.  Error checking not possible using ichar is
c   done.  Program will stop if an illegal character is passed via the 
c   string 'svalue'.
c
      implicit double precision (a-h,o-z)
#include "errquit.fh"
c
      character * 80 svalue
c
c     check for negative or positive sign (- or +)
c
      drdy_icint = 0
      istart = 1
      isign = +1
      if (svalue(1:1).eq.'+') then
         isign = +1
         istart = 2
      else if(svalue(1:1).eq.'-') then
         isign = -1
         istart = 2
      end if
c
c     now convert each character to an integer - and store in drdy_icint
c
      k = istart
      do while (svalue(k:k).ne.' ')
         j = ichar(svalue(k:k)) - ichar('0')
         if (j.lt.0.or.j.gt.9) then
            write(6,1000)svalue
            call errquit('drdy_icint: fatal error',911, UNKNOWN_ERR)
         end if
1000     format(3x,'ERROR: illegal digit in an integer: ',A80)
c
         drdy_icint = 10.d00*drdy_icint + dble(j)
         k = k+1
      end do
      drdy_icint = drdy_icint * isign
c
      return
      end
