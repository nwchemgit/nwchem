************************************************************************
*                                                                      *
      subroutine ecp_matpr (a,mlo,mhi,nlo,nhi,m1,m2,n1,n2,heading,fmt,
     &    nwid,ncol)
*                                                                      *
*   Routine to print portion of 2D real array with indices             *
*                                                                      *
*   Argument - description                                             *
*                                                                      *
*   a - real array to be printed                                       *
*   mlo, mhi - lower and upper bounds of first dimension of a          *
*   nlo, nhi - lower and upper bounds of second dimension of a         *
*   m1, m2 - lower and upper print limits of first dimension of a      *
*   n1, n2 - lower and upper print limits of second dimension of a     *
*   heading - character array with heading to be printed               *
*   fmt - F for F format, E for E format                               *
*   nwid - width of printer page in columns                            *
*   ncol - number of columns to be used in printing array              *
*                                                                      *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
* $Id$
      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer i,j,k,l,mlo,mhi,m1,m2,nlo,nhi,n1,n2,nwid,ncol,nfield,nn
      double precision a(mlo:mhi,nlo:nhi)
      character*(*) heading
      character*24 fmti,fmtr
      character*1 fmt
*
      write (LuOut,'(/1x,A,/)') heading
      if ((n2 .lt. n1) .or. (m2 .lt. m1)) then
        write (LuOut,*) ' **** NO DATA ****'
        return
      end if
      nfield = (nwid-6)/ncol
      l = ichar('G')-ichar(fmt)
      if ((l .lt. 1) .or. (l .gt. 2)) call errquit
     &    ('Illegal value for fmt in ecp_matpr: must be E or F',99,
     &       BASIS_ERR)
      if (nfield .lt. 4*l+1) call errquit 
     &    ('Too many columns in ecp_matpr',99, BASIS_ERR)
      write (fmti,'(A,i2,a,i2,a)') '(/X,',ncol,'I',nfield,')'
      nn = nfield-4*l
      write (fmtr,'(A,i1,A,i2,a,i2,a,i2,a)') 
     &    '(X,I3,2X,',l-1,'P',ncol,fmt,nfield,'.',nn,')'
      do i = n1,n2,ncol
        write (LuOut,fmti) (k,k = i,min(i+ncol-1,n2)) 
        write (LuOut,*)
        do j = m1,m2
          write (LuOut,fmtr) j,(a(j,k),k = i,min(i+ncol-1,n2))
        end do
      end do
      write (LuOut,*)
*
      return
      end
