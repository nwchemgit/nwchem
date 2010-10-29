************************************************************************
*                                                                      *
      subroutine ecp_matpi (a,mlo,mhi,nlo,nhi,m1,m2,n1,n2,heading,
     &    nwid,ncol)
*                                                                      *
*   Routine to print portion of 2D integer array with indices          *
*                                                                      *
*   Argument - description                                             *
*                                                                      *
*   a - integer array to be printed                                    *
*   mlo, mhi - lower and upper bounds of first dimension of a          *
*   nlo, nhi - lower and upper bounds of second dimension of a         *
*   m1, m2 - lower and upper print limits of first dimension of a      *
*   n1, n2 - lower and upper print limits of second dimension of a     *
*   heading - character array with heading to be printed               *
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
#include "util.fh"
#include "errquit.fh"
      integer i,j,k,mlo,mhi,m1,m2,nlo,nhi,n1,n2,nwid,ncol,nfield
      integer a(mlo:mhi,nlo:nhi)
      character*(*) heading
      character*24 fmti,fmtr
*
      write (LuOut,'(/1X,A,/)') heading
      if ((n2 .lt. n1) .or. (m2 .lt. m1)) then
        write (LuOut,*) ' **** NO DATA ****'
        return
      end if
      nfield = (nwid-6)/ncol
      if (nfield .lt. 2) call errquit ('Too many columns in ecp_matpi',
     &    99, UNKNOWN_ERR)
      write (fmti,'(A,i2,a,i2,a)') '(/6X,',ncol,'I',nfield,')'
      write (fmtr,'(A,i2,a,i2,a)') 
     &    '(X,I3,2X,',ncol,'I',nfield,')'
C      write (LuOut,*) fmti
C      write (LuOut,*) fmtr
      do i = n1,n2,ncol
        write (LuOut,fmti) (k,k = i,min(i+ncol-1,n2)) 
        write (LuOut,*)
        do j = m1,m2
          write (LuOut,fmtr) j,(a(j,k),k = i,min(i+ncol-1,n2))
        end do
      end do
*
      return
      end
