C> \brief Implementation for use with GA-4-3
C>
C> An implementation of ga_zgop to be used only in the case of compiling
C> NWChem with GA-4-3 (i.e. the environment variable OLD_GA is set). In
C> GA-5-1 and later the routine ga_zgop is part of the global array
C> distribution. In GA-4-3 this routine did not exist causing
C> linking problems. However, this routine is used in only one place
C> to sum complex numbers across processors. We can easily reproduce
C> this particular behavior using real numbers instead of complex
C> numbers. So as a fallback we will use that approach.
C>
      subroutine ga_zgop(itype,x,n,op)
      implicit none
#include "errquit.fh"
      integer itype       !< [Input] MA data type of x 
      integer n           !< [Input] length of array x
      double complex x(n) !< [In/Output] the data array
      character*(*) op    !< [Input] the operation to perform
c
      if (op.eq."+") then
        call ga_dgop(itype,x,2*n,op)
      else
        call errquit("ga_zgop: invalid operation",0,UERR)
      endif
c
      end
c $Id$
