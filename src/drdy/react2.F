C $Id: react2.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_react2 (iop)
C
C  put correct parts of hessian matrix into fsv and call rphwrt
C
      Implicit None ! Double Precision (A-H,O-Z)
C
#include "drdyP.fh"
*::passed::
      Integer iop
*::local::
      integer i, j, ixyz, jxyz
      integer indx1, indx2, jndx1, jndx2
C
      
      indx1 = 0
      do i = 1,nratom(iop)
         indx2 = 3*(iatsv(i,iop)-1)
         do ixyz = 1,3
            indx1 = indx1 + 1
            indx2 = indx2 + 1
            jndx1 = 0
            do j = 1,nratom(iop)
               jndx2 = 3*(iatsv(j,iop)-1)
               do jxyz = 1,3
                  jndx1 = jndx1 + 1
                  jndx2 = jndx2 + 1
                  fsv(jndx1,indx1) = f(jndx2,indx2)
               enddo
            enddo
         enddo
      enddo
      call drdy_rphwrt(iop)
      return
      end
