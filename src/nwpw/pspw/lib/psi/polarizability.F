     
*     *************************************
*     *                                   *
*     *  Calculate_Resta_Polarizability   *
*     *                                   *
*     *************************************
      subroutine Calculate_Resta_Polarizability(doprint,
     >                                  ispin,ne,neq,npack1,nfft3d,psi1,
     >                                  alpha)
      implicit none
      logical doprint                  
      integer ispin,ne(2),neq(2)
      integer npack1,nfft3d
      complex*16 psi1(npack1,*)
      real*8     alpha(3,3)
      
#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid,tmp_len
      parameter (MASTER=0,tmp_len=140)
     
      real*8 autoDebye
      parameter (autoDebye=2.5416d0)

      return
      end

