*
* $Id$
*
*     *******************************
*     *								*
*     *		   lcao_make_psi		*
*     *								*
*     *******************************
*
*     This routines calculates the electronic wavefunction.
*
      subroutine lcao_make_psi(coef,psi,phi)
      implicit none
      real*8 coef(*)
      complex*16 psi(*)
      complex*16 phi(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      integer j,nbasis


*     ***** external functions ****
      integer  aorbs_nbasis
      external aorbs_nbasis


      nbasis = aorbs_nbasis()
      call Pack_c_Zero(1,psi)
      do j=1,nbasis
         call aorbs_aorb(j,phi)  !**** get phi1 ****
         call Pack_cc_daxpy(1,coef(j),phi,psi)
      end do

      return
      end

