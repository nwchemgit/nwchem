*
* $Id: cpsi_KS.F,v 1.1 2003-12-02 19:16:47 bylaska Exp $
*
*
*
* $Log: not supported by cvs2svn $



************************ KS orbital Part ************************

*     ***********************************
*     *                                 *
*     *      psi_KS_update              *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine psi_KS_update(psi_number,precondition,maxerror)
      implicit none
      integer psi_number
      logical precondition
      real*8 maxerror

#include "mafdecls.fh"
#include "cpsi.fh"


*     **** local variables ****
      logical done
      integer i,j,neall,maxit_orb,maxit_orbs
      real*8 error,error_out,tim1,tim2,tim,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      neall = ne(1)+ne(2)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        do i=1,neall
          call current_second(tim1)

         !*** orthogonalize to lower orbitals  ****
         call cpsi_project_out_f_orb1(
     >           i,
     >           dcpl_mb(psi1(1)+(i-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
         call Pack_c_SMul(1,sum,
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            dcpl_mb(psi1(1) +(i-1)*npack1))



          call cpsi_KS_update_orb(psi_number,precondition,maxit_orb,
     >                         maxerror,
     >                         0.1d0,i,error_out)
          error = error+error_out
        end do
        error = error/dble(neall)

        done = ((j.gt.maxit_orbs).or.(error.lt.maxerror))
      if (.not.done) go to 2

      return
      end
