C
C $Id: twodd_cont.F,v 1.1 1999-05-06 17:29:35 windus Exp $
C
      subroutine twodd_cont(rtdb)
c
c  Add the two electron derivative (first and seconds) contributions to 
c  the hessian matrix and the RHS.
c
      implicit none
#include "hess_info.fh"
#include "bgj.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb, i
      double precision alpha
c
c  Coulomb and exact exchange derivative contributions
c  For now alpha is set to one which give a full exchange contribution
c
      alpha = 1.0
c     call twodd_coul_ex(.true., .true., alpha)
c
c  For now we send in seperate memory, k_hess_xc and g_rhs_xc, to make sure
c  that the xc contributions are coming out O.K.  Eventually, the xc routines
c  will need to modified to add their contribution to k_exy and g_rhs and
c  not zero out the matrices when they enter the routines.
c  I will copy the results into the actual "permanent" storage for now so
c  the rest of the code will work properly.
c
      if (theory.eq.'DFT') then
        if (bgj_print() .gt. 0)
     *    write(6,*)'*** In twodd_cont: calling xc_hessian'
        call xc_hessian(g_dens, log_mb(k_act), nactive, 
     *     dbl_mb(k_hess_xc))
        if (bgj_print() .gt. 0)
     *     write(6,*)'*** In twodd_cont: calling xc_cpks_rhs'
        call xc_cpks_rhs(g_dens, log_mb(k_act), nactive, g_rhs_xc)
      endif

c
c  Temporary copy to real spot (needs to be an add in the very near future)
c
      do i = 1,ncent
        call ga_copy(g_rhs_xc(i), g_rhs(i))
      enddo

      call dcopy(n3xyz2,dbl_mb(k_hess_xc),1,dbl_mb(k_exy),1)
c
      end
