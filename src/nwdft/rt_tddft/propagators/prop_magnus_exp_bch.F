C
C     prop_magnus_exp_bch.F
C
C     Steps density matrix forward in time using Baker-Campbell-Hausdorff.
C
C     P(t+dt) = e^W P(t) e^-W
C             = P + 1/1! [W,P] + 1/2![W,[W,P]] + 1/3![W,[W,[W,P]]] + ...
C      
C     where W = -iF(t+dt/2).  Note both P and W are in MO basis, but we
C     drop the prime notation for clarity.  
C     
      subroutine prop_magnus_exp_bch(params, g_zscr1, g_zscr2,
     $     g_w, g_zdens)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_w               ! W_1 = -iF(t+dt/2) dt
      integer, intent(in)           :: g_zscr1, g_zscr2  ! work space


C     == In/out ==
      integer, intent(in) :: g_zdens           !P'(t) in; P'(t+dt) out


C     == Parameters ==
      character(*), parameter     :: pname = "prop_magnus_exp_bch: "
      integer, parameter          :: kmax = 2000 !max number of BCH terms


C     == Variables ==
      double complex invfac
      double complex zinvk
      logical converged
      integer ik
      integer g_vn1             ! commutator at k-1
      logical ok
      double precision norm
      integer me
      integer num_zeroterms
      integer n
      character*100 outstring
      double precision elapsed

      if (params%prof) call prof_start (elapsed)


      me = ga_nodeid ()

      n = params%ns_mo

      ok = .true.
      ok = ok .and. ga_duplicate (g_zdens, g_vn1, "vn1")
      if (.not. ok) call errquit (pname//"alloc failed", 0, GA_ERR)


C
C     The k=0 term is just P.
C      
      ik = 0
      invfac = z1
      call ga_copy (g_zdens, g_vn1)

      
C
C     Do the summation from k=1 until kmax, or until converges.
C      
      num_zeroterms = 0
      converged = .false.
      
      do while (.not. converged)
         call ga_sync ()

         ik = ik + 1
         zinvk = z1 / dcmplx (ik)
         
C     
C     Compute 1/k*[W, vn1], where vn1 is the commutator which includes the
C     1/k prefactor from the previous step.  Then take this result and
C     store it in vn1 for use in the next k step.
C
         call ga_zero (g_zscr1)
         call ga_zero (g_zscr2)
         
         call ga_zgemm ("N", "N", n, n, n,
     $        zinvk, g_w, g_vn1, z0, g_zscr1)
         call ga_zgemm ("N", "N", n, n, n,
     $        zinvk, g_vn1, g_w, z0, g_zscr2)

         call ga_zero (g_vn1)
         call ga_add (z1, g_zscr1, zn1, g_zscr2, g_vn1)

         
C
C     Accumulate 1/k [W, vn1] (computed above) into the P(t+dt) sum.  If
C     the density matrix has not changed the sum has converged and we
C     can stop.
C     
         norm = mat_norm (g_vn1)
         
         if (norm .lt. params%tol_series)
     $        num_zeroterms = num_zeroterms + 1

         if (num_zeroterms .ge. params%terms_series)
     $        converged = .true.
         
         call ga_add (z1, g_zdens, z1, g_vn1, g_zdens)
         norm = mat_norm (g_zdens)
         
C
C     Stop if we fail to converge after kmax terms.
C
         if (ik.gt.kmax)
     $        call errquit (pname//"failed to converge", 0, 0)
      enddo

      
C
C     Clean up.
C
      ok = .true.
      ok = ok .and. ga_destroy (g_vn1)
      if (.not.ok) call errquit (pname//"failed to destroy", 0, GA_ERR)


C
C     If profiling is enabled print convergence data.
C      
      if (params%prof) then 
         write(outstring,"(a,i0,a)") "BCH converged after ",ik," terms;"
         call prof_end (elapsed, trim(outstring))
      endif
      
      end subroutine

      
c $Id$
