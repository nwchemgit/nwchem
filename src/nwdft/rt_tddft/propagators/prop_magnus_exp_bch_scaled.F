C
C     prop_magnus_exp_bch_scaled.F
C
C     Steps density matrix forward in time using
C     Baker-Campbell-Hausdorff.  The W matrix is scaled such that ||W||
C     <= 1, then the result is subsequently squared to recover the
C     correct answer:
C
C     P(t+dt) = e^W P(t) e^-W
C             = P + 1/1! [W,P] + 1/2![W,[W,P]] + 1/3![W,[W,[W,P]]] + ...
C      
C     where W = -iF(t+dt/2).  Note both P and W are in MO basis, but we
C     drop the prime notation for clarity.  
C
C
C     Approach
C     --------
C     (1) Choose m such that || W / 2^m || <= 1 (any norm method)
C     (2) Perform BCH expansion with W = W/2^m
C     (3) Square result m times to recover correct result
C
C     Proof:
C
C     [e^{W/y}]^y = e^W
C
C     => P(t+dt) = [ e^{W/y} P(t) e^{-W/y} ]^y
C                = e^W [ P(t) ]^y e^-W
C     
C     but P(t) (dens mat in MO basis) is idempotent: PP = P
C     
C     =>         = e^W P(t) e^{-W}    QED
C
C
C     XXX ALSO P(t+dt) is idempotent
C
C
C     XXXX WRONG BECAUSE (A B)^m /= A^m B^m
C
      subroutine prop_magnus_exp_bch_scaled (params, g_zscr1, g_zscr2,
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
      integer, intent(in) :: g_w               ! W_1 = -iF(t+dt/2) dt
      integer, intent(in) :: g_zscr1, g_zscr2  ! work space


C     == In/out ==
      integer, intent(in) :: g_zdens           !P'(t) in; P'(t+dt) out


C     == Parameters ==
      character(*), parameter :: pname="prop_magnus_exp_bch_scaled: "
      integer, parameter          :: kmax = 2000 !max number of BCH terms
      double precision, parameter :: thresh = 1d-12
      integer, parameter          :: reqd_zeroterms = 5
      integer, parameter          :: scale_power = 5


C     == Variables ==
      integer g_w_scale
      double precision scale_val
      double complex zscale
      double complex invfac
      double complex zinvk
      logical converged
      integer im
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
      ok = ok .and. ga_duplicate (g_w, g_w_scale, "W scale")
      if (.not. ok) call errquit (pname//"alloc failed", 0, GA_ERR)


C
C     First, scale the W matrix by 2^m
C
      scale_val = 1d0/2d0**scale_power
      zscale = dcmplx (scale_val)
      call ga_zero (g_w_scale)
      call ga_copy (g_w, g_w_scale)
      call ga_scale (g_w_scale, zscale)
      

      scale_val = 2d0**scale_power
      zscale = dcmplx (scale_val)
      

C
C     Start standard BCH
C

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
         ik = ik + 1


         zinvk = z1 / dcmplx (ik)
C         zinvk = z1 / dcmplx (ik) * zscale!note zscale
         
C     
C     Compute 1/k*[W, vn1], where vn1 is the commutator which includes the
C     1/k prefactor from the previous step.  Then take this result and
C     store it in vn1 for use in the next k step.
C
         call ga_zero (g_zscr1)
         call ga_zero (g_zscr2)
         

C     xxx perhaps take zinvk outside here? => clamp norm <= 1

         call ga_zgemm ("N", "N", n, n, n,
     $        zinvk, g_w_scale, g_vn1, z0, g_zscr1)
         call ga_zgemm ("N", "N", n, n, n,
     $        zinvk, g_vn1, g_w_scale, z0, g_zscr2)


         call ga_zero (g_vn1)
         call ga_add (z1, g_zscr1, zn1, g_zscr2, g_vn1)

         
C
C     Accumulate 1/k [W, vn1] (computed above) into the P(t+dt) sum.  If
C     the density matrix has not changed the sum has converged and we
C     can stop.
C     
         norm = mat_norm (g_vn1)

         if (me.eq.0) then
            write (luout, *) "norm = ", norm
            call util_flush (luout)
         endif
         
         if (norm .lt. thresh)
     $        num_zeroterms = num_zeroterms + 1

         if (num_zeroterms .ge. reqd_zeroterms)
     $        converged = .true.
         
         call ga_add (z1, g_zdens, z1, g_vn1, g_zdens)
C         norm = mat_norm (g_zdens)
         

C
C     Stop if we fail to converge after kmax terms.
C
         if (ik.gt.kmax)
     $        call errquit (pname//"failed to converge", 0, 0)
         
         call ga_sync ()
      enddo


      
C
C     Clean up.
C
      ok = .true.
      ok = ok .and. ga_destroy (g_vn1)
      ok = ok .and. ga_destroy (g_w_scale)
      if (.not.ok) call errquit (pname//"failed to destroy", 0, GA_ERR)


C
C     If profiling is enabled print convergence data.
C      
      if (params%prof) then 
         write(outstring,"(a,i0,a)") "Scaled BCH converged after "
     $        ,ik," terms;"
         call prof_end (elapsed, trim(outstring))
      endif
      
      end subroutine

      


C====================================================================
      subroutine magnus_exp_bch_scaled_new (params, g_zscr1, g_zscr2,
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
      integer, intent(in) :: g_w               ! W_1 = -iF(t+dt/2) dt
      integer, intent(in) :: g_zscr1, g_zscr2  ! work space


C     == In/out ==
      integer, intent(in) :: g_zdens           !P'(t) in; P'(t+dt) out


C     == Parameters ==
      character(*), parameter :: pname = "magnus_exp_bch_scaled_new: "
      integer, parameter      :: scalefac = 100


C     == Variables ==
      integer g_w_scale
      double complex zscalefac
      integer ik



C
C     Allocation
C     
      if (.not. ga_duplicate (g_w, g_w_scale, "W scaled"))
     $     call errquit ("failed to create W scaled", 0, 0)


      zscalefac = dcmplx (dble(1d0/scalefac), 0d0)

      call ga_zero (g_w_scale)
      call ga_copy (g_w, g_w_scale)
      call ga_scale (g_w_scale, zscalefac)


      do ik = 1, scalefac
         call magnus_exp_bch (params, g_zscr1, g_zscr2,
     $        g_w_scale, g_zdens)
      enddo
      

C
C     Clean up
C
      if (.not. ga_destroy (g_w_scale))
     $     call errquit ("failed to destroy W scaled", 0, 0)
      
      end subroutine
      
c $Id$
