C
C     exp_pseries.F
C
C     Brute force exponentiation of a complex operator A via a power
C     series.  Works for any square complex matrix.
C
C     e^A = \sum_{k=0}^{\infty} 1/k! A^k
C
C     We first scale A by 2^m (such that ||A|| <= 1), then square the
C     result m times to recover the answer:
C
C     [ e^{A/2^m} ]^{2^m} = e^A
C     
C     We compute the necessary m each time via the infinity norm of ||A||.
C     
      subroutine exp_pseries (params, g_za, g_zexpa)
      implicit none


#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
      

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_za 

      
C     == Outputs ==
      integer, intent(in) :: g_zexpa


C     == Parameters ==
      character(*), parameter   :: pname = "exp_pseries: "
      integer, parameter        :: maxk = 1000 !max num terms in expansion
      integer, parameter        :: max_mscale = 16 !=65536, norm should never be this big

      
C     == Variables ==
      integer ik
      integer dtype, n1, n2
      integer g_za_scaled
      integer im
      double precision scale_val
      double complex zscale
      integer g_new_term, g_prev_term
      double complex invfac
      integer me
      logical converged
      integer num_zeroterms
      character*120 outstring
      double precision elapsed
      double precision invk
      double complex zinvk
      double precision norm
      double precision ratio
      integer mscale            !divide A by 2^m, then square result m times


      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid ()

      
C
C     Check the GA.
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//" not a valid GA")
      
      call ga_inquire (g_za, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")

      if (.not. ga_duplicate(g_za, g_za_scaled, "g_za_scaled"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)
      if (.not. ga_duplicate(g_za, g_new_term, "g_new_term"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)
      if (.not. ga_duplicate(g_za, g_prev_term, "g_prev_term"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)



C
C     Determine 2^m factor and scale input matrix.
C
C     
C     We technically need to compute the scaling factor from the
C     spectral range of the W operator, but a simple norm will suffice.
C
C     If L is the largest abs eigenvalue (or the norm as we do it
C     below), then we pick scaling factor "m" such that:
C
C            2^m = 2.5*L
C     => m log 2 = log (2.5*L)
C              m = log (2.5*L) / log 2
C
C     (we multiply by 2.5 to ensure ||W|| < 1)
C
C
      call ga_norm_infinity (g_za, norm)

      ratio = dlog (2.5d0*norm) / dlog (2d0)
      mscale = max (int (ratio), 0)


      if (mscale .lt. 0)
     $     call errquit (pname//"negative mscale", 0, 0)

      if (mscale .gt. max_mscale)
     $     call errquit (pname//"mscale exceeds 16 (2^16 = 65536)",
     $     0, 0)


      scale_val = 1d0/(2d0**mscale)
      zscale = dcmplx (scale_val, 0d0)

      call ga_zero (g_za_scaled)
      call ga_copy (g_za, g_za_scaled)
      call ga_scale (g_za_scaled, zscale)


C
C     k = 0 term is just I
C
      call ga_zero (g_prev_term)
      call mat_set_ident (g_prev_term)

      call ga_zero (g_zexpa)
      call ga_copy (g_prev_term, g_zexpa)

      ik = 0
      converged = .false.
      num_zeroterms = 0
      
      do while (.not. converged)
         call ga_sync ()

         ik = ik + 1
         zinvk = z1 / dcmplx (ik)
         
         
C
C     New term is 1/k A times previous term
C
         call ga_zero (g_new_term)
         call ga_zgemm ("N", "N", n1, n1, n1,
     $        zinvk, g_prev_term, g_za_scaled, z0, g_new_term)

         

C
C     Add new term to sum and store new term as previous.
C
         call ga_add (z1, g_zexpa, z1, g_new_term, g_zexpa)
         call ga_zero (g_prev_term)
         call ga_copy (g_new_term, g_prev_term)
         
         
C
C     Check if we have converged
C
         norm = mat_norm (g_new_term)
 
         if (norm < params%tol_series) num_zeroterms = num_zeroterms + 1

         if (num_zeroterms .ge. params%terms_series)
     $        converged = .true.


C
C     Stop if we have failed to converge after max num terms.
C
         if (ik.gt.maxk)
     $        call errquit (pname//"failed to converge", 0, 0)

      enddo


C
C     Square output matrix m times to recover unscaled result.  Use
C     g_prev_term as scratch space.
C
      do im = 1, mscale
         call ga_sync ()
         call ga_zero (g_prev_term)
         call ga_zgemm ("N", "N", n1, n1, n1,
     $        z1, g_zexpa, g_zexpa, z0, g_prev_term)
         call ga_zero (g_zexpa)
         call ga_copy (g_prev_term, g_zexpa)
      enddo

      
C
C     Clean up
C
      if (.not. ga_destroy(g_new_term))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      if (.not. ga_destroy(g_prev_term))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      if (.not. ga_destroy(g_za_scaled))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      
C
C     If profiling is enabled print convergence data.
C      
      if (params%prof) then 
         write(outstring,"(a,i0,a,i0,a)")
     $        "Power series with scaling 2^",
     $        mscale, " converged after ", ik," terms;"
         call prof_end (elapsed, trim(outstring))
      endif

      end subroutine

c $Id$
