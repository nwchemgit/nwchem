C
C     prop_magnus_ord2_interpol.F
C
C     Propagate complex density matrix forward via Magnus with
C     interpolated and converged Fock matrix (very accurate but
C     potentially involves many Fock builds per time step).
C
C     This routine works for closed shell, open shell, spin-orbit, etc
C     by automatically propagating each of the density matrix GA's
C     (e.g., you have one for CS or SO, two or OS) and using the
C     supplied external "tdfock" function to compute F(t).
C
C
C     Overall approach
C     ----------------
C      
C     (1) Extrapolate F'(t+1/2*dt) from F'(t-1/2*dt) and F'(t)
C     (2) Propagate P'(t) -> P'(t+dt) using extrapolated F'(t+1/2*dt); store
C     (3) Convert P'(t+dt)->P(t+dt); build F(t+dt) from P(t+dt)
C     (4) Convert F(t+dt)->F'(t+dt); interpolate new F'(t+1/2*dt)
C     (5) Propagate P'(t)->P'(t+dt) using new F'(t+1/2*dt) from (4)
C     (6) If new P'(t+1/2*dt) is same as previous one exit, else go to (3)
C
C
C     Extrapolating F'
C     ----------------
C
C     Given F'(t-1/2*dt) and F'(t) predict F'(t+1/2*dt):
C
C     F'(t+1/2*dt) = F'(t) + m*1/2*dt,
C
C     where m is the slope between the two previous values:
C
C     m = [F'(t) - F'(t-dt)]/(1/2*dt)
C
C     Therefore,
C
C     F'(t+1/2*dt) = F'(t) + [F'(t) - F'(t-1/2*dt)]/(1/2*dt)*(1/2*dt)
C                = F'(t) + [F'(t) - F'(t-1/2*dt)]
C                = 2*F'(t) - F(t-1/2*dt)
C     
C     Similarly, 
C     
C     F'(t+dt) = F'(t) + m * dt
C              = F'(t) + [F'(t) - F'(t-dt)]/(1/2*dt)*dt
C              = F'(t) + 2*[F'(t) - F'(t-1/2*dt)]
C              = 3*F'(t) - 2*F'(t-1/2*dt)
C     (not used)
C
C
C     Propagating P'
C     --------------
C     P'(t+dt) = e^W P'(t) e^-W,
C     
C     where W = Omega_1 = -i F'(t+1/2*dt)  (2nd order Magnus).
C
C
C     
C     Interpolating F'
C     ---------------
C
C     Given F'(t) and F'(t+dt) interpolate to get F'(t+1/2*dt):
C
C     F'(t+1/2*dt) = F'(t) + m * dt/2,
C
C     where m is the slope:
C
C     m = [ F'(t+dt) - F'(t) ] / dt
C
C     Therefore,
C
C     F'(t+1/2*dt) = F'(t) + [ F'(t+dt) - F'(t) ] / dt * (dt/2)
C                  = F'(t) + 1/2* [ F'(t+dt) - F'(t) ]
C                  = 1/2*F'(t+dt) + 1/2*F'(t)
C
C
      subroutine prop_magnus_ord2_interpol (params, tt, tdfock,
     $     nmats, g_zfock_n12dt, g_zfock, energies, g_zdens_mo)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt      !current time
      logical, external             :: tdfock  !external function for computing F(t)
      integer, intent(in)           :: nmats   !num of matrices (1 for CS and SO, 2 for OS)


C     == In/out ==
      integer, intent(in) :: g_zdens_mo(nmats)        ! in = P'(t); out = P'(t+dt)
      integer, intent(in) :: g_zfock_n12dt(nmats)     ! F'(t-1/2*dt) (MO basis); update at end
      integer, intent(in) :: g_zfock(nmats)           ! F'(t) (MO basis); update at end

      
C     == Outputs ==
      type(rt_energies_t) energies


C     == Parameters ==
      character(*), parameter     :: pname =
     $     "prop_magnus_ord2_interpol: "
      integer, parameter          :: ninter_max = 200
      

C     == Variables ==
      integer me
      logical ok
      integer nmats2, nmats3
      integer imat
      character amat(nmats)

      logical converged
      integer iinter
      double precision elapsed
      double precision ttdt     ! t + dt
      character*256 outstring
      integer num_same

C     (simple way: have two GA handles and allocate one or two arrays, depending)
      integer g_zdens(nmats)           ! P'(t), alias for g_zdens_mo 
      integer g_zdens_pdt(nmats)       ! P'(t+dt)
      integer g_zdens_pdt_old(nmats)   ! P'(t+dt) from previous iteration
      integer g_zfock_p12dt(nmats)     ! F'(t+1/2*dt); first extrapolated, then interpolated)
      integer g_zfock_pdt(nmats)       ! F'(t+dt); built from P'(t+dt)
      integer g_zfock_ao(nmats), g_zdens_ao(nmats)

      double precision norm1(nmats), norm2(nmats), diff(nmats)
      double complex zjunk

      logical lskip_interpol

      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid ()


      if (.not. rtdb_get (params%rtdb, "rt_tddft:magnus_skipinterpol",
     $     mt_log, 1, lskip_interpol))
     $     lskip_interpol = .false.
      
      
C     
C     Determine number of matricies and check consistency.
C
      if ( (nmats .ne. 1).and.(nmats .ne. 2) )
     $     call errquit (pname//"Only works for 1 or 2 matrices",0,0)

      do imat = 1, nmats
         write (amat(imat), "(i1)") imat !assumes single digit number of mats (ok since max 2)
      enddo

      
C
C     Allocation
C
      ok = .true.

      do imat = 1, nmats
         ok = ok .and. ga_duplicate (g_zdens_mo(imat),
     $        g_zdens_pdt(imat), "P'(t+dt) "//amat(imat))

         ok = ok .and. ga_duplicate (g_zdens_mo(imat),
     $        g_zdens_pdt_old(imat), "P'(t+dt) old "//amat(imat))

         ok = ok .and. ga_duplicate (g_zfock(imat),
     $        g_zfock_p12dt(imat), "F'(t+1/2dt) "//amat(imat))

         ok = ok .and. ga_duplicate (g_zfock(imat),
     $        g_zfock_pdt(imat), "F'(t+dt) "//amat(imat))

         ok = ok .and. ga_create (mt_dcpl, params%ns_ao, params%ns_ao,
     $        "F in magnus inter "//amat(imat), 0, 0, g_zfock_ao(imat))

         ok = ok .and. ga_create (mt_dcpl, params%ns_ao, params%ns_ao,
     $        "P in magnus inter "//amat(imat), 0, 0, g_zdens_ao(imat))
         
      enddo
      if (.not. ok) call errquit (pname//"allocation failed",0,GA_ERR)


      ttdt = tt + params%dt

C
C     (1) Extrapolate to guess F'(t+1/2*dt) and F'(t+dt)
C
C     F'(t+1/2*dt) = 2*F'(t) - F'(t-1/2*dt)
C
      do imat = 1, nmats
         call ga_zero (g_zfock_p12dt(imat))
      
         call ga_add (z2, g_zfock(imat), zn1,
     $        g_zfock_n12dt(imat), g_zfock_p12dt(imat))
      
C
C     (2) Propagate P'(t) -> P'(t+dt)
C
         g_zdens(imat) = g_zdens_mo(imat)

         call ga_copy (g_zdens(imat), g_zdens_pdt(imat))

         call prop_magnus_ord2_step (params, params%dt,
     $        g_zfock_p12dt(imat), g_zdens_pdt(imat))
      enddo


      
C
C     Keep interpolating and propagating until result doesnt change.
C      
      converged = .false.
      iinter = 0
      num_same = 0

      do imat = 1, nmats
         diff(imat) = 999d0
      enddo
      
      do while (.not. converged)
         
         iinter = iinter + 1

C
C     (3) Build new Fock matrix (in AO basis) using newly created
C     P'(t+dt), and convert to MO basis.
C
C     Note that even though g_zfock_ao and g_zdens_ao are length 2, we
C     can send the array of GA handles to the closedshell/spinorbit fock
C     builder without issue (they will only use the first element).
C
         do imat = 1, nmats
            call ga_zero (g_zdens_ao(imat))
            call ga_zero (g_zfock_ao(imat))

            call canorg_trans (params, "P", "MO->AO",
     $           g_zdens_pdt(imat), g_zdens_ao(imat))
         enddo

         ok = tdfock (params, ttdt, g_zdens_ao, energies, g_zfock_ao)
         
         if (.not.ok) call errquit (pname//"tdfock failed", 0, 0)

         do imat = 1, nmats
            call canorg_trans (params, "F", "AO->MO",
     $           g_zfock_ao(imat), g_zfock_pdt(imat))
            

C
C     (4) Interpolate new F'(t+1/2*dt)
C
C     F'(t+1/2*dt) = 1/2*F'(t+dt) + 1/2*F'(t)
C
            call ga_zero (g_zfock_p12dt(imat))
         
            call ga_add (zhalf, g_zfock_pdt(imat), zhalf, g_zfock(imat),
     $           g_zfock_p12dt(imat))
            
         
C
C     (5) Propagate P'(t) -> P'(t+dt)
C
            call ga_copy (g_zdens(imat), g_zdens_pdt(imat))

            call prop_magnus_ord2_step (params, params%dt,
     $           g_zfock_p12dt(imat), g_zdens_pdt(imat))
         enddo

C
C     Compare old P'(t+dt) to one we just calculated--if same we are
C     done.  Note we compare alpha and beta densities separately, don't
C     compare if this is the first iteration.
C
         if (iinter .gt. 1) then
            do imat = 1, nmats 
               norm1(imat) = mat_norm (g_zdens_pdt(imat))
               norm2(imat) = mat_norm (g_zdens_pdt_old(imat))
               diff(imat) = abs (norm1(imat) - norm2(imat))
            enddo
         endif

         ok = .true.
         do imat = 1, nmats
            ok = ok .and. (diff(imat) .le. params%tol_interpol)
         enddo
         if (ok) num_same = num_same + 1
         
         do imat = 1, nmats
            if ( params%print_convergence .and. (iinter .gt. 1) ) then
               write (outstring, "(a,i0,a,i0,a)")
     $              "Magnus interpol (matrix ", imat, " of ", nmats, ")"

               call rt_tddft_print_val (params, tt,
     $              diff(imat), trim (outstring))
            endif
         enddo
         
         if (num_same .ge. params%terms_interpol) converged = .true.

C         
C     Skip interpolation.  Saves a Fock build but might be dangerous...
         if (lskip_interpol) then
            converged = .true.
            call rt_tddft_print_notice ("Skipped Magnus interpolation")
         endif

         
         if (iinter .gt. ninter_max)
     $        call errquit (pname//"failed to converge", 0, 0)


C
C     Store P'(t+dt) for comparison against next iteration.
C
         do imat = 1, nmats
            call ga_copy (g_zdens_pdt(imat), g_zdens_pdt_old(imat))
         enddo
         
      enddo


C
C     We're now effectively at t + dt: update stored Fock matrices.
C  
      do imat = 1, nmats
         
         call ga_copy (g_zfock_p12dt(imat), g_zfock_n12dt(imat))
         call ga_copy (g_zfock_pdt(imat), g_zfock(imat))
         call ga_copy (g_zdens_pdt(imat), g_zdens_mo(imat))
      

C
C     Clean up.
C
         ok = .true.
         ok = ok .and. ga_destroy (g_zdens_pdt(imat))
         ok = ok .and. ga_destroy (g_zdens_pdt_old(imat))
         ok = ok .and. ga_destroy (g_zfock_p12dt(imat))
         ok = ok .and. ga_destroy (g_zfock_pdt(imat))
         ok = ok .and. ga_destroy (g_zfock_ao(imat))
         ok = ok .and. ga_destroy (g_zdens_ao(imat))
      enddo
      
      if (.not. ok) call errquit (pname//"destroys failed",0,GA_ERR)

      
      if (params%print_convergence)
     $     call rt_tddft_print_ival (params, tt,
     $     iinter, "Magnus total interpolations")


C
C     If profiling is enabled print convergence data.
C      
      if (params%prof) then 
         write (outstring,"(a,i0,a)")
     $        "Magnus propagation converged after ", iinter,
     $        " interpolations;"
         call prof_end (elapsed, trim(outstring))
      endif

      end subroutine


         
c $Id$
