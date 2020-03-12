C
C     prop_rk4.F
C
C     Propagate MO density matrix/matricies using RK4.  Works for
C     non-Hermitian Fock matrix:
C
C     dP'/dt = -i (F'P' - P' F'^+)  Note the hermitian conjugate in last term.
C
C
      subroutine prop_rk4 (params, tt, tdfock, nmats, g_zdens_mo,
     $     energies)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  ! struct containing parameters
      double precision, intent(in)     :: tt
      logical, external                :: tdfock  !external function for computing F(t)
      integer, intent(in)              :: nmats     
      integer, intent(inout)           :: g_zdens_mo(nmats)
      type(rt_energies_t), intent(inout) :: energies


C     == Parameters ==
      character(*), parameter     :: pname       = "prop_rk4: "


      integer  g_zdens_ao (nmats), g_zfock_ao(nmats), g_zfock_mo(nmats)
      integer  g_ztmp1_mo(nmats), g_ztmp2_mo(nmats)
      double complex  zidt
      integer me
      integer g_zdens_1_mo(nmats), g_zdens_2_mo(nmats)
      integer g_zdens_3_mo(nmats), g_zdens_4_mo(nmats)
      integer g_ztmp_k1_mo(nmats), g_ztmp_k2_mo(nmats)
      integer g_ztmp_k3_mo(nmats), g_ztmp_k4_mo(nmats)
      integer g_ztemp3_mo(nmats)
      double precision elapsed
      integer imat
      logical ok

      if (params%prof) call prof_start (elapsed)

      
      me = ga_nodeid()

      if ( (nmats .ne. 1).and.(nmats .ne. 2) )
     $     call errquit (pname//"Only works for 1 or 2 matrices",0,0)


C
C     Create tempy GAs
C
      do imat=1,nmats
         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $        "zdens_ao", 0, 0, g_zdens_ao(imat)))
     $        call errquit ("failed to create zdens_ao", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $        "zfock_ao", 0, 0, g_zfock_ao(imat)))
     $        call errquit ("failed to create zfock_ao", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zfock_mo", 0, 0, g_zfock_mo(imat)))
     $        call errquit ("failed to create zfock_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp1_mo", 0, 0, g_ztmp1_mo(imat)))
     $        call errquit ("failed to create ztmp1_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k1_mo", 0, 0, g_ztmp_k1_mo(imat)))
     $        call errquit ("failed to create ztmp_k1_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k2_mo", 0, 0, g_ztmp_k2_mo(imat)))
     $        call errquit ("failed to create ztmp_k2_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k3_mo", 0, 0, g_ztmp_k3_mo(imat)))
     $        call errquit ("failed to create ztmp_k3_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k4_mo", 0, 0, g_ztmp_k4_mo(imat)))
     $        call errquit ("failed to create ztmp_k4_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp2_mo", 0, 0, g_ztmp2_mo(imat)))
     $        call errquit ("failed to create ztmp2_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_1_mo", 0, 0, g_zdens_1_mo(imat)))
     $        call errquit ("failed to create zdens_1_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_2_mo", 0, 0, g_zdens_2_mo(imat)))
     $        call errquit ("failed to create zdens_2_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_3_mo", 0, 0, g_zdens_3_mo(imat)))
     $        call errquit ("failed to create zdens_3_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_4_mo", 0, 0, g_zdens_4_mo(imat)))
     $        call errquit ("failed to create zdens_4_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztemp3_mo", 0, 0, g_ztemp3_mo(imat)))
     $        call errquit ("failed to create ztemp3_mo", 0, GA_ERR)
         
         
         call ga_zero (g_zfock_ao(imat))
         call ga_zero (g_zdens_ao(imat))
         call ga_zero (g_zfock_mo(imat))
         call ga_zero (g_ztmp1_mo(imat))
         call ga_zero (g_ztmp2_mo(imat))
         call ga_zero (g_zdens_1_mo(imat))
         call ga_zero (g_zdens_2_mo(imat))
         call ga_zero (g_zdens_3_mo(imat))
         call ga_zero (g_zdens_4_mo(imat))
         call ga_zero (g_ztmp_k1_mo(imat))
         call ga_zero (g_ztmp_k2_mo(imat))
         call ga_zero (g_ztmp_k3_mo(imat))
         call ga_zero (g_ztmp_k4_mo(imat))
         
      enddo
      
C
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
      
C
C     Compute F(t) in AO basis and convert to MO basis
C

      do imat=1,nmats
         call canorg_trans (params, "P", "MO->AO",  g_zdens_mo(imat),
     $        g_zdens_ao(imat))
         
      enddo
C     build Fock matrix in AO basis:
      if (.not. tdfock (params, tt, g_zdens_ao, energies,
     $     g_zfock_ao)) 
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      do imat=1,nmats
         call canorg_trans (params, "F", "AO->MO", g_zfock_ao(imat),
     $        g_zfock_mo(imat))
      enddo


             
       zidt = dcmplx (0d0, -1d0*params%dt)
         

c k1
C
C     Compute (dP'/dt)*dt = -i*dt(F'P' - P'F'^+) 
C     Note: hermitian conjugate of second term Fock matrix as we do not assume F is Hermitian
C
      do imat=1, nmats

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'

!     use g_ztemp3_mo(imat) temporarily for F'^+
         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))
         
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k1_mo(imat)) ! ztmp1_mo now has F'P - P'F'^+
         
         call ga_scale(g_ztmp_k1_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt
         
         
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C     
C     P = P + dt * dPdt 
         call ga_zadd(z1, g_zdens_mo(imat), z1/2, g_ztmp_k1_mo(imat), 
     $        g_zdens_1_mo(imat)) ! P = Pmo + (k1/6d0) or dp1/dt
         

c
c k2
c

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_1_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
         

!     use g_ztemp3_mo(imat) temporarily for F'^+
         
 
         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         
         
        call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_1_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))

          
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k2_mo(imat)) ! ztmp1_mo now has F'P - P'F'

         call ga_scale(g_ztmp_k2_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']
         


C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C      P = P + dt * dPdt 
         call ga_zadd(z1, g_zdens_mo(imat), z1/2, g_ztmp_k2_mo(imat), 
     $        g_zdens_2_mo(imat)) ! P = Pmo + (k1/6d0) + (k2/3d0) or dp2/dt
         


c
c k3
c

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_2_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'


         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_2_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))

         
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k3_mo(imat)) ! ztmp1_mo now has F'P - P'F'
         
         call ga_scale(g_ztmp_k3_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']
         
         
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C     P = P + dt * dPdt 
         call ga_zadd(z1, g_zdens_mo(imat), z1, g_ztmp_k3_mo(imat), 
     $        g_zdens_3_mo(imat)) ! P = Pmo + (k1/6d0) + (k2/3d0) +  + (k4/3d0) or dp3/dt
         

c
c k4
c

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_3_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
         
         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         

         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_3_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))
         
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k4_mo(imat)) ! ztmp1_mo now has F'P - P'F'
         
         call ga_scale(g_ztmp_k4_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']

C
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C     Pmo = Pmo + (k1/6d0) + (k2/3d0) + (k3/3d0) + (k4/6d0)
C       
         call ga_zadd(z1, g_zdens_mo(imat), z1/6,  g_ztmp_k1_mo(imat),
     $        g_ztmp1_mo(imat))
         call ga_zadd(z1, g_ztmp1_mo(imat), z1/3,  g_ztmp_k2_mo(imat),
     $        g_ztmp1_mo(imat))
         call ga_zadd(z1, g_ztmp1_mo(imat), z1/3,  g_ztmp_k3_mo(imat),
     $        g_ztmp1_mo(imat))
         call ga_zadd(z1, g_ztmp1_mo(imat), z1/6,  g_ztmp_k4_mo(imat),
     $        g_zdens_mo(imat))
         

c$$$  Normalizing and damping not needed for real time propagation
      
c$$$C     Re-normalize to have correct number of electrons
c$$$C     
c$$$         if (imat .eq. 1) then
c$$$            call rt_tddft_normalize (params, g_zdens_mo(imat), 
c$$$     $           g_zfock_mo(imat), "alpha")
c$$$         elseif (imat .eq. 2) then
c$$$            call rt_tddft_normalize (params, g_zdens_mo(imat), 
c$$$     $           g_zfock_mo(imat), "beta")
c$$$         else
c$$$            call errquit("invalid imat", 0, 0)
c$$$         endif
c$$$            
c$$$            
c$$$c     
c$$$C     
c$$$C     Damping: 0.8 old + 0.2 new
c$$$C     
c$$$         call ga_zadd(z1*0.75d0, g_zdens_mo(imat), z1*0.25d0,
c$$$     $        g_zdens_3_mo(imat), g_zdens_mo(imat))

      enddo

C
C     Now P is at P(t+dt), we're done!
C

C
C     Destroy arrays
C     

      do imat=1,nmats
         if (.not. ga_destroy (g_zfock_ao(imat)))
     $        call errquit (pname//"failed destruction zfock_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zfock_mo(imat)))
     $        call errquit (pname//"failed destruction zfock_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp1_mo(imat)))
     $        call errquit (pname//"failed destruction ztmp_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp2_mo(imat)))
     $        call errquit (pname//"failed destruction ztmp_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_ao(imat)))
     $        call errquit (pname//"failed destruction zfock_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_1_mo(imat)))
     $        call errquit (pname//"failed destruction zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_2_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_3_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_4_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp_k1_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k1_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp_k2_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k2_mo", 0
     $        ,0)
         if (.not. ga_destroy (g_ztmp_k3_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k3_mo", 0
     $        ,0)
         if (.not. ga_destroy (g_ztmp_k4_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k4_mo", 0
     $        ,0)

         if (.not. ga_destroy (g_ztemp3_mo(imat)))
     $        call errquit (pname//"failed destruction  ztemp3_mo", 0
     $        ,0)


         
      enddo
      
      if (params%prof) call prof_end (elapsed, "RK4 propagation")

      end subroutine
      
