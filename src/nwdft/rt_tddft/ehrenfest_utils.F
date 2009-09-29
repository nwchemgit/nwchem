c
c     Calculate D.E
c
      subroutine calc_ddotE(t,g_dipx,g_dipy,g_dipz,g_ddotE)
c
      implicit none
c
#include "global.fh"
c
      double precision t
      integer g_dipx,g_dipy,g_dipz,g_ddotE
c
      double precision eet, Et, field_w
      double precision pi,field_Emax
      double precision evec(3)
c
      pi = 3.1415926535897932384626433832795d0
c
      field_Emax = 0.5d0
      field_w = 0.13d0
c
c      field_Emax = 0.05d0
c      field_w = 0.4764d0
c
c      field_Emax = 0.07d0
c      field_w = 0.10d0
c
c      field_Emax = 0.03d0
c      field_w = 0.10d0
c
c      field_Emax = 0.05d0
c      field_w = 0.277d0
c
c     Specify the electric field (from Xiaosong Li, PCCP, 7, 233, 2005)
      eet = field_Emax * sin(field_w*t)
c
c     Electric field 1
      if (t .ge. 0.d0) then
        Et = eet
      end if
c
#if 0
c
c     Electric field 2
      if ((t.ge.0.d0).and.(t.le.(2.d0*pi/field_w))) then
         Et = (field_w*t)/(2d0*pi) * eet
      elseif ((t.ge.(2.d0*pi/field_w)).and.(t.le.(4.d0*pi/field_w)))then
         Et = eet
      elseif ((t.ge.(4.d0*pi/field_w)).and.(t.le.(6.d0*pi/field_w)))then
         Et = (3.d0 - field_w*t/(2.d0*pi)) * eet
      else
         Et = 0.d0
      end if
#endif
c     
c     Assemble D.E  : Just X component for now
      evec(1) = 1
      evec(2) = 0
      evec(3) = 0
      call ga_zero(g_ddotE)
      call ga_dadd(evec(1)*Et,g_dipx,1.d0,g_ddotE,g_ddotE)
      call ga_dadd(evec(2)*Et,g_dipy,1.d0,g_ddotE,g_ddotE)
      call ga_dadd(evec(3)*Et,g_dipz,1.d0,g_ddotE,g_ddotE)
c
c     call ga_print(g_ddotE)

      write(20,*) t,Et
c
      return 
      end
c
c     Calculate energy of the dipole
c
      subroutine energy_dipole(g_dens,g_ddotE,Edip)
c
      implicit none
c
#include "global.fh"
c
      integer g_ddotE,g_dens(2)
      double precision Edip
c
      Edip = 0.d0
      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))      ! combine density matrices
      Edip = ga_ddot(g_dens(1), g_ddotE)
      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1)) ! restore density matrices
c
c     call ga_print(g_ddotE)
c
      write(22,*) "Edip:",Edip
c
      return
      end
c
c     Calculate dP/dt = -i*(S^-1*F*P - P*F*S^-1)
c     Fock Matrix: Real = g_fockso(1), Imag = g_fockso(2)
c     Density Matrix: Real = g_densso(1), Imag = g_densso(2)
c
      subroutine calc_dpdt(nbf_mo,g_fockso,g_densso,g_sm1,g_work,
     &          g_dprdt,g_dpidt)
c
      implicit none
c
#include "global.fh"
#include "../so_dft/consts.fh"
c
      integer nbf_mo
      integer g_fockso(2) ! real and imag parts of Fock matrix
      integer g_densso(2) ! real and imag parts of the density matrix
      integer g_sm1       ! S^(-1)
      integer g_work(*)   ! work arrays
      integer g_dprdt     ! time derivative of the real part of the density matrix
      integer g_dpidt     ! time derivative of the imag part of the density matrix
c
      ! FIX ME: Need to check size of work array here
c 
      ! Real part
      ! dPr/dt = S^-1*(Fr*Pi + Fi*Pr) - (Pr*Fi + Pi*Fr)*S^(-1)
      ! S^-1*(Fr*Pi + Fi*Pr)
      call ga_zero(g_work(1))  ! initialise work arrays
      call ga_zero(g_work(2))  ! initialise work arrays
      call ga_zero(g_dprdt)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_fockso(1), g_densso(2), 0.d0, g_work(1))    ! Fr*Pi
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_fockso(2), g_densso(1), 0.d0, g_work(2))    ! Fi*Pr
      call ga_dadd(1.d0, g_work(1), 1.d0, g_work(2), g_work(1)) ! Fr*Pi + Fi*Pr
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_sm1, g_work(1), 0.d0, g_dprdt)              ! S^-1*(Fr*Pi + Fi*Pr)
      ! (Pr*Fi + Pi*Fr)*S^(-1)
      call ga_zero(g_work(1))  ! initialise work arrays
      call ga_zero(g_work(2))  ! initialise work arrays
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_densso(1), g_fockso(2), 0.d0, g_work(1))    ! Pr*Fi
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_densso(2), g_fockso(1), 0.d0, g_work(2))    ! Pi*Fr
      call ga_dadd(1.d0, g_work(1), 1.d0, g_work(2), g_work(1)) ! Pr*Fi + Pi*Fr
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_work(1), g_sm1, 0.d0, g_work(2))            ! (Pr*Fi + Pi*Fr)*S^(-1)
      !  S^-1*(Fr*Pi + Fi*Pr) - (Pr*Fi + Pi*Fr)*S^(-1)
      call ga_dadd(1.d0, g_dprdt, -1.d0, g_work(2), g_dprdt)
c
      ! Imaginary part
      ! dPi/dt = S^(-1)*(Fi*Pi - Fr*Pr) + (Pr*Fr - Pi*Fi)*S^(-1)
      ! S^(-1)*(Fi*Pi - Fr*Pr)
      call ga_zero(g_work(1))  ! initialise work arrays
      call ga_zero(g_work(2))  ! initialise work arrays
      call ga_zero(g_dpidt)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_fockso(2), g_densso(2), 0.d0, g_work(1))    ! Fi*Pi
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_fockso(1), g_densso(1), 0.d0, g_work(2))    ! Fr*Pr
      call ga_dadd(1.d0, g_work(1), -1.d0, g_work(2), g_work(1))! Fi*Pi - Fr*Pr
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_sm1, g_work(1), 0.d0, g_dpidt)              ! S^(-1)*(Fi*Pi - Fr*Pr)
      ! (Pr*Fr - Pi*Fi)*S^(-1)
      call ga_zero(g_work(1))  ! initialise work arrays
      call ga_zero(g_work(2))  ! initialise work arrays
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_densso(1), g_fockso(1), 0.d0, g_work(1))    ! Pr*Fr
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_densso(2), g_fockso(2), 0.d0, g_work(2))    ! Pi*Fi
      call ga_dadd(1.d0, g_work(1), -1.d0, g_work(2), g_work(1))! Pr*Fr - Pi*Fi
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, 1.d0,
     &            g_work(1), g_sm1, 0.d0, g_work(2))            ! (Pr*Fr - Pi*Fi)*S^(-1)
      ! S^(-1)*(Fi*Pi - Fr*Pr) + (Pr*Fr - Pi*Fi)*S^(-1)
      call ga_dadd(1.d0, g_dpidt, 1.d0, g_work(2), g_dpidt)
c
      return
      end
c
c     Propagate density matrix using a 4th order Runge Kutta
c
      subroutine prop_densmat_rk4(nbf_mo,g_sm1,g_fockso,g_densso,
     &    g_work,g_k1,g_k2,g_k3,g_k4,g_kwork,t,dt)
c
      implicit none
c
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer nbf_mo
      integer g_sm1  ! S^(-1)
      integer g_fockso(2)
      integer g_densso(2)
      double precision t
      double precision dt
c
      integer g_work(2) ! work arrays
      integer g_k1(2) ! work arrays
      integer g_k2(2) ! work arrays
      integer g_k3(2) ! work arrays
      integer g_k4(2) ! work arrays
      integer g_kwork(2) ! work arrays
c
      integer irk
      double precision xh
      double precision xh2
      double precision h2
      double precision h6
c
c     Perform 4th order Runge Kutta
c
      do irk = 1,4
         xh = t + dt
         xh2 = t + dt/2d0
         h2 = dt/2d0
         h6 = dt/6d0
         call ga_zero(g_work(1))
         call ga_zero(g_work(2))
         call ga_zero(g_kwork(1))
         call ga_zero(g_kwork(2))
         if (irk.eq.1) then
             call calc_dpdt(nbf_mo,g_fockso,g_densso,g_sm1,g_work,
     &          g_k1(1),g_k1(2))
         else if (irk.eq.2) then
             call ga_dadd(1.d0,g_densso(1),h2,g_k1(1),g_kwork(1))
             call ga_dadd(1.d0,g_densso(2),h2,g_k1(2),g_kwork(2))
             call calc_dpdt(nbf_mo,g_fockso,g_kwork,g_sm1,g_work,
     &          g_k2(1),g_k2(2))
         else if (irk.eq.3) then
             call ga_dadd(1.d0,g_densso(1),h2,g_k2(1),g_kwork(1))
             call ga_dadd(1.d0,g_densso(2),h2,g_k2(2),g_kwork(2))
             call calc_dpdt(nbf_mo,g_fockso,g_kwork,g_sm1,g_work,
     &          g_k3(1),g_k3(2))
         else  ! irk.eq.4
             call ga_dadd(1.d0,g_densso(1),dt,g_k3(1),g_kwork(1))
             call ga_dadd(1.d0,g_densso(2),dt,g_k3(2),g_kwork(2))
             call calc_dpdt(nbf_mo,g_fockso,g_kwork,g_sm1,g_work,
     &          g_k4(1),g_k4(2))
         end if
      end do
c         
c     Assemble the step
c
      ! k1+2k2
      call ga_add(1d0,g_k1(1),2d0,g_k2(1),g_kwork(1))
      call ga_add(1d0,g_k1(2),2d0,g_k2(2),g_kwork(2))
      ! (k1+2k2)+2k3
      call ga_add(1d0,g_kwork(1),2d0,g_k3(1),g_kwork(1))
      call ga_add(1d0,g_kwork(2),2d0,g_k3(2),g_kwork(2))
      ! (k1+2k2+2k3)+k4
      call ga_add(1d0,g_kwork(1),1d0,g_k4(1),g_kwork(1))
      call ga_add(1d0,g_kwork(2),1d0,g_k4(2),g_kwork(2))
      ! y_{i+1} = y_i + h/6 (k1 + 2k2 + 2k3 + k4)
      call ga_add(1d0,g_densso(1),h6,g_kwork(1),g_densso(1))
      call ga_add(1d0,g_densso(2),h6,g_kwork(2),g_densso(2))
c
      return
      end
c
c     Propagate density matrix using the Euler method
c
      subroutine prop_densmat_euler(nbf_mo,g_sm1,g_fockso,g_densso,
     &    g_work,g_dprdt,g_dpidt,t,dt)
c
      implicit none
c
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer nbf_mo
      integer g_sm1  ! S^(-1)
      integer g_fockso(2)
      integer g_densso(2)
      integer g_work(2)
      integer g_dprdt
      integer g_dpidt
      double precision t
      double precision dt
c
      integer irk
      double precision xh
      double precision xh2
      double precision h2
      double precision h6
c
c     Propagate the step using the Euler scheme
c
      call calc_dpdt(nbf_mo,g_fockso,g_densso,g_sm1,g_work,
     &          g_dprdt,g_dpidt)
      ! g_densso(1) = g_densso(1) + dt * g_dprdt
      call ga_dadd(1.d0,g_densso(1),dt,g_dprdt,g_densso(1))
      ! g_densso(2) = g_densso(2) + dt * g_dpidt
      call ga_dadd(1.d0,g_densso(2),dt,g_dpidt,g_densso(2))
c
      return
      end
c
c     Calculate the instantaneous dipole moment
c
      subroutine calc_dipole(rtdb, ao_bas_han, ipol, g_dens,dipole,t)
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb              ! runtime database handle
      integer ao_bas_han        ! basis set handle
      integer ipol              ! spins
      integer g_dens(2)         ! components of the density
      double precision t        ! time
c
c     Performs a multipole analysis of the wavefunction and store dipole
c     in the database
c
      logical oprintmultipole
c
      integer geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments) 
      double precision alphamoments (len_moments)
      double precision betamoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      double precision dipole(*)
      integer k, l, m, ltotal, ind
      integer i
      data center/3*0.0d0/
      character*32 pname
      pname = "calc_dipole:"
c
      if (.not. bas_cando_mpoles(ao_bas_han)) return
      if (.not. bas_geom(ao_bas_han, geom)) 
     &  call errquit(pname//'bas_geom failed',0, BASIS_ERR)
c
c     Initialisation
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, alphamoments, 1)
      call dfill(len_moments, 0.0d0, betamoments, 1)
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, alphamoments, 1)
      call dfill(len_moments, 0.0d0, betamoments, 1)
c
c     Calculate the multipole components of the density and nucleus
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(ao_bas_han, center, lmax, g_dens(1), 
     &     alphamoments, len_moments)
      if (ipol .gt. 1) then
         call multipole_density(ao_bas_han, center, lmax, g_dens(2), 
     &     betamoments, len_moments)
      else
         do k = 1, len_moments
            alphamoments(k) = 0.5d0 * alphamoments(k)
            betamoments(k)  = alphamoments(k)
         enddo
      endif
      do k = 1, len_moments
        totalmoments(k) = alphamoments(k) + betamoments(k) + 
     &        nuclmoments(k)
      enddo
c
c     Store the dipole moment
      if (.not. rtdb_put(rtdb,'ehrenfest:dipole', mt_dbl, 
     & 3,totalmoments(2))) 
     & call errquit(pname//'Failed to store dipole moments',0,RTDB_ERR)
c     
c     Retrieve dipoles      
      if (.not.rtdb_get(rtdb,'ehrenfest:dipole',mt_dbl,3,dipole))
     &  call errquit(pname//'Failed to get dipole moments',0,RTDB_ERR)
c
      return
      end
c
c     Perform a Mulliken Population Analysis
c
      subroutine calc_mullpop(rtdb,geom,ao_bas_han,g_s,g_dens)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "../so_dft/consts.fh"
c
      integer rtdb        ! runtime database handle
      integer geom        ! geometry handle
      integer ao_bas_han  ! basis set handle
      integer g_s         ! overlap
      integer g_dens(2)   ! density components
c
c     Consolidate the densities
c
c     Perform a Mulliken Population Analysis
c
c     Total density analysis
      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))  ! total density
      call mull_pop(geom,ao_bas_han,g_dens(1),g_s,'total')
      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1)) ! separate components
c
c     Spin density analysis
c      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1)) ! spin density
c      call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'spin')
c      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))  ! separate components
c
      return
      end
c
c     Calculate the instantaneous dipole moment
c
      subroutine calc_dipole2(rtdb, ao_bas_han, nbf_ao, ipol, g_dens, 
     &  g_dipx, g_dipy, g_dipz, g_temp, t)
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb              ! runtime database handle
      integer ao_bas_han        ! basis set handle
      integer nbf_ao
      integer ipol              ! spins
      integer g_dens(2)         ! components of the density
      integer g_dipx   
      integer g_dipy   
      integer g_dipz   
      integer g_temp
      double precision t        ! time
c
      double precision mux
      double precision muy
      double precision muz
      double precision ga_trace_diag
      external ga_trace_diag
c
      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))  ! total density
c
      call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &            g_dens(1), g_dipx, 0.d0, g_temp) 
      mux = -ga_trace_diag(g_temp)
      call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &            g_dens(1), g_dipy, 0.d0, g_temp)
      muy = -ga_trace_diag(g_temp)
      call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &            g_dens(1), g_dipz, 0.d0, g_temp)
      muz = -ga_trace_diag(g_temp)
c
      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1)) ! separate
c
      write(6,*) "Dipole:",t,mux,muy,muz
c
      return
      end
c
c     Calculate the instantaneous current
c     jvec = Sum_{mu,nu} A_(mu,nu) phi_mu del phi_nu
c     where A_(mu,nu) is the antisymmetric or imaginary part 
c     of the density matrix
c
      subroutine calc_current(ao_bas_han, nbf_ao, ipol, g_dens,
     &  g_vel, g_temp, j, t)
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer ao_bas_han  ! basis set handle
      integer nbf_ao      ! number of basis functions
      integer ipol        ! spins
      integer g_dens(2)   ! components of the imaginary part of the density
      integer g_vel       ! velocity integral
      integer g_temp      ! scratch space
      double precision j(*) ! current components
      double precision t  ! time
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer idir
c
      double precision ga_trace_diag
      external ga_trace_diag
c
c     Total imaginary part
      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))

c     Calculate the current components
      alo(1) = 1              ! range of A, first dimension
      ahi(1) = nbf_ao
      alo(2) = 1              ! range of A, second dimension
      ahi(2) = nbf_ao
      blo(1) = 1              ! range of B, first dimension, etc.
      bhi(1) = nbf_ao
      blo(2) = 1
      bhi(2) = nbf_ao
      clo(1) = 1
      chi(1) = nbf_ao
      clo(2) = 1
      chi(2) = nbf_ao

      do idir = 1,3
          blo(3) = idir         ! B is here the dipole matrix with
          bhi(3) = idir         ! dimension (nbf, nbf, 3)
          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dens(1),alo,ahi,
     &       g_vel,blo,bhi,
     &       g_temp,clo,chi)
          j(idir) = ga_trace_diag(g_temp)
      end do
c
c     Separate components
      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1))
c
      !write(6,*) "Current-x:",t,j(1)
      !write(6,*) "Current-y:",t,j(2)
      !write(6,*) "Current-z:",t,j(3)
c
      return
      end
