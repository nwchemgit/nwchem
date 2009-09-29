c
c     Build Fock Matrix: Core Spin-Free Parts
c
      subroutine build_fock_core_sf(ao_bas_han, nbf_ao, 
     &  g_dens, g_fock, oskel)
c
      implicit none 
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ao_bas_han       ! ao basis handle 
      integer nbf_ao           ! ao basis functions 
      integer g_dens(2)        ! up and down density matrices (nbf x nbf)
      integer g_fock           ! fock matrix (nbf x nbf)
      integer g_fockso(2)      ! real and imag parts of the fock matrix (2*nbf x 2*nbf)
      double precision Ecore   ! spin-free core energy
      logical oskel            ! symmetrize fock matrix
c
c     Spin-free Components
c
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'potential', oskel)
c
      return
      end
c
c     Build Fock Matrix: Core Spin-Orbit Parts
c
      subroutine build_fock_core_so(ao_bas_han, nbf_ao, 
     &  g_so, g_densso, g_fockso, g_tmp, Ecoreso, oskel)
c
      implicit none 
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ao_bas_han       ! ao basis handle 
      integer nbf_ao           ! ao basis functions 
      integer g_so(3)          ! spin-orbit matrices: x,y,z (nbf x nbf each)
      integer g_densso(2)      ! real and imag parts of the density matrix (2*nbf x 2*nbf)
      integer g_fockso(2)      ! real and imag parts of the fock matrix (2*nbf x 2*nbf)
      integer g_tmp            ! work array (nbf x nbf)
      double precision Ecoreso ! spin-orbit core energy
      logical oskel            ! symmetrize the fock matrix
c
c     Spin-orbit components
c
      call ga_zero(g_so(1))
      call ga_zero(g_so(2))
      call ga_zero(g_so(3))
      call int_1e_ga(ao_bas_han, ao_bas_han, g_so, 'so', oskel)
      call ga_scale(g_so(1),dble(0.5d0))
      call ga_scale(g_so(2),dble(0.5d0))
      call ga_scale(g_so(3),dble(0.5d0))
      call ga_fock_so(g_so, g_fockso, nbf_ao)      ! pack g_so into full complex g_fockso
c
c     Accumulate the spin-orbit energy
c
      Ecoreso = 0.d0
c
c     Z-component SO contribution
c     Re(Dz)=-Im(Daa)+Im(Dbb)
c     <Hz>=Re(Dz) dot Vz
c
      call ga_zero(g_tmp)
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'z')
      Ecoreso = ga_ddot(g_tmp, g_so(1))
c
c     Y-component SO contribution
c     Re(Dy)=Re(Dab)-Re(Dba)
c     <Hy>=Re(Dy) dot Vy
c
      call ga_zero(g_tmp)
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'y')
      Ecoreso = Ecoreso + ga_ddot(g_tmp, g_so(2))
c
c     X-component SO contribution
c     Re(Dx)=-Im(Dab)-Im(Dba)
c     <Hx>=Re(Dx) dot Vx
c
      call ga_zero(g_tmp)
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'x')
      Ecoreso = Ecoreso + ga_ddot(g_tmp, g_so(3))
c
      return
      end
c
c     Build Fock Matrix: 2e & xc Parts
c
      subroutine build_fock_2e_xc(rtdb, nbf_ao, nbf_mo, 
     &   Exc, Ecoul, nExc, 
     &   iVxc_opt, g_xcinv, g_dens, g_vxc, IOLGC, 
     &   g_wght, g_xyz, g_nq, wght_GA, rho_n, 
     &   rdens_atom, cetobfr, natoms, ipol,
     &   nbf_ao_mxnbf_ce, g_tmp, g_fockso, g_densso) 
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c
         integer rtdb
         integer nbf_ao
         integer nbf_mo
         double precision Exc(*)
         double precision Ecoul
         integer nExc
         integer iVxc_opt
         integer g_xcinv
         integer g_dens(*)
         integer g_vxc(4)
         logical IOLGC
         integer g_wght
         integer g_xyz
         integer g_nq
         logical wght_GA
         double precision rho_n
         integer ipol
         integer natoms
         double precision rdens_atom(ipol*natoms*natoms)
         integer cetobfr(2,natoms)
         integer nbf_ao_mxnbf_ce         
         integer g_tmp            ! work array (nbf x nbf)
         integer g_fockso(*)      ! real and imag parts of the fock matrix (2*nbf x 2*nbf)
         integer g_densso(*)      ! real and imag parts of the density matrix (2*nbf x 2*nbf)
c
c        Local variables
         integer nscr
         integer lscr
         integer iscr
         character*32 pname
         pname = "build_fock_2e_xc:"
c
c        Pre-compute reduced alpha and beta density matrices over atoms
         call dfill(ipol*natoms*natoms, 0.0d0, rdens_atom, 1)
         nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit(pname//'Cannot allocate scr',0, MA_ERR)
         call util_ga_mat_reduce(nbf_ao, natoms, cetobfr,
     &     g_dens, ipol, rdens_atom,
     &     'rms', dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
         if (.not.ma_pop_stack(lscr))
     &     call errquit(pname//'Cannot pop scr:',0, MA_ERR)
c
c        Compute the XC potential and energy.
         call ga_zero(g_tmp)
         g_vxc(1) = g_tmp
         call ga_zero(g_vxc(1))
         rho_n = 0.0d0
         call ga_zero(g_vxc(2))
c
         Ecoul = 0.d0
         Exc(1) = 0.d0
         Exc(2) = 0.d0
         call xc_getv(rtdb, Exc, Ecoul, nExc, 
     &    iVxc_opt, g_xcinv, g_dens,
     &    g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &    rdens_atom, cetobfr, natoms)
c
c        Add the exchange-correlation to the fock matrix
         call ga_sync()
         call ga_dadd_patch( 1.d0, g_fockso(1), 1, nbf_ao,
     &     1, nbf_ao,
     &     1.0d0, g_vxc(1),  1, nbf_ao, 1, nbf_ao,
     &     g_fockso(1), 1, nbf_ao, 1, nbf_ao)
         call ga_dadd_patch( 1.d0, g_fockso(1), 1+nbf_ao, nbf_mo,
     &     1+nbf_ao, nbf_mo,
     &     1.0d0, g_vxc(2),  1, nbf_ao, 1, nbf_ao,
     &     g_fockso(1), 1+nbf_ao, nbf_mo, 1+nbf_ao, nbf_mo)
c
c        Add the exact exchange contribution
         call xc_exso_mod(rtdb,Exc,Ecoul,nExc,g_densso,g_fockso)
c
c        call ga_print(g_fockso(1))
c        call ga_print(g_fockso(2))
c
         return
         end
c
c     Build Total Fock Matrix
c
      subroutine build_fock_total(rtdb, ao_bas_han, nbf_ao, nbf_mo,
     &    Exc, Ecoul, nExc, iVxc_opt, g_xcinv, g_dens, g_vxc, IOLGC, 
     &    g_wght, g_xyz, g_nq, wght_GA, rho_n, rdens_atom, cetobfr, 
     &    natoms, ipol, nbf_ao_mxnbf_ce, g_tmp, g_fockso, g_densso, 
     &    g_so, g_fock, oskel, Etnew) 
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer rtdb
      integer ao_bas_han
      integer nbf_ao
      integer nbf_mo
      double precision Exc(2)
      double precision Ecoul
      integer nExc
      integer iVxc_opt
      integer g_xcinv
      integer g_dens(2)
      integer g_vxc(4)
      logical IOLGC
      integer g_wght
      integer g_xyz
      integer g_nq
      logical wght_GA
      integer natoms
      integer ipol
      double precision rho_n
      double precision rdens_atom(ipol*natoms*natoms)
      integer cetobfr(2,natoms)
      integer nbf_ao_mxnbf_ce         
      integer g_tmp            ! work array (nbf x nbf)
      integer g_fockso(2)      ! real and imag parts of the fock matrix (2*nbf x 2*nbf)
      integer g_densso(2)      ! real and imag parts of the density matrix (2*nbf x 2*nbf)
      integer g_so(3)          ! spin-orbit matrices: x,y,z (nbf x nbf each)
      integer g_fock
      logical oskel
      double precision Etnew
c
      double precision Ecore
      double precision Ecoreso
c
c     Local variables
      integer nscr
      integer lscr
      integer iscr

      Ecore = 0.d0
      Ecoreso = 0.d0
c
c     Build the spin-free core parts of the Fock matrix
      call build_fock_core_sf(ao_bas_han, nbf_ao,
     &  g_dens, g_fock, oskel)
c
c     Pack g_fock into real part g_fockso(1)
      call ga_fock_sf(g_fock, g_fockso(1), nbf_ao)
c
c     Calculate the core contribution energy
      call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))  ! combine density matrices
      Ecore = ga_ddot(g_dens(1), g_fock)  ! core parts depend on the full density
      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1)) ! restore density matrices
c
c     Build the spin-orbit core parts of the Fock matrix
      call build_fock_core_so(ao_bas_han, nbf_ao,
     &  g_so, g_densso, g_fockso, g_tmp, Ecoreso, oskel)
c
c     Full core energy
      Ecore = Ecore + Ecoreso
c
c     Build two-electron and exchange-correlation parts of the Fock matrix
      Ecoul = 0.d0
      Exc(1) = 0.d0
      Exc(2) = 0.d0
      call build_fock_2e_xc(rtdb, nbf_ao, nbf_mo, 
     &     Exc, Ecoul, nExc,
     &     iVxc_opt, g_xcinv, g_dens, g_vxc, IOLGC,
     &     g_wght, g_xyz, g_nq, wght_GA, rho_n,
     &     rdens_atom, cetobfr, natoms, ipol,
     &     nbf_ao_mxnbf_ce, g_tmp, g_fockso, g_densso)
c
      !write(6,*) "2epart:",Exc(1),Ecoul
c
      Etnew = 0.d0
      if (nExc.eq.1)then
        Etnew = Ecore + Ecoul + Exc(1)
      else
        Etnew = Ecore + Ecoul + Exc(1) + Exc(2)
      endif
c
      return
      end
