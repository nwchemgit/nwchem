C> \ingroup wfn1
C> @{
C>
C> \brief Evaluate a density functional with the wfn1 density
C>
C> The density functional is evaluated by calling xc_getv. This
C> routine expects a number of arguments that we need to provide
C> somehow. Here we generate the right arguments and make the call.
C> The result feeding back is purely the exchange-correlation energy.
C>
      subroutine wfn1_xc(rtdb,nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,exc,
     &                   dd,docc)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer rtdb  !< [Input] The RTDB handle
      integer nbf   !< [Input] The number of basis functions
      integer nea   !< [Input] The number of \f$\alpha\f$ spin electrons
      integer neb   !< [Input] The number of \f$\beta\f$ spin electrons
c
      double precision r_noa(nbf,nbf) !< [Input] ] The
      !< \f$\alpha\f$-spin natural orbitals. This is a 2D array with the
      !< dimensions:
      !< - 1. The basis functions
      !< - 2. The natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] ] The
      !< \f$\beta\f$-spin natural orbitals. This is a 2D array with the
      !< same dimensions as `r_noa`.
      double precision r_ofa(nbf,nbf) !< [Input] The
      !< \f$\alpha\f$-spin occupation functions. This is a 2D array with
      !< the dimension:
      !< - 1. The natural orbitals
      !< - 2. The occupation functions
      double precision r_ofb(nbf,nbf) !< [Input] The
      !< \f$\beta\f$-spin occupation functions. This is a 2D array with
      !< the same dimensions as `r_ofa`.
      double precision exc !< [Output] The exchange-correlation energy
c
      double precision docc(nbf) !< [Scratch] The occupation numbers
      double precision dd(nbf,nbf) !< [Scratch] The density matrix
c
      integer g_dens(2) !< The density matrix handles
      integer g_vxc(2)  !< The Fock matrix handles
c
      double precision dume, dum1, dum2, rho_n
      double precision e_vxc
      integer ia, ib, ii, ir
c
      logical  xc_gotxc
      external xc_gotxc
c
      exc = 0.0d0
      if (.not.xc_gotxc()) return
c
      if (.not.ga_create(MT_DBL,nbf,nbf,"dens a",-1,-1,g_dens(1)))
     &  call errquit("wfn1_xc: allocate GA dens a failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"dens b",-1,-1,g_dens(2)))
     &  call errquit("wfn1_xc: allocate GA dens b failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"Vxc a",-1,-1,g_vxc(1)))
     &  call errquit("wfn1_xc: allocate GA Vxc a failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"Vxc b",-1,-1,g_vxc(2)))
     &  call errquit("wfn1_xc: allocate GA Vxc b failed",0,GA_ERR)
c
c     Compute alpha density matrix
c
      do ii = 1, nbf
        docc(ii) = 0.0d0
      enddo
      do ir = 1, nea
        do ii = 1, nbf
          docc(ii) = docc(ii) + r_ofa(ii,ir)*r_ofa(ii,ir)
        enddo
      enddo
      call dfill(nbf*nbf,0.0d0,dd,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            dd(ia,ib) = dd(ia,ib)
     +                + r_noa(ia,ii)*r_noa(ib,ii)*docc(ii)
          enddo
        enddo
      enddo
      call ga_put(g_dens(1),1,nbf,1,nbf,dd,nbf)
c
c     Compute beta density matrix
c
      do ii = 1, nbf
        docc(ii) = 0.0d0
      enddo
      do ir = 1, neb
        do ii = 1, nbf
          docc(ii) = docc(ii) + r_ofb(ii,ir)*r_ofb(ii,ir)
        enddo
      enddo
      call dfill(nbf*nbf,0.0d0,dd,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            dd(ia,ib) = dd(ia,ib)
     +                + r_nob(ia,ii)*r_nob(ib,ii)*docc(ii)
          enddo
        enddo
      enddo
      call ga_put(g_dens(2),1,nbf,1,nbf,dd,nbf)
c
      call ga_zero(g_vxc(1))
      call ga_zero(g_vxc(2))
c
c     Do the deed
c
      call xc_getv(rtdb,exc,dume,1,0,0,g_dens,g_vxc,.false.,0,0,0,
     &             .false.,rho_n,dum1,dum2,0)
      if (abs(rho_n-dble(nea+neb)).gt.1.0d-4) then
        write(*,*)'No electrons problem: ',nea+neb,rho_n
      endif
      e_vxc = ga_ddot(g_dens(1),g_vxc(1))
     &      + ga_ddot(g_dens(2),g_vxc(2))
c
c     Clean up memory
c
      if (.not.ga_destroy(g_vxc(2)))
     &  call errquit("wfn1_xc: could not deallocate GA g_Vxc b",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_vxc(1)))
     &  call errquit("wfn1_xc: could not deallocate GA g_Vxc a",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_dens(2)))
     &  call errquit("wfn1_xc: could not deallocate GA g_dens b",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_dens(1)))
     &  call errquit("wfn1_xc: could not deallocate GA g_dens a",
     &               0,GA_ERR)
c
      end
C>
C> \brief Evaluate a density functional with the wfn1 density
C>
C> The density functional is evaluated by calling xc_getv. This
C> routine expects a number of arguments that we need to provide
C> somehow. Here we generate the right arguments and make the call.
C> The result feeding back is purely the exchange-correlation energy,
C> as well as the corresponding \f$\alpha\f$ and \f$\beta\f$-Fock 
C> matrices.
C>
      subroutine wfn1_xc_f(rtdb,nbf,nea,neb,densa,densb,exc,
     &                     focka,fockb,tfock)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer rtdb  !< [Input] The RTDB handle
      integer nbf   !< [Input] The number of basis functions
      integer nea   !< [Input] The number of alpha electrons
      integer neb   !< [Input] The number of beta electrons
c
      double precision densa(nbf,nbf) !< [Input] ] The
      !< \f$\alpha\f$-spin density matrix
      double precision densb(nbf,nbf) !< [Input] ] The
      !< \f$\beta\f$-spin density matrix
      double precision focka(nbf,nbf) !< [In/Output] The
      !< \f$\alpha\f$-spin Fock matrix
      double precision fockb(nbf,nbf) !< [In/Output] The
      !< \f$\beta\f$-spin Fock matrix
      double precision exc !< [Output] The exchange-correlation energy
      double precision tfock(nbf,nbf) !< [Scratch] Temporary matrix
c
      integer g_dens(2) !< The density matrix handles
      integer g_vxc(2)  !< The Fock matrix handles
c
      double precision dume, dum1, dum2, rho_n
      double precision e_vxc
      integer ia, ib, ii, ir
c
      logical  xc_gotxc
      external xc_gotxc
c
      exc = 0.0d0
      if (.not.xc_gotxc()) return
c
      if (.not.ga_create(MT_DBL,nbf,nbf,"dens a",-1,-1,g_dens(1)))
     &  call errquit("wfn1_xc: allocate GA dens a failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"dens b",-1,-1,g_dens(2)))
     &  call errquit("wfn1_xc: allocate GA dens b failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"Vxc a",-1,-1,g_vxc(1)))
     &  call errquit("wfn1_xc: allocate GA Vxc a failed",0,GA_ERR)
      if (.not.ga_create(MT_DBL,nbf,nbf,"Vxc b",-1,-1,g_vxc(2)))
     &  call errquit("wfn1_xc: allocate GA Vxc b failed",0,GA_ERR)
c
      call ga_put(g_dens(1),1,nbf,1,nbf,densa,nbf)
      call ga_put(g_dens(2),1,nbf,1,nbf,densb,nbf)
c
      call ga_zero(g_vxc(1))
      call ga_zero(g_vxc(2))
c
c     Do the deed
c
      call xc_getv(rtdb,exc,dume,1,0,0,g_dens,g_vxc,.false.,0,0,0,
     &             .false.,rho_n,dum1,dum2,0)
      if (abs(rho_n-dble(nea+neb)).gt.1.0d-4) then
        write(*,*)'No electrons problem: ',nea+neb,rho_n
      endif
      e_vxc = ga_ddot(g_dens(1),g_vxc(1))
     &      + ga_ddot(g_dens(2),g_vxc(2))
      call ga_get(g_vxc(1),1,nbf,1,nbf,tfock,nbf)
      call daxpy(nbf*nbf,1.0d0,tfock,1,focka,1)
      call ga_get(g_vxc(2),1,nbf,1,nbf,tfock,nbf)
      call daxpy(nbf*nbf,1.0d0,tfock,1,fockb,1)
c
c     Clean up memory
c
      if (.not.ga_destroy(g_vxc(2)))
     &  call errquit("wfn1_xc: could not deallocate GA g_Vxc b",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_vxc(1)))
     &  call errquit("wfn1_xc: could not deallocate GA g_Vxc a",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_dens(2)))
     &  call errquit("wfn1_xc: could not deallocate GA g_dens b",
     &               0,GA_ERR)
      if (.not.ga_destroy(g_dens(1)))
     &  call errquit("wfn1_xc: could not deallocate GA g_dens a",
     &               0,GA_ERR)
c
      end
C>
C> @}
