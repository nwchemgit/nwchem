      subroutine gw_cdgw(pars)
c
c     Contour Deformation GW
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"

      type(gw_params_t) :: pars
      character(len=9), parameter :: pname = 'gw_cdgw: '

      double precision gw_fermi
      external gw_fermi
c
c     -----------------------------------------------------------------
c
      nri = nbf_cd
      nmo = pars%nmo(1)
      nocc = pars%nocc(1)
      nvir = pars%nvir(1)
      npoles = nocc*nir
c
c     Get fermi energy
c
      efermi = gw_fermi(dbl_mb(pars%k_evals),nocc)
c
c     Shift eigenvalues
c
      do i=0,nmo-1
        dbl_mb(k_eold+i) = dbl_mb(pars%k_evals+i)
      enddo
c
c     Copy occupied-virtual block of ERIs into a contiguous block
c
      lo1 = pars%k_3ceri + npoles
      lo2 = k_eriia
      nblk = nvir*nri
      do i=0,nocc-1
        call dcopy(nblk,dbl_mb(lo1),1,dbl_mb(lo2),1)
        lo1 = lo1 + nmo*nri - i*nri
        lo2 = lo2 + nblk
      enddo
c
c     Obtain Sigma_x
c      
      write(luout,9000)
 9000 format(10x,'Computing Sigma_x')

      call gw_sigmax(dbl_mb(pars%k_3ceri),dbl_mb(k_sigmax),nmo,nri,nocc)
c
c     Obtain Vxc
c
      write(luout,9010)
 9010 format(10x,'Computing V_xc')

      call gw_vxc(dbl_mb(k_vxc),nmo,pars)
c
c     Obtain Screened Coulomb matrices on imaginary axis
c
c      call gw_buldiw
c
      do iqp=1,pars%noqp+pars%nvqp
        do i=1,pars%qpiter
c          call gw_buildin
c          call gw_buildrn
        enddo
      enddo
c
      end subroutine
