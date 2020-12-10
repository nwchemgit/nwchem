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

      integer nri, nmo, nocc, nvir, npoles, nstates
      integer nri_me, istart, iend, lo1, ld
      logical main, oactive

      integer k_eriij,k_eriia,k_sigma,k_sigmax,k_eold,k_enew,k_sign
      integer l_eriij,l_eriia,l_sigma,l_sigmax,l_eold,l_enew,l_sign
      integer k_omega,k_eia
      integer l_omega,l_eia

      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4

      double precision efermi
c
c     -----------------------------------------------------------------
c
      nri = nbf_cd
      nmo = pars%nmo(1)
      nocc = pars%nocc(1)
      nvir = pars%nvir(1)
      npoles = nocc*nvir
      nri_me = pars%nbf_cd_me
      nstates = pars%nstates
      main = pars%me .eq. 0
      oactive = pars%nbf_cd_me .gt. 0
      lo1 = nocc - pars%noqp
      istart = pars%istart
      iend = istart + nri_me - 1
c
c     Obtain Vxc
c
      if (main) write(luout) 9000
 9000 format(10x,'Computinv Vxc')

      if (.not.ma_push_get(mt_dbl,2*nmo,'Sigma',l_sigma,k_sigma))
     $  call errquit(pname//'cannot allocate Sigma',0,MA_ERR)

      call gw_vxc(dbl_mb(k_sigma),nmo,pars)
c
c     Obtain Sigma_x
c
      if (main) write(luout,9010)
 9010 format(10x,'Computing Sigma_x')

      if (.not.ma_push_get(mt_dbl,nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)

      call ga_access(pars%g_moeris,1,nmo*nocc,istart,iend,k_eriij,ld)
      call gw_sigmax(dbl_mb(k_eriij),dbl_mb(k_sigmax),nmo,nri_me,nocc)
c
c     Allocations
c
      if (.not.ma_push_get(mt_dbl,nri_me*npoles,'eri_ia',l_eriia,
     $  k_eriia))
     $  call errquit(pname//'cannot allocate eri_ia',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals old',l_eold,k_eold))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals new',l_enew,k_enew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'signs',l_sign,k_sign))
     $  call errquit(pname//'failed to allocate sign',0,MA_ERR)
       if (.not.ma_push_get(mt_dbl,npoles,'eval_ia',l_eia,k_eia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
c     Get fermi energy
c
      efermi = (dbl_mb(pars%k_evals+nocc-1) +
     $         dbl_mb(pars%k_evals+nocc))/2.0d0
c
c     Shift eigenvalues
c
      do i=0,nmo-1
        dbl_mb(k_eold+i) = dbl_mb(pars%k_evals+i) - efermi
      enddo
c
  100 continue ! evGW calculations
c
      if (pars%evgw.and.main) then
        write(luout,9020) eviter,eviter
      endif
 9020 format(/,2x,"*** G",I1,"W",I1)     
c
c     Obtain eigenvalue difference
c
      ipole = 0
      do i=0,nocc-1
        do j=nocc,nmo-1
          dbl_mb(k_eia+ipole) = dbl_mb(k_eold+j) - dbl_mb(k_eold+i)
          ipole = ipole + 1
        enddo
      enddo
c
c     Copy occupied-virtual into contiguous block
c
      lo2 = k_eriia
      do ibf=istart,iend
        do i=0,nocc-1
          ilo = nmo*i + 1
          ihi = ilo + nmo - 1
          call ga_access(pars%g_moeris,ilo,ihi,ibf,ibf,lo4,ld)
          call dcopy(nvir,dbl_mb(lo4+nocc),1,dbl_mb(lo2),1)
          lo2 = lo2 + nvir
        enddo
      enddo
c
c     Obtain Screened Coulomb matrices on imaginary axis
c
      if (main) write(luout,9030)
 9030 format(10x,'Computing W(iw) on the imaginary grid')
c
      call gw_buldiw(dbl_mb(k_eriia),dbl_mb(k_eia),dbl_mb(pars%k_glx),
     $   dbl_mb(pars%k_glw),npoles,nri_me,pars%ngl)
c
  200 continue ! evGW0
c
      if (pars%evgw0.and.main) then
        write(luout,9040) eviter
      endif
 9040 format(/,2x,"*** G",I1,"W0")     

      do iqp=1,nstates
        do i=1,pars%qpiter
c          call gw_buildin
c          call gw_buildrn
        enddo
      enddo
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildiw(eriia,eia,x2,w,npoles,nocc,nri,ngl)
      implicit none
      integer,intent(in) :: npoles,nocc,nri,ngl
      double precision,intent(in) :: eriia(npoles,nri),eia(npoles)
      double precision,intent(in) :: x2(ngl),w(ngl)

      logical ok
      integer g_a,g_w(ngl)
      double precision factor(npoles),buf(npoles)

c      call ga_create(mt_dbl,npoles,nri_tot,'eriia',npoles,-1,g_a)
c      call ga_distribution(g_a,me,ilo,ihi,jlo,jhi)
c
c      do igl=1,ngl
c        do ipoles=1,npoles
c          factor(ipoles) = dsqrt(eia(ipoles)/(eia(ipoles)**2 + x2(igl)))
c        enddo
c        do iri=1,nri
c          buf(:) = eriia(:,iri)*factor(:)
c          call ga_put(g_a,1,npoles,jlo+iri-1,jlo+iri-1,buf,1)
c        enddo
c        call ga_dgemm('t','n',nri_tot,nri_tot,npoles,2d0,g_a,g_a,0d0,
c     $                 g_w(igl))
c      enddo
c
c
c          
c      do iri=1,nri
c        factor = eia(:)/(eia**2 + x2(igl))
c      enddo

c      do g=1,ngrid
c        factor = ieta - eai
c        factor = 2.0d0*sqrt(factor/(factor**2 + glw2(g)))
c
c      enddo
      end subroutine

c     *****************************************************************

      subroutine gw_buildpi
      end subroutine

c     *****************************************************************

      subroutine gw_buildin
      end subroutine

c     *****************************************************************

      subroutine gw_buildrn
      end subroutine
