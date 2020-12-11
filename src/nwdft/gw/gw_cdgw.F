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
      integer g_epsilon(pars%ngl)
      integer igl,ngl

      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4

      double precision efermi
      double precision Iqp, Rqp
c
c     -----------------------------------------------------------------
c
      ngl = pars%ngl
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
      if (main) write(luout,9000)
 9000 format(10x,'Computing Vxc')

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

      call gw_sigmax(pars%g_erioo,pars%g_eriov,dbl_mb(k_sigmax),nri,nmo,
     $               nocc,nvir,pars%me)

c
c     Allocations
c
      if (.not.ga_create(mt_dbl,nri,nri,'Epsilon',-1,-1,g_epsilon(1)))
     $  call errquit(pname//'can''t create global array',0,GA_ERR)

      do igl=2,ngl
        if (.not.ga_duplicate(g_epsilon(1),g_epsilon(igl),'Eps(iw)'))
     $    call errquit(pname//'can''t create global array',0,GA_ERR)   
      enddo

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
      eviter = 0
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
c     Obtain Screened Coulomb matrices on imaginary axis
c
      if (main) write(luout,9030)
 9030 format(10x,'Computing W(iw) on the imaginary grid')
c
      call gw_buildiw(pars%g_eriov,dbl_mb(k_eia),dbl_mb(pars%k_glx),
     $   g_epsilon,npoles,nri,pars%ngl)
c
  200 continue ! evGW0
c
      if (pars%evgw0.and.main) then
        write(luout,9040) eviter
      endif
 9040 format(/,2x,"*** G",I1,"W0")     

      do iqp=1,nstates
        do i=1,pars%qpiter
          call gw_buildin(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    g_epsilon,dbl_mb(k_eia),Iqp,lo1+iqp,npoles,
     $                    nocc,nvir,nri,pars%ngl)     
c          call gw_buildrn
        enddo
      enddo
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildiw(ov,eia,x,g_epsilon,npoles,nri,ngl)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer,intent(in) :: npoles,nri,ngl
      integer,intent(in) :: ov,g_epsilon(ngl)
      double precision,intent(in) :: eia(npoles),x(ngl)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      integer g_a,g_b,g_c,g_d,g_e
      double precision factor(npoles),buf(npoles)


      integer ilo,ihi,me,k_ov,ld,mynpoles,igl,ipole,info
      double precision w2

      integer,external :: ga_nodeid
      logical,external :: nga_create,ga_destroy,ga_create,ga_duplicate

      me = ga_nodeid()

      if (.not.nga_create(mt_dbl,1,npoles,'factor',-1,g_a))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,nri,'one',-1,g_d))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,nri,'-one',-1,g_e))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.ga_duplicate(ov,g_c,'temp'))
     $  call errquit(pname//'can''t duplicate array',0,GA_ERR) 
      if (.not.ga_duplicate(g_epsilon(1),g_b,'temp'))
     $  call errquit(pname//'can''t duplicate array',0,GA_ERR) 

      call nga_distribution(g_a,me,ilo,ihi)
      call ga_access(g_c,1,nri,ilo,ihi,k_ov,ld)
      call ga_zero(g_d)
      call ga_zero(g_e)
      call ga_add_constant(g_d,1d0)
      call ga_add_constant(g_e,-1d0)

      mynpoles = ihi-ilo+1

      do igl=1,ngl
        w2 = x(igl)**2
        do ipole=ilo,ihi
          factor(ipole) = eia(ipole)/(eia(ipole)**2 + w2)
        enddo
        call nga_put(g_a,ilo,ihi,factor(ilo),1)
        call ga_sync()
        call ga_copy(ov,g_c)
        call ga_scale_cols(g_c,g_a)
        call ga_dgemm('n','t',nri,nri,npoles,4d0,g_c,ov,0d0,g_b)
        call ga_add_diagonal(g_b,g_d)
        call ga_chol_inv_seq(g_b,g_epsilon(igl),info)
        if (info.ne.0) call errquit(pname//'Cholesky failed',0,0)
      enddo

      if ((.not.ga_destroy(g_a)).or.(.not.ga_destroy(g_c)).or.
     $    (.not.ga_destroy(g_d)).or.(.not.ga_destroy(g_e)))
     $  call errquit(pname//'can''t destroy global arrays',0,GA_ERR) 
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildpi
      end subroutine

c     *****************************************************************

      subroutine gw_buildin(oo,ov,vv,eps,eig,Iqp,imo,npoles,nocc,nvir,
     $                      nri,ngl)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer,intent(in) :: npoles,nocc,nvir,nri,ngl
      integer,intent(in) :: oo,ov,vv,eps(ngl)
      integer,intent(in) :: imo
      double precision, intent(in) :: Iqp
      double precision, intent(in), dimension(nocc+nvir) :: eig

      integer ilo,ihi,nmo
      integer g_a
      double precision vmo(nocc+nvir),wmo(nocc+nvir),factor(nocc+nvir)

      nmo = nocc + nvir

      if (.not.ga_create(mt_dbl,nri,nmo,'temp',nri,-1,g_a))
     $  call errquit('gw_buildin: cannot create array',0,GA_ERR)
      if (.not.ga_duplicate(g_a,g_c,'temp'))
     $  call errquit('',0,GA_ERR)


      if (imo.le.nocc) then
        ilo = (imo-1)*nocc + 1
        ihi = imo*nocc
        call ga_copy_patch('n',oo,1,nri,ilo,ihi,g_a,1,nri,1,nocc)
        ilo = (imo-1)*(nmo-nocc) + 1
        ihi = imo*(nmo-nocc)
        call ga_copy_patch('n',ov,1,nri,ilo,ihi,g_a,1,nri,nocc+1,nmo)
      else
        do iocc=1,nocc
          ilo = (iocc-1)*(nmo-nocc) + imo - nocc
          call ga_copy_patch('n',ov,1,nri,ilo,ilo,g_a,1,nri,iocc,iocc)
        enddo
        ilo = (imo-nocc-1)*(nmo-nocc) + 1
        ihi = (imo-nocc)*(nmo-nocc)
        call ga_copy_patch('n',vv,1,nri,ilo,ihi,g_a,1,nri,nocc+1,nmo)
      endif

      call ga_distribution(g_a,ilo,ihi,jlo,jhi)
      call ga_access(g_a,1,nri,jlo,jhi,k_eri,ld)
      call ga_access(g_c,1,nri,jlo,jhi,k_tmp,ld)

      mynmo = jhi - jlo + 1
      do jmo=1,mynmo
        idx = k_eri + (jmo-1)*nri
        vmo(jmo) =  -ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
        factor(jmo) = omega - eig(jmo+jlo-1)
      enddo

      do igl=1,ngl
        call ga_dgemm('n','n',nri,nmo,nri,1d0,eps(igl),g_a,0d0,g_c)
        w2 = x(igl)**2
        do jmo=1,mynmo
          idx = k_tmp + (jmo-1)*nri
          wmo(jmo) = ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1) + vmo(jmo)
          factor2(jmo) =  w(igl)*factor(jmo)/(factor(jmo)**2 + w2)
        enddo
        Iqp = Iqp - ddot(mynmo,factor2,1,wmo,1)
      enddo
      call ga_dgop('+')
      

      end subroutine

c     *****************************************************************

      subroutine gw_buildrn
      end subroutine
