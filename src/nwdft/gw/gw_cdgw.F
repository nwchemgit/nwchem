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
      integer k_vmo, l_vmo
      integer g_epsilon(pars%ngl)
      integer igl,ngl,nqp

      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4

      double precision efermi
      double precision I_n, R_n, dI_n, dR_n, dI_np, dI_nm, dR_np, dR_nm
      double precision d2I_n, d2R_n, dummy
      double precision sigma, sigma_c, dsigma_c, d2sigma_c
      double precision eps_in, eps_out, z, func, linear, dfunc, d2func
      double precision sigma_x, sigma_old

      logical converged
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
      nqp = pars%noqp + pars%nvqp
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
      if (.not.ma_push_get(mt_dbl,nqp*nmo,'Vmo',l_vmo,k_vmo))
     $  call errquit(pname//'cannot allocate Vmo',0,MA_ERR)

      call gw_vmo(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $            dbl_mb(k_vmo),dbl_mb(k_sigmax+lo1),nqp,nmo,nocc,
     $            nvir,nri,lo1+1,pars%me,pars%nvqp)
c
c     Allocations
c
      if (.not.ga_create(mt_dbl,nmo,nqp,'Epsilon',-1,nqp,g_epsilon(1)))
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
      call gw_buildiw(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $    dbl_mb(k_eia),dbl_mb(pars%k_glx),g_epsilon,dbl_mb(k_vmo),
     $    npoles,nmo,nocc,nvir,nri,pars%ngl,nqp,lo1+1)
c
  200 continue ! evGW0
c
      if (pars%evgw0.and.main) then
        write(luout,9040) eviter
      endif
 9040 format(/,2x,"*** G",I1,"W0")     

      if (main) write(luout,9050)
 9050 format(10x,'Solving quasiparticle equations')

      if (main) write(luout,9060)
 9060 format(/,10x,'State',5x,'Energy (eV)',5x,'Z',8x,'Elin (eV)',
     $       /,10x,45('-'))

      do iqp=1,pars%noqp+pars%nvqp
        eps_out = dbl_mb(k_eold+lo1+iqp-1)
        sigma_x = dbl_mb(k_sigmax+lo1+iqp-1)
        sigma_old = dbl_mb(k_sigma+lo1+iqp-1)
        do i=1,pars%qpiter
          eps_in = eps_out
          call gw_buildin(g_epsilon,dbl_mb(k_eold),dbl_mb(pars%k_glx),
     $                    dbl_mb(pars%k_glw),eps_in,I_n,dI_n,iqp,
     $                    npoles,nocc,nvir,nri,pars%ngl,pars%me)
          call gw_buildin(g_epsilon,dbl_mb(k_eold),dbl_mb(pars%k_glx),
     $                    dbl_mb(pars%k_glw),eps_in+0.001d0,dI_np,dummy,
     $                    iqp,npoles,nocc,nvir,nri,pars%ngl,pars%me)
          call gw_buildin(g_epsilon,dbl_mb(k_eold),dbl_mb(pars%k_glx),
     $                    dbl_mb(pars%k_glw),eps_in-0.001d0,dI_nm,dummy,
     $                    iqp,npoles,nocc,nvir,nri,pars%ngl,pars%me)
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in,R_n,
     $                    lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in+0.001d0,
     $                    dR_np,lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in-0.001d0,
     $                    dR_nm,lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          sigma_c = R_n + I_n
          dsigma_c = (dR_np - dR_nm + dI_np - dI_nm)/0.002d0
          sigma = sigma_c + sigma_x

          z = 1.0d0/(1.0d0-dsigma_c)
            
          func = dbl_mb(k_eold+lo1+iqp-1) + sigma - sigma_old - eps_in
          eps_out = eps_in + 0.667*z*func

          if (i.eq.1) linear = eps_in + z*func

          if (main) converged = abs(eps_out-eps_in).lt.pars%thresh
          call ga_brdcst(1001,converged,ma_sizeof(mt_log,1,mt_byte),0)
          call ga_brdcst(1002,eps_out,ma_sizeof(mt_dbl,1,mt_byte),0)
          if (converged) exit
        enddo
        dbl_mb(k_enew+lo1+iqp-1) = eps_out
        dbl_mb(k_sigma+lo1+iqp-1) = sigma_c + sigma_x
        if (main) write(luout,9070) lo1+iqp,(eps_out+efermi)*ha2ev,z,
     $                              (linear+efermi)*ha2ev     
      enddo

      if (pars%evgw0.or.pars%evgw) then
        call ga_brdcst(1003,dbl_mb(k_enew),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call ga_brdcst(1004,dbl_mb(k_sigma),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call dcopy(nmo,dbl_mb(k_enew),1,dbl_mb(k_eold),1)
        eviter = eviter + 1      
        if (eviter.lt.pars%eviter) then
          if (pars%evgw0) goto 200
          if (pars%evgw)  goto 100
        endif  
      endif

 9070 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3)

      do igl=1,pars%ngl
        if (.not.ga_destroy(g_epsilon(igl)))
     $    call errquit(pname//'cannot destroy global arrays',0,GA_ERR)
      enddo
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildiw(oo,ov,vv,eia,x,g_iw,vmo,npoles,nmo,nocc,
     $                      nvir,nri,ngl,nqp,llmo)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer,intent(in) :: npoles,nmo,nocc,nvir,nri,ngl,nqp,llmo
      integer,intent(in) :: oo,ov,vv,g_iw(ngl)
      double precision,intent(in) :: eia(npoles),x(ngl)
      double precision,intent(in) :: vmo(nmo,nqp)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      integer g_a,g_b,g_c,g_d,g_e,g_f
      double precision factor(npoles),buf(npoles)


      integer ilo,ihi,me,k_ov,ld,mynpoles,igl,ipole,info,k_c
      integer pilo,pihi,k_pi,l_pi,idx,idy,mynri,offset,l_tmp,tmp
      integer oolo,oohi,ovlo,ovhi,vvlo,vvhi,jmo,kmo,imo
      double precision w2,wmo

      double precision,external :: ddot
      integer,external :: ga_nodeid
      logical,external :: nga_create,ga_destroy,ga_create,ga_duplicate

      me = ga_nodeid()

      if (.not.nga_create(mt_dbl,1,npoles,'factor',-1,g_a))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.ga_create(mt_dbl,nri,nri,'Pi',nri,-1,g_b))
     $  call errquit(pname//'can''t duplicate array',0,GA_ERR) 
      if (.not.ga_create(mt_dbl,nri,nri,'EPS',nri,-1,g_f))
     $  call errquit(pname//'can''t duplicate array',0,GA_ERR)
      if (.not.ma_push_get(mt_dbl,nri*nri,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 

      call ga_distribution(oo,me,ilo,ihi,oolo,oohi)
      call ga_distribution(ov,me,ilo,ihi,ovlo,ovhi)
      call ga_distribution(g_b,me,ilo,ihi,pilo,pihi)
      call nga_distribution(g_a,me,ilo,ihi)
      call ga_access(ov,1,nri,ilo,ihi,k_ov,ld)
      call ga_access(g_b,1,nri,pilo,pihi,idx,ld)

      mynri = pihi - pilo + 1
      mynpoles = ihi-ilo+1
      offset = (pilo-1)*nri

      do igl=1,ngl

        w2 = x(igl)
        do ipole=ilo,ihi
          factor(ipole) = dsqrt(eia(ipole)/(eia(ipole)**2 + w2))
        enddo
        call gw_buildpi(dbl_mb(k_ov),dbl_mb(k_pi),factor(ilo),
     $                  nri,mynpoles,pilo,pihi)
        call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')
        call dcopy(nri*mynri,dbl_mb(k_pi+offset),1,dbl_mb(idx),1)
        call ga_sync()
        call ga_chol_inv_seq(g_b,g_f,info)
        call ga_get(g_f,1,nri,1,nri,dbl_mb(k_pi),nri)

        do imo=llmo,llmo+nqp-1
          jmo = imo - llmo + 1
          if (imo.le.nocc) then
            do kmo=1,nocc
              idy = (imo-1)*nocc + kmo
              if (idy.lt.oolo) cycle
              if (idy.gt.oohi) exit
              call ga_get(oo,1,nri,idy,idy,dbl_mb(tmp),nri)
              call dtrmv('l','n','n',nri,dbl_mb(k_pi),nri,dbl_mb(tmp),1)
              wmo = ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1) + vmo(kmo,jmo)
              call ga_put(g_iw(igl),kmo,kmo,jmo,jmo,wmo,1)
            enddo
            do kmo=1,nvir
              idy = (imo-1)*nvir + kmo
              if (idy.lt.ovlo) cycle
              if (idy.gt.ovhi) exit
              call ga_get(ov,1,nri,idy,idy,dbl_mb(tmp),nri)
              call dtrmv('l','n','n',nri,dbl_mb(k_pi),nri,dbl_mb(tmp),1)
              wmo = ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1) + 
     $                   vmo(kmo+nocc,jmo)
              call ga_put(g_iw(igl),kmo+nocc,kmo+nocc,jmo,jmo,wmo,1)
            enddo
          else
            do kmo=1,nocc
              idy = (kmo-1)*nvir + imo - nocc
              if (idy.lt.ovlo) cycle
              if (idy.gt.ovhi) exit
              call ga_get(ov,1,nri,idy,idy,dbl_mb(tmp),nri)
            enddo
            do kmo=1,nvir
              idy = (kmo-1)*nvir + kmo
              if (idy.lt.vvlo) cycle
              if (idy.gt.vvhi) exit
              call ga_get(vv,1,nri,idy,idy,dbl_mb(tmp),nri)
            enddo
          endif
        enddo

      enddo

      if ((.not.ga_destroy(g_a)).or.
     $    (.not.ga_destroy(g_b)).or.(.not.ga_destroy(g_f)))
     $  call errquit(pname//'can''t destroy global arrays',0,GA_ERR)
      if (.not.ma_chop_stack(l_pi))
     $  call errquit(pname//'can''t pop stack',0,MA_ERR) 
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildin(eps,eig,x,w,omega,Iqp,dIqp,imo,npoles,
     $                      nocc,nvir,nri,ngl,me)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: npoles,nocc,nvir,nri,ngl
      integer,intent(in) :: eps(ngl)
      integer,intent(in) :: imo, me
      double precision, intent(in), dimension(nocc+nvir) :: eig
      double precision, intent(in), dimension(ngl) :: x, w
      double precision, intent(in)  :: omega
      double precision, intent(out) :: Iqp,dIqp

      integer ilo,ihi,jlo,jhi
      double precision wmo(nocc+nvir),factor(nocc+nvir)
      double precision factor2(nocc+nvir),dfactor2(nocc+nvir),w2

      integer igl,mynmo,jmo

      double precision,external :: ddot

      call ga_distribution(eps(1),me,jlo,jhi,ilo,ihi)

      mynmo = jhi - jlo + 1
      do jmo=1,mynmo
        factor(jmo) = omega - eig(jmo+jlo-1)
      enddo

      Iqp = 0d0
      dIqp = 0d0
      do igl=1,ngl
        w2 = x(igl)
        call ga_get(eps(igl),jlo,jhi,imo,imo,wmo,mynmo)
        do jmo=1,mynmo
          factor2(jmo) =  factor(jmo)/(factor(jmo)**2 + w2)
          dfactor2(jmo) = 1d0/(factor(jmo)**2 + w2) - 
     $                    2d0*factor2(jmo)**2
        enddo
        Iqp = Iqp - w(igl)*ddot(mynmo,factor2,1,wmo,1)
        dIqp = dIqp - w(igl)*ddot(mynmo,dfactor2,1,wmo,1)
      enddo
      call ga_dgop(mt_dbl,Iqp,1,'+')
      call ga_dgop(mt_dbl,dIqp,1,'+')
      
      end subroutine

c     *****************************************************************

      subroutine gw_buildrn(oo,ov,vv,eig,eia,omega,R_n,imo,npoles,nocc,
     $                      nvir,nri,me)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer,intent(in) :: oo, ov, vv
      integer,intent(in) :: imo, npoles, nocc, nvir, nri, me
      double precision, intent(in) :: omega
      double precision, intent(in), dimension(nocc+nvir) :: eig
      double precision, intent(in), dimension(npoles) :: eia
      double precision, intent(out) :: R_n

      character(len=12), parameter :: pname = 'gw_buildrn: '

      integer g_a, g_b, g_c, g_d, g_e, g_epsilon, k_eps
      integer nmo, jmo, ipole, ilo, ihi, jlo, jhi, info, ld, mynri, idx
      double precision fac, arg, w
      double precision factor(npoles), vector(nri), vector2(nri)

      double precision, parameter :: shift = 0.0001d0
      double precision, external :: ddot
      logical, external :: is_near, nga_create, ga_duplicate, ga_destroy
      logical, external :: ga_create

      nmo = nocc + nvir

      if (.not.nga_create(mt_dbl,1,npoles,'factor',-1,g_a))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,nri,'one',-1,g_d))
     $  call errquit(pname//'can''t create global array',0,GA_ERR) 
      if (.not.ga_duplicate(ov,g_c,'temp'))
     $  call errquit(pname//'can''t duplicate array',0,GA_ERR) 
      if (.not.ga_create(mt_dbl,nri,nri,'epsilon',nri,-1,g_epsilon))
     $  call errquit(pname//'can''t create array',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nri,nri,'Pi',nri,-1,g_b))
     $  call errquit(pname//'can''t create array',0,GA_ERR)


      call ga_distribution(g_epsilon,me,ilo,ihi,jlo,jhi)
      call nga_distribution(g_a,me,ilo,ihi)
      call ga_zero(g_d)
      call ga_add_constant(g_d,1d0)
      call ga_access(g_epsilon,1,nri,jlo,jhi,k_eps,ld)

      mynri = jhi - jlo + 1
      R_n = 0d0

      if (omega.lt.0d0) then
        fac = -1d0
        do jmo=1,nocc
          if (eig(jmo).lt.omega-shift) cycle
          if (is_near(eig(jmo),omega,shift)) fac = -0.5d0

          arg = (omega - eig(jmo))**2
          do ipole=ilo,ihi
            factor(ipole) = eia(ipole)/(eia(ipole)**2 - arg)
          enddo
          call nga_put(g_a,ilo,ihi,factor(ilo),1)
          call ga_sync()
          call ga_copy(ov,g_c)
          call ga_scale_cols(g_c,g_a)
          call ga_dgemm('n','t',nri,nri,npoles,4d0,g_c,ov,0d0,g_b)
          call ga_add_diagonal(g_b,g_d)
          call dft_invdiag(g_b, g_epsilon, nri)

          idx = (imo-1)*nocc + jmo
          call ga_get(oo,1,nri,idx,idx,vector,nri)
          call dgemv('t',nri,mynri,1d0,dbl_mb(k_eps),ld,vector,1,0d0,
     $                vector2(jlo),1)

          vector2(jlo:jhi) = vector2(jlo:jhi) - vector(jlo:jhi)
          w = ddot(mynri,vector(jlo),1,vector2(jlo),1)
          R_n = R_n + fac*w
        enddo
      else if (omega.gt.0d0) then
        fac = 1d0
        do jmo=nocc+1,nmo
          if (eig(jmo).gt.omega+shift) exit
          if (is_near(eig(jmo),omega,shift)) fac = 0.5d0
        enddo
      end if
      call ga_dgop(mt_dbl,R_n,1,'+')

      if ((.not.ga_destroy(g_a)).or.(.not.ga_destroy(g_b)).or.
     $    (.not.ga_destroy(g_c)).or.(.not.ga_destroy(g_d)).or.
     $    (.not.ga_destroy(g_epsilon)))
     $  call errquit(pname//'can''t destroy global array',0,GA_ERR) 
c     
      return
      end subroutine

      logical function is_near(test,targ,tol)
      implicit none
      double precision, intent(in) :: test, targ, tol

      if (abs(test-targ).lt.tol) then
        is_near = .true.
      else
        is_near = .false.
      endif

      return
      end function

      subroutine gw_vmo(oo,ov,vv,vmo,sigma_x,nqp,nmo,nocc,nvir,nri,
     $                  llmo,me,nvqp)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer,intent(in) :: oo, ov, vv, me
      integer,intent(in) :: nqp, nmo, nocc, nvir, nri, llmo, nvqp
      double precision,intent(out),dimension(nmo,nqp) :: vmo
      double precision,intent(out),dimension(nqp) :: sigma_x

      double precision,external :: ddot

      integer ilo,ihi,oolo,oohi,ovlo,ovhi,vvlo,vvhi,idx
      integer l_tmp, tmp
      integer imo,jmo,kmo

      call ga_distribution(oo,me,ilo,ihi,oolo,oohi)
      call ga_distribution(ov,me,ilo,ihi,ovlo,ovhi)
      if (nvqp.gt.0) call ga_distribution(vv,me,ilo,ihi,vvlo,vvhi)

      if (.not.ma_push_get(mt_dbl,nri,'temp',l_tmp,tmp))
     $  call errquit('gw_vmo: could not allocate',0,MA_ERR) 

      vmo = 0d0
      do imo=llmo,llmo+nqp-1
        jmo = imo - llmo + 1
        if (imo.le.nocc) then
          do kmo=1,nocc
            idx = (imo-1)*nocc + kmo
            if (idx.lt.oolo) cycle
            if (idx.gt.oohi) exit
            call ga_get(oo,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
          do kmo=1,nvir
            idx = (imo-1)*nvir + kmo
            if (idx.lt.ovlo) cycle
            if (idx.gt.ovhi) exit
            call ga_get(ov,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo+nocc,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
        else
          do kmo=1,nocc
            idx = (kmo-1)*nvir + imo - nocc
            if (idx.lt.ovlo) cycle
            if (idx.gt.ovhi) exit
            call ga_get(ov,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
          do kmo=1,nvir
            idx = (kmo-1)*nvir + kmo
            if (idx.lt.vvlo) cycle
            if (idx.gt.vvhi) exit
            call ga_get(vv,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo+nocc,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
        endif
      enddo

      if (.not.ma_pop_stack(l_tmp))
     $  call errquit('gw_vmo: could not pop stack',0,MA_ERR) 

      call dgop(mt_dbl,vmo,nqp*nmo,'+')

      do imo=llmo,llmo+nqp-1
        sigma_x(imo) = sum(vmo(1:nocc,imo))
      enddo

      end subroutine

      subroutine gw_buildpi(ov,pi,factor,nri,npoles,ilo,ihi)
      implicit none
      integer,intent(in) :: nri,npoles,ilo,ihi
      double precision,intent(in),dimension(npoles) :: factor
      double precision,intent(inout),dimension(nri,npoles) :: ov
      double precision,intent(out),dimension(nri,nri) :: pi

      integer ipole,iri

      do ipole=1,npoles
        ov(1:nri,ipole) = ov(1:nri,ipole)*factor(ipole)
      enddo

      call dsyrk('l','n',nri,npoles,4d0,ov,nri,0d0,pi,nri)

      do iri=1,nri-1
        pi(iri,iri+1:) = pi(iri+1:,iri)
      enddo

      do iri=ilo,ihi
        pi(iri,iri) = pi(iri,iri) + 1d0
      enddo

      do ipole=1,npoles
        ov(1:nri,ipole) = ov(1:nri,ipole)/factor(ipole)
      enddo
      
      end subroutine
