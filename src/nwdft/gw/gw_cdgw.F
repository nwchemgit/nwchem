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
#include "util.fh"      
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
      integer mynpoles,idx_oo,idx_ov,idx_vv

      integer nprocs
      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4, mynmo
      integer oolo,oohi,jlo,jhi,vvlo,vvhi,k_iw,l_iw,k_iwidx,l_iwidx

      double precision efermi
      double precision I_n, R_n, dI_n, dR_n, dI_np, dI_nm, dR_np, dR_nm
      double precision d2I_n, d2R_n, dummy
      double precision sigma, sigma_c, dsigma_c, d2sigma_c
      double precision eps_in, eps_out, z, func, linear, dfunc, d2func
      double precision sigma_x, sigma_old

      double precision timer_gw,timer_vmo,timer_iw,timer_in,timer_rn
      double precision timer_vxc,timer

      logical converged

      double precision, parameter :: fdstep = 0.001d0
      double precision, external :: gw_efermi
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
      mynpoles = pars%mynpoles
      nprocs = ga_nnodes()

      call ga_distribution(pars%idx_oo,pars%me,oolo,oohi,jlo,jhi)
      call ga_access(pars%idx_oo,oolo,oohi,1,2,idx_oo,ld)
      call ga_access(pars%idx_ov,pars%polelo,pars%polehi,1,2,idx_ov,ld)
      if (nstates.gt.nocc) then
        call ga_distribution(pars%idx_vv,pars%me,vvlo,vvhi,jlo,jhi)
        call ga_access(pars%idx_vv,vvlo,vvhi,1,2,idx_vv,ld)
      else
        idx_vv = 1
      endif
c
c     Obtain Vxc
c
      if (main) write(luout,9000)
 9000 format(10x,'Computing Vxc')

      timer = util_wallsec()

      if (.not.ma_push_get(mt_dbl,2*nmo,'Sigma',l_sigma,k_sigma))
     $  call errquit(pname//'cannot allocate Sigma',0,MA_ERR)

      call gw_vxc(dbl_mb(k_sigma),nmo,pars)

      pars%timer%vxc = pars%timer%vxc + util_wallsec() - timer
c
c     Obtain Sigma_x
c
      if (main) write(luout,9010)
 9010 format(10x,'Computing Sigma_x')

      timer = util_wallsec()

      if (.not.ma_push_get(mt_dbl,nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nqp*nmo,'Vmo',l_vmo,k_vmo))
     $  call errquit(pname//'cannot allocate Vmo',0,MA_ERR)

      call gw_vmo(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $            dbl_mb(k_vmo),dbl_mb(k_sigmax+lo1),nqp,nmo,nocc,
     $            nvir,nri,lo1+1,pars%me,pars%nvqp,int_mb(idx_oo),
     $            int_mb(idx_ov),int_mb(idx_vv))

      pars%timer%vmo = pars%timer%vmo + util_wallsec() - timer
c
c     Allocations
c
      mynmo = nmo/nprocs + 2
      if (.not.ma_alloc_get(mt_dbl,mynmo*nqp*ngl,'W(igl)',l_iw,k_iw))
     $  call errquit(pname//'can''t create iw array',0,MA_ERR)
      if (.not.ma_alloc_get(mt_int,mynmo*nqp+nqp,'iW idx',l_iwidx,
     $  k_iwidx))
     $  call errquit(pname//'can''t create iw index',0,MA_ERR)

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
      efermi = gw_efermi(dbl_mb(pars%k_evals),nocc)
c
c     Shift eigenvalues
c
      call gw_shift(dbl_mb(pars%k_evals),dbl_mb(k_eold),efermi,nmo)
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
      call gw_get_eia(dbl_mb(k_eold),dbl_mb(k_eia),int_mb(idx_ov),
     $                mynpoles)
c
c     Obtain Screened Coulomb matrices on imaginary axis
c
      if (main) write(luout,9030)
 9030 format(10x,'Computing W(iw) on the imaginary grid')
c
      timer = util_wallsec()

      call gw_get_iwidx(int_mb(k_iwidx),nmo,nocc,nvir,mynmo,nqp,lo1+1)

      call gw_buildiw(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $    dbl_mb(k_eia),dbl_mb(pars%k_glx),dbl_mb(k_iw),
     $    int_mb(k_iwidx),dbl_mb(k_vmo),
     $    npoles,nmo,nocc,nvir,nri,pars%ngl,nqp,lo1+1,mynmo)

      pars%timer%iw = pars%timer%iw + util_wallsec() - timer
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
        sigma_x = dbl_mb(k_sigmax+lo1+iqp-1)
        sigma_old = dbl_mb(k_sigma+lo1+iqp-1)
        eps_out = dbl_mb(k_eold+lo1+iqp-1)
        if (eviter.eq.0) eps_out = eps_out + 0.17d0*(sigma_x-sigma_old)
c        
        do i=1,pars%qpiter
          eps_in = eps_out
          timer = util_wallsec()
          call gw_buildin(dbl_mb(k_iw),dbl_mb(k_eold),
     $       dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),int_mb(k_iwidx),
     $       eps_in,I_n,dI_n,fdstep,iqp,npoles,nocc,nvir,nri,pars%ngl,
     $       pars%me,mynmo,nqp)
          pars%timer%in = pars%timer%in + util_wallsec() - timer
          timer = util_wallsec()
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in,R_n,
     $                    lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in+fdstep,
     $                    dR_np,lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          call gw_buildrn(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     $                    dbl_mb(k_eold),dbl_mb(k_eia),eps_in-fdstep,
     $                    dR_nm,lo1+iqp,npoles,nocc,nvir,nri,pars%me)
          pars%timer%rn = pars%timer%rn + util_wallsec() - timer
          sigma_c = R_n + I_n
          dsigma_c = (dR_np - dR_nm)/(2d0*fdstep) + dI_n
          sigma = sigma_c + sigma_x

          z = 1d0/(1d0-dsigma_c)
            
          func = dbl_mb(k_eold+lo1+iqp-1) + sigma - sigma_old - eps_in
          eps_out = eps_in + z*func

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

      if (.not.(ma_free_heap(l_iw).and.ma_free_heap(l_iwidx)))
     $  call errquit(pname//'cannot destroy iW arrays',0,MA_ERR)

c
      end subroutine

c     *****************************************************************

      subroutine gw_buildiw(oo,ov,vv,eia,x,iw,iwidx,vmo,npoles,nmo,nocc,
     $                      nvir,nri,ngl,nqp,llmo,mynmo)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer,intent(in) :: npoles,nmo,nocc,nvir,nri,ngl,nqp,llmo,mynmo
      integer,intent(in) :: oo,ov,vv
      double precision,intent(in) :: eia(npoles),x(ngl)
      double precision,intent(in) :: vmo(nmo,nqp)
      integer,intent(in),dimension(0:mynmo,nqp) :: iwidx
      double precision,intent(out),dimension(mynmo,ngl,nqp) :: iw


      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok
      integer g_a,g_b,g_c,g_d,g_e,g_f
      double precision factor(npoles),buf(npoles)


      logical main
      integer molo,mohi,lmo
      integer ilo,ihi,me,k_ov,ld,mynpoles,igl,ipole,info,k_c
      integer pilo,pihi,k_pi,l_pi,idx,idy,mynri,offset,l_tmp,tmp
      integer oolo,oohi,ovlo,ovhi,vvlo,vvhi,jmo,kmo,imo,nprocs
      double precision w2,wmo(nmo),vector(nri)

      double precision,external :: ddot
      integer,external :: ga_nodeid, ga_nnodes
      logical,external :: nga_create,ga_destroy,ga_create,ga_duplicate

      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0

      if (.not.ma_push_get(mt_dbl,nri*nri,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri*mynmo,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate tmp',0,MA_ERR) 
      call dfill(nri*nri,0d0,dbl_mb(k_pi),1)

      call ga_distribution(oo,me,ilo,ihi,oolo,oohi)
      call ga_distribution(ov,me,ilo,ihi,ovlo,ovhi)
      call ga_access(ov,1,nri,ovlo,ovhi,k_ov,ld)
      if (llmo+nqp-1.gt.nocc) then
        call ga_distribution(vv,me,ilo,ihi,vvlo,vvhi)
      endif
      mynpoles = ovhi-ovlo+1

      do igl=1,ngl
c
        w2 = x(igl)
        do ipole=1,mynpoles
          factor(ipole) = eia(ipole)/(eia(ipole)**2 + w2)
        enddo
c        
        call gw_buildpi('w',dbl_mb(k_ov),dbl_mb(k_pi),factor,
     $                  nri,mynpoles,main)
        call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')
        !
        ! All cores invert the same matrix in order to avoid
        ! communication. This is good for small matrices but
        ! will be very slow for large ones.
        !
        call dpotrf('l',nri,dbl_mb(k_pi),nri,info)
        if (info.ne.0) then
          call errquit(pname//'Cholesky decomposition failed',0,0)
          !call dft_invdiag
        endif
        !
        do imo=llmo,llmo+nqp-1
          jmo = imo - llmo + 1
          if (imo.le.nocc) then

            do kmo=1,iwidx(0,jmo)
              lmo = iwidx(kmo,jmo)
              idx = tmp + (kmo-1)*nri
              if (lmo.le.nocc) then
                idy = ((imo-1)*nocc+lmo-1)/nprocs + oolo
                call ga_get(oo,1,nri,idy,idy,dbl_mb(idx),nri)
              else
                idy = ((imo-1)*nvir+lmo-nocc-1)/nprocs + ovlo
                call ga_get(ov,1,nri,idy,idy,dbl_mb(idx),nri)
              endif
            enddo


          else

            do kmo=1,iwidx(0,jmo)
              lmo = iwidx(kmo,jmo)
              idx = tmp + (kmo-1)*nri
              if (lmo.le.nocc) then
                idy = ((lmo-1)*nvir+imo-nocc-1)/nprocs + ovlo
                call ga_get(ov,1,nri,idy,idy,dbl_mb(idx),nri)
              else
                idy = ((imo-nocc-1)*nvir+lmo-nocc-1)/nprocs + ovlo
                call ga_get(ov,1,nri,idy,idy,dbl_mb(idx),nri)
              endif
            enddo

          endif
c
          call dtrtrs('l','n','n',nri,iwidx(0,jmo),dbl_mb(k_pi),
     $                nri,dbl_mb(tmp),nrii,info)
c
          if (info.ne.0) then
            call errquit(pname//' triangular solve failed',0,0)
          endif
c
          do kmo=1,iwidx(0,jmo)
            idx = tmp + (kmo-1)*nri
            lmo = iwidx(kmo,jmo)
            iw(kmo,igl,jmo) = ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1) +
     $                        vmo(lmo,jmo)
          enddo
c
        enddo
      enddo

      if (.not.ma_chop_stack(l_pi))
     $  call errquit(pname//'can''t pop stack',0,MA_ERR) 
c
      end subroutine

c     *****************************************************************

      subroutine gw_buildin(iw,eig,x,w,iwidx,omega,Iqp,dIqp,fdstep,imo,
     $                      npoles,nocc,nvir,nri,ngl,me,myld,nqp)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: npoles,nocc,nvir,nri,ngl,myld,nqp
      integer,intent(in) :: imo, me
      integer,intent(in),dimension(0:myld,nqp) :: iwidx
      double precision, intent(in) :: fdstep
      double precision, intent(in), dimension(nocc+nvir) :: eig
      double precision, intent(in), dimension(ngl) :: x, w
      double precision, intent(in)  :: omega
      double precision, intent(in), dimension(myld,ngl,nqp) :: iw
      double precision, intent(out) :: Iqp,dIqp

      integer ilo,ihi,jlo,jhi
      double precision wmo(nocc+nvir),factor(nocc+nvir)
      double precision factor2(nocc+nvir),dfactor2(nocc+nvir),w2

      integer igl,mynmo,jmo

      double precision,external :: ddot

      mynmo = iwidx(0,imo)
      do jmo=1,mynmo
        factor(jmo) = omega - eig(iwidx(jmo,imo))
      enddo

      Iqp = 0d0
      dIqp = 0d0
      do igl=1,ngl
        w2 = x(igl)
        do jmo=1,mynmo
          factor2(jmo) =  factor(jmo)/(factor(jmo)**2 + w2)
          dfactor2(jmo) = (w2 - factor(jmo)**2+fdstep**2)/
     $  (((factor(jmo)+fdstep)**2 + w2)*((factor(jmo)-fdstep)**2+w2))
        enddo
        Iqp = Iqp - w(igl)*ddot(mynmo,factor2,1,iw(1,igl,imo),1)
        dIqp = dIqp - w(igl)*ddot(mynmo,dfactor2,1,iw(1,igl,imo),1)
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

      integer mynpoles, root2
      integer root, oolo, oohi, vvlo, vvhi, nrisize, ovlo, ovhi
      integer g_a, g_b, g_c, g_d, g_e, g_epsilon, k_eps, k_sol, nprocs
      integer nmo, jmo, ipole, ilo, ihi, jlo, jhi, info, ld, mynri, idx
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, k_ov, lscr
      double precision fac, arg, w
      double precision factor(npoles), vector(nri), vector2(nri)

      logical main
      double precision, parameter :: shift = 0.0001d0
      double precision, external :: ddot
      logical, external :: is_near, nga_create, ga_duplicate, ga_destroy
      logical, external :: ga_create
      integer,external :: ga_nodeid,ga_nnodes

      main = me.eq.0
      nmo = nocc + nvir
      nprocs = ga_nnodes()

      if (.not.ma_push_get(mt_dbl,nri**2,'pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate pi matri',0,MA_ERR) 
      if (.not.ma_push_get(mt_int,nri,'ipiv',l_ipiv,k_ipiv))
     $  call errquit(pname//'can''t allocate ipiv',0,MA_ERR) 
      call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $            factor,-1,info)
      lscr = int(factor(1))
      if (.not.ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $  call errquit(pname//'can''t allocate scratch',0,MA_ERR) 

      call ga_distribution(oo,me,ilo,ihi,oolo,oohi)
      if (omega.gt.0d0) call ga_distribution(vv,me,ilo,ihi,vvlo,vvhi)
      call ga_distribution(ov,me,ilo,ihi,ovlo,ovhi)
      call ga_access(ov,1,nri,ovlo,ovhi,k_ov,ld)

      mynpoles = ovhi - ovlo + 1

      R_n = 0d0

      if (omega.lt.0d0) then
        do jmo=1,nocc
          if (eig(jmo).lt.omega-shift) cycle
          fac = -1d0
          if (is_near(eig(jmo),omega,shift)) fac = -0.5d0

          arg = (omega - eig(jmo))**2
          do ipole=1,mynpoles
            factor(ipole) = eia(ipole)/(eia(ipole)**2 - arg)
          enddo

          call gw_buildpi('r',dbl_mb(k_ov),dbl_mb(k_pi),factor,
     $                    nri,mynpoles,main)
          call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')

          root = mod((imo-1)*nocc+jmo-1,nprocs)
          if (me.ne.root) cycle

          idx = ((imo-1)*nocc+jmo-1)/nprocs + oolo
          call ga_get(oo,1,nri,idx,idx,vector,nri)
          vector2(:) = vector(:)
          !
          ! Only one core is working at a time in order to
          ! reduce commuication. Will be very slow for large
          ! systems
          !
          call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $                 dbl_mb(k_scr),lscr,info)
          call dsytrs('l',nri,1,dbl_mb(k_pi),nri,int_mb(k_ipiv),vector2,
     $                 nri,info)     
          vector2(:) = vector2(:) - vector(:)

          w = ddot(nri,vector,1,vector2,1)
          R_n = R_n + fac*w
        enddo

      else if (omega.gt.0d0) then

        do jmo=nocc+1,nmo
          if (eig(jmo).gt.omega+shift) exit

          fac = 1d0
          if (is_near(eig(jmo),omega,shift)) fac = 0.5d0

          arg = (omega - eig(jmo))**2
          do ipole=1,mynpoles
            factor(ipole) = eia(ipole)/(eia(ipole)**2 - arg)
          enddo

          call gw_buildpi('r',dbl_mb(k_ov),dbl_mb(k_pi),factor,
     $                    nri,mynpoles,main)
          call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')

          root = mod((imo-nocc-1)*nvir+jmo-nocc-1,nprocs)
          if (me.ne.root) cycle

          idx = ((imo-nocc-1)*nvir+jmo-nocc-1)/nprocs + vvlo
          call ga_get(vv,1,nri,idx,idx,vector,nri)
          vector2(:) = vector(:)

          call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $                 dbl_mb(k_scr),lscr,info)
          call dsytrs('l',nri,1,dbl_mb(k_pi),nri,int_mb(k_ipiv),vector2,
     $                 nri,info)     
          vector2(:) = vector2(:) - vector(:)

          w = ddot(nri,vector,1,vector2,1)
          R_n = R_n + fac*w
        enddo

      end if
      call ga_dgop(mt_dbl,R_n,1,'+')

      if (.not.ma_chop_stack(l_pi))
     $   call errquit(pname//'can''t destroy global array',0,MA_ERR)
c     
      return
      end subroutine

c     *****************************************************************

      logical function is_near(test,targ,tol)
      implicit none
      double precision, intent(in) :: test, targ, tol

      if (abs(test-targ).le.tol) then
        is_near = .true.
      else
        is_near = .false.
      endif

      return
      end function

c     *****************************************************************

      subroutine gw_vmo(oo,ov,vv,vmo,sigma_x,nqp,nmo,nocc,nvir,nri,
     $                  llmo,me,nvqp,idx_oo,idx_ov,idx_vv)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer,intent(in) :: oo, ov, vv, me
      integer,intent(in) :: nqp, nmo, nocc, nvir, nri, llmo, nvqp
      integer,intent(in),dimension(nocc**2,2) :: idx_oo
      integer,intent(in),dimension(nvir*nocc,2) :: idx_ov
      integer,intent(in),dimension(nvir*max(1,nvqp),2) :: idx_vv
      double precision,intent(out),dimension(nmo,nqp) :: vmo
      double precision,intent(out),dimension(nqp) :: sigma_x

      double precision,external :: ddot

      integer ilo,ihi,oolo,oohi,ovlo,ovhi,vvlo,vvhi,idx
      integer l_tmp, tmp, nprocs
      integer imo,jmo,kmo

      integer,external :: ga_nnodes

      nprocs = ga_nnodes()

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
            if (me.ne.mod((imo-1)*nocc+kmo-1,nprocs)) cycle
            idx = ((imo-1)*nocc+kmo-1)/nprocs + oolo
            call ga_get(oo,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-1)*nvir+kmo-1,nprocs)) cycle
            idx = ((imo-1)*nvir+kmo-1)/nprocs + ovlo
            call ga_get(ov,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo+nocc,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo

        else

          do kmo=1,nocc
            if (me.ne.mod((kmo-1)*nvir+imo-nocc-1,nprocs)) cycle
            idx = ((kmo-1)*nvir+imo-nocc-1)/nprocs + ovlo
            call ga_get(ov,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-nocc-1)*nvir+kmo-1,nprocs)) cycle
            idx = ((imo-nocc-1)*nvir+kmo-1)/nprocs + vvlo
            call ga_get(vv,1,nri,idx,idx,dbl_mb(tmp),nri)
            vmo(kmo+nocc,jmo) = -ddot(nri,dbl_mb(tmp),1,dbl_mb(tmp),1)
          enddo

        endif
      enddo

      if (.not.ma_pop_stack(l_tmp))
     $  call errquit('gw_vmo: could not pop stack',0,MA_ERR) 

      call dgop(mt_dbl,vmo,nqp*nmo,'+')

      do imo=1,nqp
        sigma_x(imo) = sum(vmo(1:nocc,imo))
      enddo

      end subroutine

c     *****************************************************************

      subroutine gw_buildpi(which,ov,pi,factor,nri,npoles,main)
      implicit none
#include "mafdecls.fh"      
      logical,intent(in) :: main
      integer,intent(in) :: nri,npoles
      character(len=1),intent(in) :: which
      double precision,intent(in),dimension(npoles) :: factor
      double precision,intent(inout),dimension(nri,npoles) :: ov
      double precision,intent(out),dimension(nri,nri) :: pi

      double precision, dimension(nri) :: vector, scaling
      integer ipole,iri,ma_dbl,l_scr,k_scr,maxpoles,mynpoles,llpole
      integer ulpole,nbatch,ibatch
      logical ok

      if (which.eq.'w') then
        do ipole=1,npoles
          ov(1:nri,ipole) = ov(1:nri,ipole)*dsqrt(factor(ipole))
        enddo
        call dsyrk('l','n',nri,npoles,4d0,ov,nri,0d0,pi,nri)
        do ipole=1,npoles
          ov(1:nri,ipole) = ov(1:nri,ipole)/dsqrt(factor(ipole))
        enddo
      else if (which.eq.'r') then

        ma_dbl = ma_inquire_avail(mt_dbl)/nri
        nbatch = npoles/ma_dbl + 1
        maxpoles = min(npoles/nbatch + 1,npoles)
        ok = ma_alloc_get(mt_dbl,nri*maxpoles,'scratch',l_scr,k_scr)

        call dcopy(nri*maxpoles,ov,1,dbl_mb(k_scr),1)
        do ipole=1,maxpoles
          call dscal(nri,factor(ipole),dbl_mb(k_scr+(ipole-1)*nri),1)
        enddo
        call dsyr2k('l','n',nri,maxpoles,2d0,ov,nri,dbl_mb(k_scr),nri,
     $               0d0,pi,nri)

        do ibatch=2,nbatch
          llpole = (ibatch-1)*maxpoles + 1
          ulpole = min(ibatch*maxpoles,npoles)
          mynpoles = ulpole - llpole + 1
          call dcopy(nri*mynpoles,ov(1,llpole),1,dbl_mb(k_scr),1)
          do ipole=1,mynpoles
            call dscal(nri,factor(llpole+ipole-1),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call dsyr2k('l','n',nri,mynpoles,2d0,ov(1,llpole),nri,
     $                 dbl_mb(k_scr),nri,1d0,pi,nri)     
        enddo

        ok = ma_free_heap(l_scr)

      endif

      if (main) then
        do iri=1,nri
          pi(iri,iri) = pi(iri,iri) + 1d0
        enddo
      endif


      end subroutine

c     *****************************************************************

      double precision function gw_efermi(eig,nocc)
      implicit none
      integer,intent(in) :: nocc
      double precision,intent(in),dimension(*) :: eig
      gw_efermi = (eig(nocc) + eig(nocc+1))/2d0
      return
      end function

      subroutine gw_shift(eig,shifted,efermi,nmo)
      implicit none
      integer,intent(in) :: nmo
      double precision,intent(in) :: efermi,eig(nmo)
      double precision,intent(out) :: shifted(nmo)
      integer imo
      do imo=1,nmo
        shifted(imo) = eig(imo) - efermi
      enddo
      return
      end subroutine

      subroutine gw_get_eia(eig,eia,idx,npoles)
      implicit none
      integer,intent(in) :: npoles
      integer,intent(in),dimension(npoles,2) :: idx
      double precision,intent(in),dimension(*) :: eig
      double precision,intent(out),dimension(npoles) :: eia
      integer ipole,imo,jmo
      do ipole=1,npoles
        imo = idx(ipole,1)
        jmo = idx(ipole,2)
        eia(ipole) = eig(jmo) - eig(imo)
      enddo
      return
      end subroutine

      subroutine gw_get_iwidx(idx,nmo,nocc,nvir,mynmo,nqp,llmo)
      implicit none
      integer,intent(in) :: nmo, nocc, nvir, mynmo, nqp, llmo
      integer,intent(out),dimension(0:mynmo,nqp) :: idx

      integer me, imo, jmo, kmo, lmo, nprocs
      integer,external :: ga_nodeid, ga_nnodes

      me = ga_nodeid()
      nprocs = ga_nnodes()

      do imo=llmo,llmo+nqp-1
        jmo = imo - llmo + 1
        lmo = 0

        if (imo.le.nocc) then
          do kmo=1,nocc
            if (me.ne.mod((imo-1)*nocc+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        else
          do kmo=1,nocc
            if (me.ne.mod((kmo-1)*nvir+imo-nocc-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-nocc-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        endif
        idx(0,jmo) = lmo
      enddo
      return
      end subroutine
