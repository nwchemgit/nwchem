      subroutine gw_run(pars)
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      character(len=8), parameter :: pname = 'gw_run: '
c
      integer nmo, nri, nocc, nvir, iocc, ivir, iqp, i, j, nri_me
      integer npoles, ipole, lo1, lo2, eviter, ilo, ihi, ibf, lo3
      integer nstates
      integer l_sigmax, k_sigmax
      integer l_sigma, k_sigma
      integer l_eia, k_eia
      integer l_amb, k_amb
      integer l_omega, k_omega
      integer l_xpy, k_xpy
      integer l_scr, k_scr
      integer l_eriia,k_eriia
      integer l_vxc, k_vxc
      integer l_wmn, k_wmn
      integer l_smomega, k_smomega
      integer l_eritmp, k_eritmp
      integer l_eold, k_eold
      integer l_enew, k_enew
      integer l_sign, k_sign
      integer l_eriij, k_eriij
      integer lscr,lo4,ld, istart, iend
      logical main,oactive
c
      double precision sigma_c, dsigma_c, eps_in, eps_out, sigma_old
      double precision z, sigma, func, linear, efermi, sigma_x
c
      integer ilaenv
      external ilaenv
c
c     -----------------------------------------------------------------
c
c     synchronize parallel processes
c
      call ga_sync()
c
c     initialize useful variables
c
      nmo  = pars%nmo(1)
      nocc = pars%nocc(1)
      nvir = pars%nvir(1)
      nri  = nbf_cd
      npoles = nocc*nvir
      lscr = ilaenv(1,'dsytrd','u',npoles,-1,-1,-1)
      lscr = (lscr+2)*npoles
      nri_me = pars%nbf_cd_me
      nstates = pars%nstates
      main = pars%me.eq.0
      oactive = pars%nbf_cd_me .gt. 0
      lo1 = nocc - pars%noqp
      istart = pars%istart
      iend = pars%istart + pars%nbf_cd_me - 1
c
c     Get fermi energy
c
      efermi = (dbl_mb(pars%k_evals+nocc-1) +
     $         dbl_mb(pars%k_evals+nocc))/2.0d0
c
c     Obtain Vxc and store it in Sigma
c
      if (main) write(luout,9010)
 9010 format(10x,'Computing V_xc')

      if (.not.ma_push_get(mt_dbl,2*nmo,'Sigma',l_sigma,k_sigma))
     $  call errquit(pname//'cannot allocate Sigma',0,MA_ERR)

      call gw_vxc(dbl_mb(k_sigma),nmo,pars)
c
c     Obtain Sigma_x
c
      if (main) write(luout,9000)
 9000 format(10x,'Computing Sigma_x')
      
      if (.not.ma_push_get(mt_dbl,nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)
c
      call ga_access(pars%g_moeris,1,nmo*nocc,istart,iend,k_eriij,ld)
      call gw_sigmax(dbl_mb(k_eriij),dbl_mb(k_sigmax),nmo,nri_me,nocc)
c
c     Allocation
c
      if (.not.ma_push_get(mt_dbl,nri_me*npoles,'eri_ia',
     $    l_eriia,k_eriia))
     $  call errquit(pname//'cannot allocate eri_ia',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals old',l_eold,k_eold))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals new',l_enew,k_enew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'signs',l_sign,k_sign))
     $  call errquit(pname//'failed to allocate sign',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,npoles,'omega',l_omega,k_omega))
     $  call errquit(pname//'cannot allocate Omega',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,npoles,'eval_ia',l_eia,k_eia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
c     Shift eigenvalues
c
      do i=0,nmo-1
        dbl_mb(k_eold+i) = dbl_mb(pars%k_evals+i) - efermi
      enddo
c
      eviter = 0
c
  100 continue  ! evGW calculations recompute RPA polarizability
c
      if (pars%evgw.and.main) then
        write(luout,*)
        write(luout,'("  *** G",I1,"W",I1)') eviter,eviter
        write(luout,*)
      endif
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
c     RPA polarizability
c
      if (main) write(luout,9020)
 9020 format(10x,'Computing RPA polarizability')

      if(.not.ma_push_get(mt_dbl,npoles*npoles,'XpY',l_xpy,k_xpy))
     $  call errquit(pname//'cannot allocate XpY',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lscr,'work',l_scr,k_scr))
     $  call errquit('gw_rpa: could not allocate work',0,MA_ERR)
c
      call gw_rpa(dbl_mb(k_eia),dbl_mb(k_omega),
     $  dbl_mb(k_xpy),dbl_mb(k_eriia),dbl_mb(k_scr),npoles,nri_me,lscr)
c
      if (.not.ma_chop_stack(l_xpy))
     $  call errquit(pname//'could not free stack',0,MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,nmo*npoles,'wmn',l_wmn,k_wmn))
     $  call errquit(pname//'failed to allocate wnm',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri_me*nmo,'eritemp',
     $   l_eritmp,k_eritmp))
     $  call errquit(pname//'failed to allocate eritmp',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'sm omega',l_smomega,k_smomega))
     $  call errquit(pname//'failed to allocate small omega',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,4*nmo,'scratch',l_scr,k_scr))
     $  call errquit(pname//'failed to allocate scratch',0,MA_ERR) 
c
  200 continue ! evGW0 calculations recompute the rest
c      
      if (pars%evgw0.and.main) then
        write(luout,*)
        write(luout,'("  *** G",I1,"W0")') eviter
        write(luout,*)
      endif
c      
      if (main) write(luout,9030)
 9030 format(10x,'Solving quasiparticle equations')

      if (main) write(luout,9040)
 9040 format(/,10x,'State',5x,'Energy (eV)',5x,'Z',8x,'Elin (eV)',
     $       /,10x,45('-'))

      dbl_mb(k_sign:k_sign+nmo-1) =
     $  sign(1.0d0,dbl_mb(k_eold:k_eold+nmo-1))
c
      do iqp=1,pars%noqp+pars%nvqp

        call gw_get_wmn(dbl_mb(k_wmn),dbl_mb(k_eritmp),dbl_mb(k_eriia),
     $                  pars%g_moeris,iqp+lo1,npoles,nri_me,nmo,
     $                  pars%istart)
c      
        if (pars%graph) then
        elseif (pars%newton) then
c
          eps_out = dbl_mb(k_eold+lo1+iqp-1)
          sigma_x = dbl_mb(k_sigmax+lo1+iqp-1)
          sigma_old = dbl_mb(k_sigma+lo1+iqp-1)
          do i=1,pars%qpiter
            eps_in = eps_out
            dbl_mb(k_smomega:k_smomega+nmo-1) = eps_in - 
     $                                      dbl_mb(k_eold:k_eold+nmo-1)
c
            call gw_sigmac(sigma_c,dsigma_c,dbl_mb(k_smomega),
     $       dbl_mb(k_wmn),dbl_mb(k_omega),dbl_mb(k_scr),dbl_mb(k_sign),
     $       0.001d0**2,npoles,nri_me,nmo)
c
            z = 1.0d0/(1.0d0-dsigma_c)
            sigma = sigma_c + sigma_x
            
            func = dbl_mb(k_eold+lo1+iqp-1) + sigma - sigma_old - eps_in
            eps_out = eps_in + 0.667*z*func

            if (i.eq.1) linear = eps_in + z*func

            if (abs(eps_out-eps_in).lt.pars%thresh) exit

          enddo
          dbl_mb(k_enew+lo1+iqp-1) = eps_out
          dbl_mb(k_sigma+lo1+iqp-1) = sigma_c + sigma_x
          if (main) write(luout,9050) lo1+iqp,(eps_out+efermi)*ha2ev,z,
     $                                (linear+efermi)*ha2ev     
        endif
      enddo
c
      if (pars%evgw0.or.pars%evgw) then
        call dcopy(nmo,dbl_mb(k_enew),1,dbl_mb(k_eold),1)
        eviter = eviter + 1      
        if (eviter.lt.pars%eviter) then
          if (pars%evgw0) goto 200
          if (pars%evgw) then
            if (.not.ma_chop_stack(l_wmn))
     $        call errquit(pname//'failed to chop stack',0,MA_ERR)
            goto 100
          endif
        endif  
      endif

 9050 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3)
c
      if (.not.ma_chop_stack(l_sigma))
     $  call errquit(pname//'could not free stack',0,MA_ERR)

      end subroutine gw_run
c
c     *****************************************************************
c
      subroutine gw_sigmax(eri,sigmax,nmo,nri,nocc)
c
c     computes the exchange part of the self-energy
c
      implicit none
#include "mafdecls.fh"
      integer,intent(in) ::  nri, nmo, nocc
      double precision,intent(in) :: eri(nocc*nmo,nri)
      double precision,intent(out) :: sigmax(nmo)
c      
      double precision temp
      integer imo,jmo,ij
c
      sigmax(:) = 0.0d0
c
      ij = 0
      do imo=1,nocc
        do jmo=1,nmo
          ij = ij + 1
          sigmax(jmo) = sigmax(jmo) - 
     $    dot_product(eri(ij,:),eri(ij,:))
        enddo
      enddo
      call ga_dgop(mt_dbl,sigmax,nmo,'+')
c
      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_rpa(eia,omega,xpy,eriia,scr,
     $                  npoles,nri,lscr)
c
c     builds and solves for the RPA polarizability
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
      integer,intent(in) :: npoles,nri,lscr
      double precision eia(npoles),amb(npoles),omega(npoles)
      double precision xpy(npoles,npoles),eriia(npoles,nri)
      double precision scr(lscr)
c
      integer iocc,ipole,g_a,g_b,g_c,g_d
      integer info,iri,me,ilo,ihi,k_index,ld
      logical main

      logical ga_create,nga_create
      external ga_create,nga_create

      integer ga_nodeid
      external ga_nodeid
c
      me = ga_nodeid()
      main = me .eq. 0
c      
      if (.not.ga_create(mt_dbl,npoles,npoles,'x+y',npoles,-1,g_c))
     $  call errquit('gw_rpa: can''t create global array',0,GA_ERR)
      if (.not.ga_create(mt_dbl,npoles,npoles,'rpa',npoles,-1,g_a))
     $  call errquit('gw_rpa: can''t create global array',0,GA_ERR)
      if (.not.nga_create(mt_dbl,1,npoles,'scratch',-1,g_b))
     $  call errquit('gw_rpa: can''t create global array',0,GA_ERR)
      if (.not.nga_create(mt_dbl,1,npoles,'eia',-1,g_d))
     $  call errquit('gw_rpa: can''t create global array',0,GA_ERR)
c
      call ga_zero(g_a)
      do ipole=1,npoles
        call dgemv('n',npoles,nri,4d0,eriia,npoles,eriia(ipole,1),
     $              npoles,0d0,omega,1)
        call ga_acc(g_a,1,npoles,ipole,ipole,omega,1,1.0d0)
      enddo
      call ga_sync()

      call nga_distribution(g_d,me,ilo,ihi)
      call nga_put(g_d,ilo,ihi,eia(ilo),1)
      call nga_put(g_b,ilo,ihi,dsqrt(eia(ilo:ihi)),1)
      call ga_add_diagonal(g_a,g_d)
      call ga_destroy(g_d)

      call ga_scale_cols(g_a,g_b)
      call ga_scale_rows(g_a,g_b)
c
      call ga_diag_std(g_a,g_c,omega)
      call ga_destroy(g_a)

      omega = dsqrt(omega)
c
      call ga_scale_rows(g_c,g_b)
      call nga_put(g_b,ilo,ihi,1d0/dsqrt(omega(ilo:ihi)),1)
      call ga_scale_cols(g_c,g_b)
      call ga_destroy(g_b)
c
      call ga_get(g_c,1,npoles,1,npoles,xpy,npoles)

      do iri=1,nri
        scr(1:npoles) = eriia(:,iri)
        call dgemv('t',npoles,npoles,1d0,xpy,npoles,scr,1,0d0,
     $              eriia(1,iri),1)
      enddo

      call ga_destroy(g_c)
c
      return
      end subroutine
c
c     *****************************************************************

      subroutine gw_get_wmn(wnm,erin,erim,g_a,imo,npoles,nri,nmo,istart)
c
c     computes the matrix representation of the screened Coulomb op.
c
      implicit none
#include "mafdecls.fh"      
      integer imo, npoles, nri, nmo, idx, jmo, g_a, ilo, ihi, iri
      integer istart
      double precision wnm(npoles,nmo),erim(npoles,nri),erin(nmo,nri)

      ilo = (imo-1)*nmo + 1
      ihi = imo*nmo

      do iri=istart,istart+nri-1
        call ga_get(g_a,ilo,ihi,iri,iri,erin(1,iri-istart+1),nmo)
      enddo

      call dgemm('n','t',npoles,nmo,nri,1.0d0,erim,npoles,erin,nmo,
     $            0.0d0,wnm,npoles) 

      call ga_dgop(mt_dbl,wnm,npoles*nmo,'+')

      wnm = 2.0d0*wnm**2

      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_sigmac(sigma_c,dsigma_c,omegan,wmn,Omegam,scratch,
     $                     signs,eta2,npoles,nri,nmo)
c 
c     computes the correlation part of the self-energy
c
      implicit none
      integer, intent(in) :: npoles,nri,nmo
      double precision, intent(out) ::  sigma_c, dsigma_c
      double precision, intent(in) :: omegan(nmo), wmn(npoles,nmo),
     $                                Omegam(npoles), eta2, signs(nmo)
      double precision scratch(*)

      integer ipole

      sigma_c = 0.0d0
      dsigma_c = 0.0d0

      do ipole=1,npoles
        scratch(1:nmo) = omegan - Omegam(ipole)*signs
        scratch(nmo+1:2*nmo) = 1.0d0/(scratch(1:nmo)**2 + eta2)
        scratch(2*nmo+1:3*nmo) = scratch(1:nmo)*scratch(nmo+1:2*nmo)
        scratch(3*nmo+1:4*nmo) = scratch(nmo+1:2*nmo) - 2.0d0*
     $                           scratch(2*nmo+1:3*nmo)**2
        sigma_c = sigma_c + dot_product(wmn(ipole,:),
     $                                  scratch(2*nmo+1:3*nmo))
        dsigma_c = dsigma_c + dot_product(wmn(ipole,:),
     $                                    scratch(3*nmo+1:4*nmo))
      enddo 

      return
      end subroutine
