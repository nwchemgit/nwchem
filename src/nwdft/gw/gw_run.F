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
      integer lscr,lo4,ld, istart, iend, g_erim, g_wmn
      logical main,oactive,converged
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
      if (main) write(luout,9000)
 9000 format(10x,'Computing V_xc')

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
c
      call gw_sigmax(pars%g_erioo,pars%g_eriov,dbl_mb(k_sigmax),nri,nmo,
     $               nocc,nvir,pars%me)
c
c     Allocation
c
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
      if (.not.ga_duplicate(pars%g_eriov,g_erim,'transformed eri'))
     $  call errquit(pname//'failed to duplicate array',0,GA_ERR)
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
c     RPA polarizability
c
      if (main) write(luout,9020)
 9020 format(10x,'Computing RPA polarizability')
c
      call gw_rpa(pars%g_eriov,g_erim,dbl_mb(k_eia),dbl_mb(k_omega),
     $            npoles,nri)
c
      if (.not.ga_create(mt_dbl,nmo,npoles,'wmn',nmo,-1,g_wmn))
     $  call errquit(pname//'cannot create array',0,GA_ERR)
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

        call gw_get_wmn(g_wmn,g_erim,pars%g_erioo,pars%g_eriov,
     $                  pars%g_erivv,iqp+lo1,npoles,nri,nmo,nocc)
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
     $       g_wmn,dbl_mb(k_omega),dbl_mb(k_scr),dbl_mb(k_sign),
     $       0.001d0**2,npoles,nmo,pars%me)
c
            z = 1.0d0/(1.0d0-dsigma_c)
            sigma = sigma_c + sigma_x
            
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
          if (main) write(luout,9050) lo1+iqp,(eps_out+efermi)*ha2ev,z,
     $                                (linear+efermi)*ha2ev     
        endif
      enddo
c
      if (pars%evgw0.or.pars%evgw) then
        call ga_brdcst(1003,dbl_mb(k_enew),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call ga_brdcst(1004,dbl_mb(k_sigma),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call dcopy(nmo,dbl_mb(k_enew),1,dbl_mb(k_eold),1)
        eviter = eviter + 1      
        if (eviter.lt.pars%eviter) then
          if (pars%evgw0) goto 200
          if (pars%evgw) then
            if (.not.ma_chop_stack(l_smomega))
     $        call errquit(pname//'failed to chop stack',0,MA_ERR)
            if (.not.ga_destroy(g_wmn))
     $        call errquit(pname//'failed to destroy array',0,GA_ERR)
            goto 100
          endif
        endif  
      endif

 9050 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3)
c
      if (.not.ga_destroy(g_wmn))
     $  call errquit(pname//'failed to destroy array',0,GA_ERR)
      if (.not.ga_destroy(g_erim))
     $  call errquit(pname//'failed to destroy array',0,GA_ERR)
      if (.not.ma_chop_stack(l_sigma))
     $  call errquit(pname//'could not free stack',0,MA_ERR)

      end subroutine gw_run
c
c     *****************************************************************
c
      subroutine gw_sigmax(oo,ov,sigmax,nri,nmo,nocc,nvir,me)
c
c     computes the exchange part of the self-energy
c
      implicit none
#include "mafdecls.fh"
      integer,intent(in) :: oo, ov, me
      integer,intent(in) :: nri, nmo, nocc, nvir
      double precision,intent(out) :: sigmax(nmo)
c      
      logical ok
      double precision temp
      integer imo,jmo,ij,ilo,ihi,jlo,jhi,k_tmp,l_tmp

      double precision ddot,ga_ddot_patch
      external ddot,ga_ddot_patch

c
      sigmax(:) = 0.0d0


      do imo=1,nocc
        ilo = (imo-1)*nocc + 1
        ihi = imo*nocc
        sigmax(imo) = -ga_ddot_patch(oo,'n',1,nri,ilo,ihi,
     $                               oo,'n',1,nri,ilo,ihi)
      enddo
c
      call ga_distribution(ov,me,ilo,ihi,jlo,jhi)

      ok =  (.not.ma_push_get(mt_dbl,nri,'temp',l_tmp,k_tmp))

      ij = 0
      do imo=1,nocc
        do jmo=nocc+1,nmo
          ij = ij + 1
          if (ij.lt.jlo) cycle
          if (ij.gt.jhi) exit
          call ga_get(ov,1,nri,ij,ij,dbl_mb(k_tmp),nri)
          sigmax(jmo) = sigmax(jmo) -
     $                  ddot(nri,dbl_mb(k_tmp),1,dbl_mb(k_tmp),1)
        enddo
      enddo

      ok = ma_chop_stack(l_tmp)

      call ga_dgop(mt_dbl,sigmax(nocc+1),nvir,'+')
c
      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_rpa(ov,erim,eia,omega,npoles,nri)
c
c     builds and solves for the RPA polarizability
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
      integer,intent(in) :: ov, erim
      integer,intent(in) :: npoles,nri
      double precision eia(npoles),amb(npoles),omega(npoles)
c
      integer iocc,ipole,g_a,g_b,g_c,g_d
      integer info,iri,me,ilo,ihi,k_index,ld
      logical main

      logical ga_create,nga_create,ga_destroy
      external ga_create,nga_create,ga_destroy

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
      call ga_dgemm('t','n',npoles,npoles,nri,4d0,ov,ov,0d0,g_a)
c      
      call nga_distribution(g_d,me,ilo,ihi)
      call nga_put(g_d,ilo,ihi,eia(ilo),1)
      call nga_put(g_b,ilo,ihi,dsqrt(eia(ilo:ihi)),1)

      call ga_sync()
      call ga_add_diagonal(g_a,g_d)
      call ga_scale_cols(g_a,g_b)
      call ga_scale_rows(g_a,g_b)
c
      call ga_diag_std(g_a,g_c,omega)
      if(.not.ga_destroy(g_a))
     $  call errquit('',0,GA_ERR) 

      omega = dsqrt(omega)
c
      call ga_scale_rows(g_c,g_b)
      call nga_put(g_b,ilo,ihi,1d0/dsqrt(omega(ilo:ihi)),1)

      call ga_sync()
      call ga_scale_cols(g_c,g_b)

      call ga_dgemm('n','n',nri,npoles,npoles,1d0,ov,g_c,0d0,erim)
c
      if ((.not.ga_destroy(g_b)).or.(.not.ga_destroy(g_c)).or.
     $    (.not.ga_destroy(g_d)))
     $  call errquit('',0,GA_ERR) 
c
      return
      end subroutine
c
c     *****************************************************************

      subroutine gw_get_wmn(wmn,erim,oo,ov,vv,imo,npoles,nri,nmo,nocc)
c
c     computes the matrix representation of the screened Coulomb op.
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer,intent(in) :: erim,oo,ov,vv,wmn
      integer,intent(in) :: imo,npoles,nri,nmo,nocc
      integer ilo, ihi, iocc, g_a

      logical ga_create,ga_duplicate,ga_destroy
      external ga_create,ga_duplicate,ga_destroy

      if (.not.ga_create(mt_dbl,nri,nmo,'temp',nri,-1,g_a))
     $  call errquit('gw_get_wmn: cannot create array',0,GA_ERR) 
      call ga_zero(g_a)

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

      call ga_dgemm('t','n',nmo,npoles,nri,1d0,g_a,erim,0d0,wmn)
      if(.not.ga_destroy(g_a))
     $ call errquit('',0,GA_ERR) 

      if(.not.ga_duplicate(wmn,g_a,'temp'))
     $  call errquit('',0,GA_ERR)
      call ga_copy(wmn,g_a)
      call ga_scale(g_a,2d0)
      call ga_elem_multiply(wmn,g_a,wmn)
      if (.not.ga_destroy(g_a))
     $ call errquit('',0,GA_ERR) 

      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_sigmac(sigma_c,dsigma_c,omegan,wmn,Omegam,scratch,
     $                     signs,eta2,npoles,nmo,me)
c 
c     computes the correlation part of the self-energy
c
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: wmn
      integer, intent(in) :: npoles,nmo,me
      double precision, intent(out) ::  sigma_c, dsigma_c
      double precision, intent(in) :: omegan(nmo), 
     $                                Omegam(npoles), eta2, signs(nmo)
      double precision scratch(*)

      integer ipole,ilo,ihi,jlo,jhi

      double precision ddot
      external ddot

      sigma_c = 0.0d0
      dsigma_c = 0.0d0

      call ga_distribution(wmn,me,ilo,ihi,jlo,jhi)
      do ipole=jlo,jhi
        call dcopy(nmo,omegan,1,scratch,1)
        call daxpy(nmo,-Omegam(ipole),signs,1,scratch,1)
        scratch(nmo+1:2*nmo) = 1.0d0/(scratch(1:nmo)**2 + eta2)
        scratch(2*nmo+1:3*nmo) = scratch(1:nmo)*scratch(nmo+1:2*nmo)
        scratch(3*nmo+1:4*nmo) = scratch(nmo+1:2*nmo) - 2.0d0*
     $                           scratch(2*nmo+1:3*nmo)**2
        call ga_get(wmn,1,nmo,ipole,ipole,scratch,nmo)
        sigma_c = sigma_c + ddot(nmo,scratch,1,scratch(2*nmo+1),1)
        dsigma_c = dsigma_c + ddot(nmo,scratch,1,scratch(3*nmo+1),1)
      enddo
      call ga_dgop(mt_dbl,sigma_c,1,'+')
      call ga_dgop(mt_dbl,dsigma_c,1,'+')

      return
      end subroutine
