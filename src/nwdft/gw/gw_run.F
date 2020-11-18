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
      integer nmo, nri, nocc, nvir, iocc, ivir, iqp, i, j
      integer npoles, ipole, lo1, lo2, eviter
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
      integer lscr
c
      double precision sigma_c, dsigma_c, eps_in, eps_out, sigma_old
      double precision z, sigma, func, linear, efermi, sigma_x
c
      integer ilaenv
      external ilaenv
c
c     -----------------------------------------------------------------
c
      nmo  = pars%nmo(1)
      nocc = pars%nocc(1)
      nvir = pars%nvir(1)
      nri  = nbf_cd
      npoles = nocc*nvir
      lscr = ilaenv(1,'dsytrd','u',npoles,-1,-1,-1)
      lscr = (lscr+2)*npoles
c
c     Allocation
c
      if (.not.ma_push_get(mt_dbl,nri*npoles,'eri_ia',l_eriia,k_eriia))
     $  call errquit(pname//'cannot allocate eri_ia',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'Vxc',l_vxc,k_vxc))
     $  call errquit(pname//'cannot allocate Vxc',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,2*nmo,'Sigma',l_sigma,k_sigma))
     $  call errquit(pname//'cannot allocate Sigma',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals old',l_eold,k_eold))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'evals new',l_enew,k_enew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'signs',l_sign,k_sign))
     $  call errquit(pname//'failed to allocate sign',0,MA_ERR)
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
c     Copy occupied-virtual block of ERIs into a contiguous block
c 
      lo1 = pars%k_3ceri + npoles
      lo2 = k_eriia
      do i=0,nocc-1
        call dcopy(nvir*nri,dbl_mb(lo1),1,dbl_mb(lo2),1)
        lo1 = lo1 + nmo*nri - i*nri
        lo2 = lo2 + nvir*nri
      enddo
c
c     Set first quasiparticle state index
c
      lo1 = nocc - pars%noqp
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
c     Copy Vxc into old Sigma
c
      call dcopy(nmo,dbl_mb(k_vxc),1,dbl_mb(k_sigma),1)

      if(.not.ma_push_get(mt_dbl,npoles,'omega',l_omega,k_omega))
     $  call errquit(pname//'cannot allocate Omega',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,npoles,'eval_ia',l_eia,k_eia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
      eviter = 0
c
  100 continue  ! evGW calculations recompute RPA polarizability
c
      if (pars%evgw) then
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
      write(luout,9020)
 9020 format(10x,'Computing RPA polarizability')

      if(.not.ma_push_get(mt_dbl,npoles,'AmB',l_amb,k_amb))
     $  call errquit(pname//'cannot allocate AmB',0,MA_ERR) 
      if(.not.ma_push_get(mt_dbl,npoles*npoles,'XpY',l_xpy,k_xpy))
     $  call errquit(pname//'cannot allocate XpY',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lscr,'work',l_scr,k_scr))
     $  call errquit('gw_rpa: could not allocate work',0,MA_ERR)
c
      call gw_rpa(dbl_mb(k_eia),dbl_mb(k_amb),dbl_mb(k_omega),
     $  dbl_mb(k_xpy),dbl_mb(k_eriia),dbl_mb(k_scr),npoles,nri,lscr)
c
      if (.not.ma_chop_stack(l_amb))
     $  call errquit(pname//'could not free stack',0,MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,nmo*npoles,'wmn',l_wmn,k_wmn))
     $  call errquit(pname//'failed to allocate wnm',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,nri*nmo,'eritemp',l_eritmp,k_eritmp))
     $  call errquit(pname//'failed to allocate eritmp',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo,'sm omega',l_smomega,k_smomega))
     $  call errquit(pname//'failed to allocate small omega',0,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,4*nmo,'scratch',l_scr,k_scr))
     $  call errquit(pname//'failed to allocate scratch',0,MA_ERR) 
c
  200 continue ! evGW0 calculations recompute the rest
c      
      if (pars%evgw0) then
        write(luout,*)
        write(luout,'("  *** G",I1,"W0")') eviter
        write(luout,*)
      endif
c      
      write(luout,9030)
 9030 format(10x,'Solving quasiparticle equations')

      write(luout,9040)
 9040 format(/,10x,'State',5x,'Energy (eV)',5x,'Z',8x,'Elin (eV)',
     $       /,10x,45('-'))

      dbl_mb(k_sign:k_sign+nmo-1) =
     $  sign(1.0d0,dbl_mb(k_eold:k_eold+nmo-1))
c
      do iqp=1,pars%noqp+pars%nvqp

        call gw_get_wmn(dbl_mb(k_wmn),dbl_mb(k_eritmp),dbl_mb(k_eriia),
     $                  dbl_mb(pars%k_3ceri),iqp+lo1,npoles,nri,nmo)
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
     $       0.001d0**2,npoles,nri,nmo)
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
          write(luout,9050) lo1+iqp,(eps_out+efermi)*27.2114d0,z,
     $                                (linear+efermi)*27.2114d0     
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
      if (.not.ma_chop_stack(l_eriia))
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
      integer,intent(in) ::  nri, nmo, nocc
      double precision,intent(in) :: eri(nri,*)
      double precision,intent(out) :: sigmax(nmo)
c      
      double precision temp
      integer imo,jmo,ij
c
      sigmax(:) = 0.0d0
c
      ij = 0
      do imo=1,nmo
        do jmo=imo,nmo
          ij = ij + 1
          temp = dot_product(eri(:,ij),eri(:,ij))
          if (jmo.le.nocc) sigmax(imo) = sigmax(imo) - temp
          if (imo.eq.jmo) cycle
          if (imo.le.nocc) sigmax(jmo) = sigmax(jmo) - temp
        enddo
      enddo
c
      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_rpa(eia,amb,omega,xpy,eriia,scr,
     $                  npoles,nri,lscr)
c
c     builds and solves for the RPA polarizability
c
      implicit none
      integer,intent(in) :: npoles,nri,lscr
      double precision eia(npoles),amb(npoles),omega(npoles)
      double precision xpy(npoles,npoles),eriia(nri,npoles)
      double precision scr(lscr)
c
      integer iocc,ipole
      integer info,iri
c
      call dsyrk('u','t',npoles,nri,4.0d0,eriia,nri,0.0d0,xpy,npoles)
c
      amb = dsqrt(eia)
c      
      do ipole=1,npoles
        xpy(ipole,ipole) = xpy(ipole,ipole) + eia(ipole)
        xpy(1:ipole,ipole) = xpy(1:ipole,ipole)*amb(1:ipole)*amb(ipole)
      enddo
c

      call dsyev('v','u',npoles,xpy,npoles,omega,scr,lscr,info)
c
      omega = dsqrt(omega)
c
      do ipole=1,npoles
        xpy(1:npoles,ipole) = xpy(1:npoles,ipole)*amb(1:npoles)/
     $                        dsqrt(omega(ipole))
      enddo
c
      do iri=1,nri
        scr(1:npoles) = eriia(iri,:)
        call dgemv('t',npoles,npoles,1.0d0,xpy,npoles,scr,1,0.0d0,
     $              scr(npoles+1),1)
        eriia(iri,:) = scr(npoles+1:2*npoles)
      enddo
c
      return
      end subroutine
c
c     *****************************************************************

      subroutine gw_get_wmn(wnm,erin,erim,eris,imo,npoles,nri,nmo)
c
c     computes the matrix representation of the screened Coulomb op.
c
      implicit none
      integer imo, npoles, nri, nmo, idx, jmo
      double precision wnm(npoles,nmo),erim(nri,npoles),erin(nri,nmo)
      double precision eris(nri,*)

      do jmo=1,nmo
        if (imo.lt.jmo) then
          idx = imo + (jmo*(jmo-1))/2
        else
          idx = jmo + (imo*(imo-1))/2
        endif
        erin(:,jmo) = eris(:,idx)
      enddo

      call dgemm('t','n',npoles,nmo,nri,1.0d0,erim,nri,erin,nri,0.0d0,
     $            wnm,npoles) 

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
