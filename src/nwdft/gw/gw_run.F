      subroutine gw_run(pars)
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      character(len=8), parameter :: pname = 'gw_run: '
c
      integer nmo, nri, nocc, nvir, iocc, ivir
      integer npoles, ipole, lo1, lo2
      integer l_sigmax, k_sigmax
      integer l_eia, k_eia
      integer l_amb, k_amb
      integer l_omega, k_omega
      integer l_xpy, k_xpy
      integer l_scr, k_scr
      integer l_eriia,k_eriia
      integer lscr
c
      integer ilaenv
      external ilaenv
c
c     -----------------------------------------------------------------
c
      nmo  = pars%nmo(1)
      nocc = pars%nocc(1)
      nvir = pars%nvir(1)
      nri  = pars%nbf_ri
      npoles = nocc*nvir
      lscr = ilaenv(1,'dsytrd','u',npoles,-1,-1,-1)
      lscr = (lscr+2)*npoles
c
c     Eigenvalue difference
c
      if (.not.ma_push_get(mt_dbl,npoles,'eval_ia',l_eia,k_eia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
      ipole = 0
      do iocc=0,nocc-1
        do ivir=0,nvir-1
          dbl_mb(k_eia+ipole) = dbl_mb(pars%k_evals+ivir) - 
     $                          dbl_mb(pars%k_evals+iocc)
          ipole = ipole + 1
        enddo
      enddo
c
c     Copy occupied-virtual block of ERIs into a contiguous block
c 
      if (.not.ma_push_get(mt_dbl,nri*npoles,'eri_ia',l_eriia,k_eriia))
     $  call errquit(pname//'cannot allocate eri_ia',0,MA_ERR)
c
      lo1 = pars%k_3ceri + nocc*nvir
      lo2 = k_eriia
      do iocc=1,nocc
        call dcopy(nvir*nri,dbl_mb(lo1),1,dbl_mb(lo2),1)
        lo1 = lo1 + nvir*nri + (nocc-iocc+1)*nri
        lo2 = lo2 + nvir*nri
      enddo
c
c     Obtain Sigma_x
c
      write(luout,9000)
 9000 format(10x,'Computing Sigma_x')

      if (.not.ma_push_get(mt_dbl,pars%nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)
      call gw_sigmax(dbl_mb(pars%k_3ceri),dbl_mb(k_sigmax),nmo,nri,nocc)
c
c     Obtain Vxc
c
      call gw_vxc(pars)
c
      if(.not.ma_push_get(mt_dbl,npoles,'omega',l_omega,k_omega))
     $  call errquit(pname//'cannot allocate Omega',0,MA_ERR) 
c
  100 continue  ! evGW calculations recompute RPA polarizability
c
c     RPA polarizability
c
c
      write(luout,9010)
 9010 format(10x,'Computing RPA polarizability')

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

  200 continue ! evGW0 calculations recompute the rest
c


      if (.not.ma_chop_stack(l_eia))
     $  call errquit(pname//'could not free stack',0,MA_ERR)

      end subroutine gw_run
c
      subroutine gw_sigmax(eri,sigmax,nmo,nri,nocc)
      implicit none
      integer nri, nmo, nocc
      double precision eri(nri,*),sigmax(nmo)
      double precision temp
c
      integer iocc,imo,jmo,ij
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


      subroutine gw_rpa(eia,amb,omega,xpy,eriia,scr,
     $                  npoles,nri,lscr)
      implicit none
      integer npoles,nri,lscr
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
