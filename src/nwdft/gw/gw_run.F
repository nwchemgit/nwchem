      subroutine gw_run(pars)
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"
#include "util.fh"      
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
      integer lo4,ld, istart, iend, g_erim, g_wmn
      logical main,oactive,converged
c
      double precision sigma_c, dsigma_c, eps_in, eps_out, sigma_old
      double precision z, sigma, func, linear, efermi(2), sigma_x

      double precision cputime,walltime,timer,timer_wm,timer_sigmac
      double precision timer_vxc,timer_sigmax,timer_rpa
c
      double precision, external :: gw_efermi

      integer maxpoles,maxnpoles,isp,lSize

      logical ok
c
c     -----------------------------------------------------------------
c
c     synchronize parallel processes
c
      call ga_sync()
c
c     initialize useful variables
c
      nmo  = pars%nmo
      nri  = nbf_cd
      nstates = pars%nstates(1)
      main = pars%me.eq.0
      maxpoles = maxval(pars%mynpoles(1:ipol))
      maxnpoles = maxval(pars%npoles(1:ipol))
c
c     Get fermi energy
c
      efermi = gw_efermi(dbl_mb(pars%k_evals),pars%nocc(1))
c
c     Obtain Vxc and store it in Sigma
c
      if (main) write(luout,9000)
 9000 format(10x,'Computing V_xc ... ')
 9001 format(F12.1,'s')     

      timer = util_wallsec()
c      
      if (.not.ma_push_get(mt_dbl,2*nmo,'Sigma',l_sigma,k_sigma))
     $  call errquit(pname//'cannot allocate Sigma',0,MA_ERR)

      call gw_vxc(dbl_mb(k_sigma),nmo,pars)

      if (.not.ma_verify_allocator_stuff()) stop ' VXC'

      pars%timer%vxc = pars%timer%vxc + util_wallsec() - timer
c
c     Obtain Sigma_x
c
      if (main) write(luout,9010)
 9010 format(10x,'Computing Sigma_x')
      
      timer = util_wallsec()

      if (.not.ma_push_get(mt_dbl,ipol*nmo,'Sigma_x',l_sigmax,k_sigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)
c
      call gw_sigmax(pars,dbl_mb(k_sigmax),nri,nmo)

      if (.not.ma_verify_allocator_stuff()) stop ' Sigmax'
      pars%timer%sigmax = pars%timer%sigmax + util_wallsec() - timer
c
c     Allocation
c
      if (.not.ma_push_get(mt_dbl,ipol*nmo,'evals old',l_eold,k_eold))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,ipol*nmo,'evals new',l_enew,k_enew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,ipol*maxnpoles,'omega',l_omega,
     $                    k_omega))
     $  call errquit(pname//'cannot allocate Omega',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,ipol*maxpoles,'eval_ia',l_eia,k_eia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
c     Shift eigenvalues
c
      call gw_shift(dbl_mb(pars%k_evals),dbl_mb(k_eold),efermi,nbf_ao,
     $              pars%nmo,ipol) 
      if (.not.ma_verify_allocator_stuff()) stop ' Shift'
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
      do isp=1,ipol
        lSize = (isp-1)*maxpoles
        call gw_get_eia(dbl_mb(k_eold+(isp-1)*pars%nmo),
     $                  dbl_mb(k_eia+lSize),pars%nocc(isp),
     $                  pars%nvir(isp))
      enddo
      if (.not.ma_verify_allocator_stuff()) stop ' Eia'
c
c     RPA polarizability
c
      if (main) write(luout,9020)
 9020 format(10x,'Computing RPA polarizability')

      timer = util_wallsec()
c
      call gw_analytic_rpa(pars,dbl_mb(k_eia),dbl_mb(k_omega),maxpoles,
     $                     maxnpoles,nri)

      if (.not.ma_verify_allocator_stuff()) stop ' RPA'
      pars%timer%rpa = pars%timer%rpa + util_wallsec() - timer
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

      if (pars%graph) then
        call gw_analytic_graph(pars,dbl_mb(k_sigmax),dbl_mb(k_sigma),
     $                              dbl_mb(k_eold),dbl_mb(k_enew),
     $                              dbl_mb(k_omega),efermi,nri,
     $                              maxnpoles)
      else if (pars%newton) then
        call gw_analytic_newton(pars,dbl_mb(k_sigmax),dbl_mb(k_sigma),
     $                               dbl_mb(k_eold),dbl_mb(k_enew),
     $                               dbl_mb(k_omega),efermi,nri,
     $                               maxnpoles,eviter)
      else
        call gw_analytic_linear(pars,dbl_mb(k_sigmax),dbl_mb(k_sigma),
     $                               dbl_mb(k_eold),dbl_mb(k_enew),
     $                               dbl_mb(k_omega),efermi,nri,
     $                               maxnpoles)
      endif

c
      if (pars%evgw0.or.pars%evgw) then
        call gw_scissor(dbl_mb(k_Enew),dbl_mb(k_Eold),pars%nmo,
     $                  pars%nocc,pars%nvir,pars%noqp,pars%nvqp,
     $                  pars%ipol,pars%me)

        call ga_brdcst(mt_dbl,dbl_mb(k_enew),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call ga_brdcst(mt_dbl,dbl_mb(k_sigma),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call dcopy(nmo,dbl_mb(k_enew),1,dbl_mb(k_eold),1)
        eviter = eviter + 1      
        if (eviter.lt.pars%eviter) then
          if (pars%evgw0) goto 200
          if (pars%evgw)  goto 100
        endif  
      endif

 9050 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3)
c
      if (.not.ma_chop_stack(l_sigma))
     $  call errquit(pname//'could not free stack',0,MA_ERR)
c

      end subroutine gw_run
