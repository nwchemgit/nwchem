      subroutine gw_analytic_linear(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                      Efermi,nri,maxpoles)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri,maxpoles
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw,l_scr,k_scr,dbls
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual
      double precision quadratic, dquadratic, d2quadratic
      double precision quadraticstep, linearstep, step

      if (.not.ma_alloc_get(mt_dbl,5*pars%nmo,'scratch',l_scr,k_scr))
     $  call errquit('can''t allocate scratch',0,MA_ERR)

      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)

      do isp=1,pars%ipol
        totalqp = pars%noqp(isp) + pars%nvqp(isp)
        if (totalqp.lt.1) cycle

        if ((isp.eq.1).and.main) then
          write(luout,9000)
        elseif(main) then
          write(luout,9010)
        endif

        llmo = pars%llmo(isp)

        if (.not.ga_create(mt_dbl,pars%nmo,pars%npoles(isp),'wmn',
     $                     pars%nmo,0,g_wmn))
     $    call errquit('can''t create array',0,GA_ERR)

        do imo=1,pars%nmo
          if (Eold(imo,isp).lt.0d0) then
            plus = imo
            exit
          endif
        enddo

        do iqp=1,totalqp

          Eout     = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)
          residualOld = 1d100 

          timer = util_wallsec()

          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)

          pars%timer%wm = pars%timer%wm - timer + util_wallsec()
c        
          Ein = Eout

          call ga_sync()
          timer = util_wallsec()

          call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,Eold,
     $                            omega(1,isp),g_wmn,
     $                            0.001d0**2,pars%npoles(isp),
     $                            pars%nmo,plus,pars%me)

          pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


          Z = min(1d0,1d0/(1d0-dSigmaC))

          residual = Eold(llmo+iqp-1,isp) + 
     $               SigmaC + Sigma_X - SigmaOld - Ein

          linearstep    = Z*residual

          step = linearstep

          Eout = Ein + step

          Linear = Ein + linearstep

          residualOld = residual

          call ga_brdcst(mt_dbl,Eout,dbls,0)

          Enew(llmo+iqp-1,isp)  = Eout
          Sigma(llmo+iqp-1,isp) = SigmaC + Sigma_X

          if (main) then
            write(luout,9020) llmo+iqp-1,
     $                       (Eout+efermi(isp))*ha2ev,
     $                        z,
     $                       (linear+efermi(isp))*ha2ev,
     $                        residual*ha2ev
          endif
        enddo
        if (.not.ga_destroy(g_wmn))
     $    call errquit('can''t destroy array',0,GA_ERR)   
      enddo
      if (.not.ma_free_heap(l_scr))
     $  call errquit('can''t deallocate srcatch',0,MA_ERR) 

 9000 format(/,31x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9010 format(/,31x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9020 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3,3x,F8.3)

      end subroutine
