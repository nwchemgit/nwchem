      subroutine gw_analytic_newton(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                          Efermi,nri,maxpoles,EViter)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"      
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri,maxpoles,EViter
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged, bracket, warning
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw,dbls,logs
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC
      double precision d2SigmaC, dSigmaC
      double precision residual, residualOld, dresidual, d2residual
      double precision Etmp, lower, upper, resupper, reslower
      double precision rms3(3), const

      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)
      pars%ecgm = 0d0

      do isp=1,pars%ipol
        totalqp = pars%noqp(isp) + pars%nvqp(isp)
        if (totalqp.lt.1) cycle

        if ((isp.eq.1).and.main) then
          write(luout,9000)
        elseif(main) then
          write(luout,9010)
        endif

        llmo = pars%llmo(isp)
        warning = .false.

        if (.not.ga_create(mt_dbl,pars%nmo,pars%npoles(isp),'wmn',
     $                     pars%nmo,0,g_wmn))
     $    call errquit('can''t create array',0,GA_ERR)

        do imo=1,pars%nmo
          if (Eold(imo,isp).gt.0d0) then
            plus = imo - 1
            exit
          endif
        enddo

        do iqp=1,totalqp

          Eout     = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)

          bracket = .false.
          rms3 = 1d100
          residualOld = 0d0
          const = Eold(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = 0d0
          upper = 0d0

          if (EViter.eq.0) Eout = Eout + 0.075d0*(Sigma_X - SigmaOld)

          timer = util_wallsec()
          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)
          pars%timer%wm = pars%timer%wm - timer + util_wallsec()

          if (llmo+iqp-1.le.pars%nocc(isp)) then
            call gw_analytic_ecgm(pars%ecgm,Eold(llmo+iqp-1,isp),Eold,
     $                          omega(1,isp),g_wmn,0.003d0**2,
     $                          pars%npoles(isp),pars%nmo,plus,pars%me)
          endif
c        
          do iter=1,pars%qpiter

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()
            call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,Eold,
     $                              omega(1,isp),g_wmn,
     $                              0.003d0**2,pars%npoles(isp),
     $                              pars%nmo,plus,pars%me)
            pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC


            ! Find out if we have bracketed the solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,ein,residual,etmp,residualold,
     $                            lower,upper,reslower,resupper)
            endif

 1234       format(5(F15.8))

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  sum(rms3)/3d0.lt.pars%thresh .or.
     $                  (bracket .and. abs(upper-lower).lt.pars%thresh)


            ! Exit if converged
            call ga_brdcst(mt_log,converged,logs,0)
            if (converged) then
              Eout = Ein
              exit
            endif

            ! Take next newton step
            call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,pars%thresh)


            ! Save information
            etmp = ein
            residualold = residual

            rms3(1) = rms3(2)
            rms3(2) = rms3(3)
            rms3(3) = abs(eout-ein)

          enddo
          Enew(llmo+iqp-1,isp)  = Eout
          Sigma(llmo+iqp-1,isp) = SigmaC + Sigma_X

          if (main) then
            if (converged) then
              write(luout,9020) llmo+iqp-1,
     $                         (Eout+efermi(isp))*ha2ev
            else
              warning = .true.
              write(luout,9030) llmo+iqp-1,
     $                         (Eout+efermi(isp))*ha2ev
            endif
          endif
        enddo
        if (.not.(ga_destroy(g_wmn)))
     $    call errquit('can''t destroy Wmn',0,GA_ERR)   
        if ( .not.pars%evgw0 ) then
          if (.not.(ga_destroy(pars%g_erim(isp))))
     $      call errquit('can''t destroy ERIm',0,GA_ERR)   
        endif
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Result did not converge'
          write(luout,*)
        endif
      enddo

 9000 format(/,15x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       /,10x,25('-'))
 9010 format(/,15x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       /,10x,25('-'))
 9020 format(10x,I3,7x,F8.3)
 9030 format(10x,I3,7x,F8.3,' *')

      end subroutine
