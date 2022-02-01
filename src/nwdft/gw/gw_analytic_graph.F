      subroutine gw_analytic_graph(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                             efermi,nri,maxpoles,Eviter)

      implicit none
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, maxpoles, eviter
      double precision,intent(in) :: efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(inout) :: Enew(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles*pars%ipol)

      logical main, converged, bracket, warning, fixed
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw,dbls,logs,iloc,ipoints,npoints, npoles
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC
      double precision d2SigmaC, dSigmaC, radius
      double precision residual, residualOld, dresidual, d2residual
      double precision Etmp, lower, upper, resupper, reslower
      double precision rms3(3), const, nspace, der, interval
      double precision resmat(8192),sigmamat(8192),omegamat(8192)

      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)
      pars%ecgm = 0d0

c     Create screened-Coulomb array      
      npoles = sum(pars%npoles(1:pars%ipol))
      if (.not.ga_create(mt_dbl,pars%nmo,npoles,'wmn',pars%nmo,0,g_wmn))
     $   call errquit('can''t create Wmnarray',0,GA_ERR)

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

c         Obtain screened-Coulomb matrix elements           
          timer = util_wallsec()
          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)
          pars%timer%wm = pars%timer%wm - timer + util_wallsec()

c         Compute Galiitski-Migdal correlation energy
          if (llmo+iqp-1.le.pars%nocc(isp)) then
            call gw_analytic_ecgm(pars%ecgm,Eold(llmo+iqp-1,isp),Eold,
     $                          omega,g_wmn,pars%eta**2,
     $                          npoles,pars%nmo,plus,pars%me)
          endif
c          
          ! Try just six steps and then use graphical solver
          do iter=1,6

            Ein = Eout

            call ga_sync()

c           Compute self-energy element            
            timer = util_wallsec()
            call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,Eold,
     $           omega,g_wmn,pars%eta**2,npoles,pars%nmo,plus,pars%me)
            pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC

            radius = min(0.1d0,0.1d0/residual**2)

            ! Find out if we have bracketed the solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,ein,residual,etmp,residualold,
     $                            lower,upper,reslower,resupper)
            endif

 1234       format(5(F15.8))

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  (bracket.and.abs(upper-lower).lt.pars%thresh)


            ! Exit if converged
            call ga_brdcst(mt_log,converged,logs,0)
            if (converged) then
              Eout = Ein
              exit
            endif

            ! Take next newton step
            call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,pars%thresh,
     $                       radius,iter)


            ! Save information
            etmp = ein
            residualold = residual

          enddo

          ! Try graphical solver
          nspace   = sign(1d-5,eout)
          npoints  = 8192
          fixed = .false.


          omegamat(1) = eout - sign(0.040d0,eout)
          do ipoints=2,npoints
            omegamat(ipoints) = omegamat(ipoints-1) + nspace
          enddo

          timer = util_wallsec()
          call gw_analytic_sigmac_graph(sigmamat,omegamat,Eold,
     $      omega,g_wmn,0.003d0**2,npoles,pars%nmo,plus,pars%me,npoints)
          pars%timer%sigmac = pars%timer%sigmac - timer + util_wallsec()
          resmat(:) = const + sigmamat(:) - omegamat(:)

          do ipoints=1,npoints

            ! If residual is already small, the solution has been found
            if (abs(resmat(ipoints)).lt.pars%thresh) then
              enew(llmo+iqp-1,isp) = omegamat(ipoints)
              sigma(llmo+iqp-1,isp) = sigmamat(ipoints) + sigma_x
              fixed = .true.
              exit
            endif

            ! if residual just changed sign, bisect solution
            if (ipoints.gt.1) then
              if (resmat(ipoints-1)*resmat(ipoints).lt.0d0) then
                enew(llmo+iqp-1,isp) = omegamat(ipoints) - nspace/2d0
                sigma(llmo+iqp-1,isp) = sigma_x +   
     $             (sigmamat(ipoints)+sigmamat(ipoints-1))/2d0
                fixed = .true.
                exit
              endif
            endif

          enddo

          if (main.and.fixed) then
             write(luout,9020) llmo+iqp-1,
     $                         (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                         1d-5*ha2ev/2d0
           elseif(main) then
             warning = .true.
             iloc = minloc(abs(resmat),dim=1)
             enew(llmo+iqp-1,isp) = omegamat(iloc)
             sigma(llmo+iqp-1,isp) = sigmamat(iloc) + sigma_x
             write(luout,9030) llmo+iqp-1,
     $                         (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                          abs(resmat(iloc))*ha2ev
          endif
        enddo
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Graphical solution was not found'
          write(luout,*)
        endif
      enddo

c     Destroy GAs      
      if (.not.ga_destroy(g_wmn))
     $  call errquit('can''t destroy Wmn',0,GA_ERR)   
      if ( .not.pars%evgw0 ) then
        if (.not.ga_destroy(pars%g_erim))
     $    call errquit('can''t destroy ERIm',0,GA_ERR)   
      endif

 9000 format(/,15x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,36('-'))
 9010 format(/,15x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,36('-'))
 9020 format(10x,I3,3x,F12.3,4x,F12.3)
 9030 format(10x,I3,3x,F12.3,4x,F12.3,' *')

      end
