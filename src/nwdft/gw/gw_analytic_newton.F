      subroutine gw_analytic_newton(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                          Efermi,nri,maxpoles,EViter)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri,maxpoles,EViter
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual
      double precision quadratic, dquadratic, d2quadratic
      double precision quadraticstep, linearstep, step

      logical,external :: ga_create, ga_destroy

      main = pars%me .eq. 0

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
          if (Eold(imo,isp).gt.0d0) then
            plus = imo - 1
            exit
          endif
        enddo

        do iqp=1,totalqp

          Eout     = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)
          residualOld = 1d100 

          if (EViter.eq.0) Eout = Eout + 0.075d0*(Sigma_X - SigmaOld)

          timer = util_wallsec()

          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)
          if (.not.ma_verify_allocator_stuff()) stop ' WMN'

          pars%timer%wm = pars%timer%wm - timer + util_wallsec()
c        
          do iter=1,pars%qpiter

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()

            call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,Eold,
     $                              omega(1,isp),g_wmn,
     $                              0.001d0**2,pars%npoles(isp),
     $                              pars%nmo,plus,pars%me)
            if (.not.ma_verify_allocator_stuff()) stop ' Sigmac'

            pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


            Z = 1d0/(1d0-dSigmaC)

            residual = Eold(llmo+iqp-1,isp) + 
     $                 SigmaC + Sigma_X - SigmaOld - Ein
            quadratic = residual**2
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC
            dquadratic = 2d0*residual*dresidual
            d2quadratic = 2d0*dresidual**2 + 2d0*residual*d2residual

            linearstep    = Z*residual
            quadraticstep = -dquadratic/d2quadratic

            if (d2quadratic.gt.1d-3) then
              step = sign(min(0.15d0,abs(quadraticstep)),quadraticstep)
              if (abs(step).lt.0.1d0 .and. abs(residual).gt.0.1d0) then
                step = 0.3d0*linearstep
              endif
            else
              step = 0.6180d0*linearstep
            endif

            Eout = Ein + step

            if (iter.eq.1) Linear = Ein + linearstep

            residualOld = residual

            if (main) converged = abs(residual).lt.pars%thresh
            call ga_brdcst(mt_log,converged,mitob(1),0)
            call ga_brdcst(mt_dbl,Eout,mdtob(1),0)
            if (converged) exit
          enddo
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
        if (.not.(ga_destroy(g_wmn)))
     $    call errquit('can''t destroy Wmn',0,GA_ERR)   
        if (pars%evgw) then
          if (.not.(ga_destroy(pars%g_erim(isp))))
     $      call errquit('can''t destroy ERIm',0,GA_ERR)   
        endif
      enddo

 9000 format(/,31x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9010 format(/,31x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9020 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3,3x,F8.3)

      end subroutine
