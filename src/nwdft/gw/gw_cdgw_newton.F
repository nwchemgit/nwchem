      subroutine gw_cdgw_newton(pars,Sigmax,Sigma,Eold,Enew,Eia,iW,
     $                          iWidx,Efermi,nri,mynmo,maxqp,maxpoles,
     $                          EViter)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, mynmo, maxqp, maxpoles, EViter
      integer,intent(in) :: iWidx(0:mynmo,maxqp,pars%ipol)
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eia(maxpoles,pars%ipol)
      double precision,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged
      integer isp, totalqp, llmo, iqp, iter
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual
      double precision quadratic, dquadratic, d2quadratic
      double precision quadraticstep, linearstep, step
      double precision I,dI,d2I,R,dR,d2R

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

        do iqp=1,totalqp

          Eout     = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)
          residualOld = 1d100 

          if (EViter.eq.0) Eout = Eout + 0.10d0*(Sigma_X - SigmaOld)
c        
          do iter=1,pars%qpiter

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()
            call gw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                      dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                      iWidx(0,iqp,isp),Ein,I,dI,d2I,
     $                      pars%fdstep,pars%ngl,mynmo)
            pars%timer%in = pars%timer%in + util_wallsec() - timer

            call ga_sync()
            timer = util_wallsec()
            call gw_buildrn(pars,Eold(1,isp),Eia,Ein,R,dR,d2R,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            SigmaC   = I + R
            dSigmaC  = dI + dR
            d2SigmaC = d2I + d2R

            Z = min(1d0,1d0/(1d0-dSigmaC))

            residual = Eold(llmo+iqp-1,isp) + 
     $                 SigmaC + Sigma_X - SigmaOld - Ein
            quadratic = residual**2
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC
            dquadratic = 2d0*residual*dresidual
            d2quadratic = 2d0*d2residual*dresidual**2

            linearstep    = Z*residual
            quadraticstep = -dquadratic/d2quadratic

            step = 0.6180d0*linearstep

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
      enddo

 9000 format(/,31x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9010 format(/,31x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9020 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3,3x,F8.3)

      end subroutine
