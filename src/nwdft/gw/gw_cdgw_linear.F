      subroutine gw_cdgw_linear(pars,Sigmax,Sigma,Eref,Eold,Enew,Eia,iW,
     $                          iWidx,efermi,nri,mynmo,maxqp,maxpoles)

      implicit none
#include "errquit.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, mynmo, maxqp, maxpoles
      integer,intent(in) :: iWidx(0:mynmo,maxqp,pars%ipol)
      double precision,intent(in) :: efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eref(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eia(maxpoles,pars%ipol)
      double precision,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)

      logical main, converged
      integer isp, totalqp, llmo, iqp
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual
      double precision quadratic, dquadratic, d2quadratic
      double precision quadraticstep, linearstep, step
      double precision I,dI,d2I,R,dR,d2R
      double precision Rp,Rm

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

          Eout     = Eref(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)
          residualOld = 1d100 
c        
          Ein = Eout

          call ga_sync()
          timer = util_wallsec()
          call gw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                    dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                    iWidx(0,iqp,isp),Ein,pars%eta,I,dI,
     $                    pars%ngl,mynmo)
          call ga_sync()
          pars%timer%in = pars%timer%in + util_wallsec() - timer

          timer = util_wallsec()
          call gw_buildrn(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                    llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                    pars%ipol,isp)
          call ga_sync()
          pars%timer%rn = pars%timer%rn + util_wallsec() - timer


          SigmaC   = I + R
          dSigmaC  = dI + dR

          Z = min(1d0,1d0/(1d0-dSigmaC))

          residual = Eold(llmo+iqp-1,isp) + 
     $               SigmaC + Sigma_X - SigmaOld - Ein

          linearstep    = Z*residual

          step = linearstep

          Eout = Ein + linearstep

          call ga_brdcst(mt_dbl,Eout,ma_sizeof(mt_dbl,1,mt_byte),0)

          Enew(llmo+iqp-1,isp)  = Eout
c          Sigma(llmo+iqp-1,isp) = SigmaC + Sigma_X

          if (main) then
            write(luout,9020) llmo+iqp-1,(Eout+efermi(isp))*ha2ev
          endif
        enddo
      enddo

 9000 format(/,15x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       /,10x,25('-'))
 9010 format(/,15x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       /,10x,25('-'))
 9020 format(10x,I3,3x,F12.3)

      end subroutine
