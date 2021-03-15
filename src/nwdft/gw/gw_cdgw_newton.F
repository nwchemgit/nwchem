      subroutine gw_cdgw_newton(pars,Sigmax,Sigma,Eref,Eold,Enew,Eia,iW,
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
      double precision,intent(in) :: Eref(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eia(maxpoles,pars%ipol)
      double precision,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged, bracket, warning
      integer isp, totalqp, llmo, iqp, iter
      double precision timer
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual,Etmp
      double precision quadratic, dquadratic, d2quadratic, dquadraticold
      double precision quadraticstep, linearstep, step
      double precision I,dI,d2I,R,dR,d2R

      integer idiis
      integer,parameter :: nbrack = 3
      double precision upper,lower,resupper,reslower

      logical,external :: ga_create, ga_destroy

      main = pars%me .eq. 0


      do isp=1,pars%ipol
        warning = .false.

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
          residualOld = 0d0

          if (eviter.eq.0) then
            Eout = Eout + 0.30d0*(Sigma_X - SigmaOld)
          else
            Eout = Eout + sign(0.001d0,Eout)
          endif

          call ga_zero(pars%g_sols)
          bracket = .false.
c        
          do iter=1,pars%qpiter+5

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

            if (iter.eq.1) then
              call gw_buildrn(pars,Eold(1,isp),Eia,Ein,R,dR,d2R,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            else
              call gw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,d2R,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            endif

            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            SigmaC   = I + R
            dSigmaC  = dI + dR
            d2SigmaC = d2I + d2R

            Z = 1d0/(1d0-dSigmaC)

            residual = Eref(llmo+iqp-1,isp) + 
     $                 SigmaC + Sigma_X - SigmaOld - Ein
            dresidual = dSigmaC - 1d0

            if (.not.bracket) then
              if (residual*residualold.lt.0d0) then
                bracket = .true.
                if (ein.gt.etmp) then
                  upper = ein
                  lower = etmp
                  resupper = residual
                  reslower = residualold
                else
                  upper = etmp
                  lower = ein
                  resupper = residualold
                  reslower = residual
                endif
              if (main) write(*,'(4F16.8)')lower,upper,reslower,resupper
              endif
            else
              if (residual*resupper.lt.0d0) then
                lower = ein
                reslower = residual
              elseif (residual*reslower.lt.0d0) then
                upper = ein
                resupper = residual
              endif
              if (main) write(*,'(4F16.8)')lower,upper,reslower,resupper
            endif


            if (abs(dresidual).lt.1d-8) then
              linearstep = 0d0
            else
              linearstep = Z*residual
            endif

            if (iter.eq.1) Linear = Ein + linearstep

            step = 0.61803d0*linearstep

            converged = abs(residual).lt.pars%thresh .or. 
     $                  abs(step/Ein).lt.1d-14 .or. 
     $                  (bracket .and. abs(upper-lower).lt.pars%thresh)

            call ga_brdcst(mt_log,converged,mitob(1),0)
            call ga_brdcst(mt_dbl,Eout,mdtob(1),0)
            if (converged) exit

            if (iter.gt.3) step = linearstep

            Eout = Ein + step

            if (bracket) then
              if (eout.gt.upper) eout = upper - (upper-lower)/1.61803d0
              if (eout.lt.lower) eout = lower + (upper-lower)/1.61803d0
            endif

            if (main) write(*,*) eout

            Etmp = Ein
            residualOld = residual
            dquadraticold = dquadratic

          enddo
          Enew(llmo+iqp-1,isp)  = Eout
c          Sigma(llmo+iqp-1,isp) = SigmaC + Sigma_X

          if (main) then
            if (bracket .and. abs(residual).gt.pars%thresh .and. 
     $          converged) then
              warning = .true.
              write(luout,9030) llmo+iqp-1,
     $                         (Eout+efermi(isp))*ha2ev,
     $                          z,
     $                         (linear+efermi(isp))*ha2ev,
     $                          residual*ha2ev
            else
              write(luout,9020) llmo+iqp-1,
     $                         (Eout+efermi(isp))*ha2ev,
     $                          z,
     $                         (linear+efermi(isp))*ha2ev,
     $                          residual*ha2ev
            endif
          endif
        enddo
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Residual is still large but solution'
          write(luout,*) '   was bracketed within 10 meV.'
        endif
      enddo


 9000 format(/,31x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9010 format(/,31x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Z',8x,'Elin (eV)',6x,'Res.',/,10x,55('-'))
 9020 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3,3x,F8.3)
 9030 format(10x,I3,6x,F8.3,6x,F6.3,6x,F8.3,3x,F8.3,' *')

      end subroutine

      subroutine gw_diis(error,values,step,ndiis)
      implicit none
      integer,intent(inout) :: ndiis
      double precision,intent(in) :: error(ndiis),values(ndiis)
      double precision,intent(out) :: step
      double precision :: a(ndiis+1,ndiis+1),w(ndiis+1),b(ndiis+1)
      double precision :: z(ndiis+1,ndiis+1)

      integer :: ipiv(ndiis+1)
      double precision :: work(3*(ndiis+1))

      integer idiis,jdiis,info

      a(1:ndiis,ndiis+1) = 1d0
      a(ndiis+1,1:ndiis) = 1d0
      a(ndiis+1,ndiis+1) = 0d0

      do idiis=1,ndiis
        do jdiis=1,ndiis
          a(jdiis,idiis) = error(jdiis)*error(idiis)
        enddo
      enddo

      b(1:ndiis) = 0d0
      b(ndiis+1) = 1d0

      call dsyev('v','l',ndiis+1,a,ndiis+1,w,work,3*(ndiis+1),info)

      if (abs(w(idiis)).lt.1d-3) then
        ndiis = 0
        write(*,*) 'ill-defined'
        return
      endif

      a(1:ndiis,ndiis+1) = 1d0
      a(ndiis+1,1:ndiis) = 1d0
      a(ndiis+1,ndiis+1) = 0d0

      do idiis=1,ndiis
        do jdiis=1,ndiis
          a(jdiis,idiis) = error(jdiis)*error(idiis)
        enddo
      enddo

      call dgesv(ndiis+1,1,a,ndiis+1,ipiv,b,ndiis+1,info)

      step = 0d0
      do idiis=1,ndiis
        step = step + b(idiis)*values(idiis)
      enddo

      end subroutine







