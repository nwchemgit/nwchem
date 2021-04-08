#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

      subroutine gw_cdgw_diis(pars,Sigmax,Sigma,Eref,Eold,Enew,Eia,iW,
     $                          iWidx,Efermi,nri,mynmo,maxqp,maxpoles,
     $                          EViter)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, mynmo, maxqp, maxpoles, EViter
      integer,intent(in) :: iWidx(0:mynmo,maxqp,pars%ipol)
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eref(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eia(maxpoles,pars%ipol)
      GWTYP,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged, bracket, warning, fixed(maxqp)
      integer isp, totalqp, llmo, iqp, iter, logs, dbls
      double precision timer, timer2, rms3(3), m
      double precision Ein, Eout, SigmaOld, Sigma_X, Z
      double precision Linear
      double precision residual, residualOld, dresidual, Etmp
      double precision dsigmacold
      GWTYP :: I, dI, R, dR, SigmaC, dSigmaC

      integer nclusters, clusters(pars%nmo),myll,myul
      integer llqp,ulqp,icluster
      double precision upper,lower,resupper,reslower,const,step

      integer,parameter :: maxdiis = 15
      integer           :: ndiis,idiis
      double precision  :: values(maxdiis),errors(maxdiis)

      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

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
        fixed(:) = .false.

        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.05d0)

        ulqp = 0
        do icluster=1,nclusters

          llqp = ulqp + 1
          ulqp = ulqp + clusters(icluster)

          myul = ulqp
          myll = llqp

 100      continue

          if (llmo+llqp-1.le.pars%nocc(isp)) then
            iqp  = myul
          else
            iqp  = myll
          endif
            
          Eout     = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X  = Sigmax(llmo+iqp-1,isp)

          if (myul.lt.ulqp) then
            Eout = Enew(llmo+iqp,isp) + 4d0*pars%thresh
          elseif (myll.gt.llqp) then
            Eout = Enew(llmo+iqp-2,isp) - 4d0*pars%thresh
          elseif (eviter.eq.0) then
            Eout = Eout + 0.10d0*(sigma_x - sigmaold)
          else
            Eout = Eout + sign(4d0*pars%thresh,Eout)
          endif

          call ga_zero(pars%g_sols)
          bracket = .false.
          rms3 = 1d100
          residualOld = 0d0
          dsigmacold = 0d0
          const = Eref(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = 0d0
          upper = 0d0
          ndiis = 0
c        
          do iter=1,pars%qpiter+5
            timer2 = util_wallsec()

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()
            call gw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                      dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                      iWidx(0,iqp,isp),Ein,I,dI,
     $                      pars%ngl,mynmo)
            pars%timer%in = pars%timer%in + util_wallsec() - timer

            call ga_sync()
            timer = util_wallsec()

            call gw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)

            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            SigmaC   = I + R
            dSigmaC  = dI + dR

            values(iter) = Ein
            errors(iter) = const + sigmac - ein
            residual = errors(iter)

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  sum(rms3)/3d0.lt.pars%thresh 


            ! Exit if converged
            call ga_brdcst(mt_log,converged,logs,0)
            if (converged) then
              Eout = Ein
              exit
            endif

            call gw_diis(errors,values,step,ndiis,1,iter)

            Eout = Ein + step

#ifdef GWDEBUG
            if (main) write(luout,9901) iter,iqp,util_wallsec()-timer2,
     $                                  Ein,Eout
 9901       format('  ## Finished iter',I3,' of qp ',I3,' in ',F12.2,'s'
     $             ,/,'     Ein: ',F12.6,' Eout: ',F12.6)
#endif

            rms3(1) = rms3(2)
            rms3(2) = rms3(3)
            rms3(3) = abs(eout-ein)

          enddo

          ! Save last energy
          Enew(llmo+iqp-1,isp)  = Eout

          ! check convergence info
          if (converged) fixed(iqp) = .true.

          if (llmo+iqp-1.le.pars%nocc(isp)) then
            myul = myul - 1
          else
            myll = myll + 1
          endif

          if (myll.le.myul) goto 100

          if (main) then
            do iqp=llqp,ulqp
              if (fixed(iqp)) then
                write(luout,9020) llmo+iqp-1,
     $                       (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev
              else
                warning = .true.
                write(luout,9030) llmo+iqp-1,
     $                       (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev
              endif
            enddo
          endif

        enddo
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







