#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

      subroutine gw_cdgw_newton(pars,Sigmax,Sigma,Eref,Eold,Enew,Eia,iW,
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

      logical,target  :: lvars(2)
      logical,pointer :: converged, bracket
      logical main, warning, fixed(maxqp)
      integer isp, totalqp, llmo, iqp, iter, logs, dbls, irms
      double precision timer, timer2, rms6(6), m, thresh
      double precision Ein, Eout, SigmaOld, Sigma_X
      double precision esterror(pars%nmo)
      double precision residual, dresidual, maxradius
      GWTYP :: I, dI, R, dR, SigmaC, dSigmaC

      integer nclusters, clusters(pars%nmo),myll,myul
      integer llqp,ulqp,icluster
      double precision upper,lower,resupper,reslower,const
      double precision values(25),errors(25),radius


      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

      thresh = pars%thresh
      if (eviter.gt.2) thresh = thresh/10d0

      converged => lvars(1)
      bracket   => lvars(2)

      if (eviter.gt.4) then 
        maxradius = 0.010d0
      else
        maxradius = 0.030d0
      endif

      do isp=1,pars%ipol
        warning = .false.

        totalqp = pars%noqp(isp) + pars%nvqp(isp)
        if (totalqp.lt.1) cycle

        esterror(:) = 0d0

        if ((isp.eq.1).and.main) then
          write(luout,9000)
        elseif(main) then
          write(luout,9010)
        endif

        llmo = pars%llmo(isp)
        fixed(:) = .false.

        !
        ! Find eigenvalues tightly clustered. The solution will be
        ! attempted starting from the value closer to the Fermi energy
        ! and using its result as guess for the subsequent orbitals.
        !
        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.05d0)
        call ga_brdcst(1038,nclusters,ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(1038,clusters,
     $                 ma_sizeof(mt_int,nclusters,mt_byte),0)

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

          if (eviter.lt.3) then
            if (myul.lt.ulqp) then
              Eout = Enew(llmo+iqp,isp) + 0.001d0
            elseif (myll.gt.llqp) then
              Eout = Enew(llmo+iqp-2,isp) - 0.001d0
            endif
          endif

          call ga_zero(pars%g_sols)
          bracket = .false.
          rms6 = 1d100
          const = Eref(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = -1d100
          upper = 1d100
c        
          do iter=1,pars%qpiter+5
            timer2 = util_wallsec()

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()
            call gw_cdgw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                      dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                      iWidx(0,iqp,isp),Ein,pars%eta,I,dI,
     $                      pars%ngl,mynmo)
            pars%timer%in = pars%timer%in + util_wallsec() - timer

            call ga_sync()
            timer = util_wallsec()

            call gw_cdgw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)

            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            SigmaC   = I + R
            dSigmaC  = dI + dR


            residual = const + dble(SigmaC) - Ein
            dresidual = (dble(dSigmaC) - 1d0)

            values(iter) = Ein
            errors(iter) = residual

            radius = min(maxradius,0.1d0*maxradius/residual**2)

            ! Find out if we have bracketed the solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,values(iter),errors(iter),
     $                            values(iter-1),errors(iter-1),
     $                            lower,upper,reslower,resupper)
            endif

            ! Check convergence of QP equation
            converged = abs(residual).lt.thresh .or. 
     $                  sum(rms6)/6d0.lt.thresh .or.
     $                  (bracket .and. abs(upper-lower).lt.thresh)

            ! Exit if converged
            call ga_brdcst(1038,lvars,2*logs,0)
            if (converged) then
              Eout = Ein
#ifdef GWDEBUG
              if (main) call gw_newton_debug(residual,eout,ein,lower,
     $                  upper,efermi(isp),util_wallsec()-timer2,
     $                  bracket,llmo+iqp-1,iter)
#endif
              exit
            endif

            !
            ! Take a step
            !
            if (main) then
              call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,thresh,
     $                       radius,iter)

#ifdef GWDEBUG
              call gw_newton_debug(residual,eout,ein,lower,
     $                  upper,efermi(isp),util_wallsec()-timer2,
     $                  bracket,llmo+iqp-1,iter)
#endif
            endif
            call ga_brdcst(1038,eout,ma_sizeof(mt_dbl,1,mt_byte),0)

            ! Save information
            do irms=6,2,-1
              rms6(irms) = rms6(irms-1)
            enddo
            rms6(1) = abs(eout-ein)

          enddo

          ! Save last energy
          Enew(llmo+iqp-1,isp)  = Eout

          ! check convergence info
          if (converged) then
            fixed(iqp) = .true.
            if (bracket) then
              esterror(llmo+iqp-1) = min(upper-lower,abs(residual))
            else
              esterror(llmo+iqp-1) = abs(residual)
            endif
          else
            if (bracket) then
              esterror(llmo+iqp-1) = upper - lower
            else
              esterror(llmo+iqp-1) = abs(residual)
            endif
          endif

          if (llmo+iqp-1.le.pars%nocc(isp)) then
            myul = myul - 1
          else
            myll = myll + 1
          endif

          if (myll.le.myul) goto 100

          if (main) call gw_print_energy(llmo,llqp,ulqp,
     $                       enew(llmo+llqp-1,isp),efermi(isp),
     $                       esterror(llmo+llqp-1),fixed(llqp),
     $                       warning)
        enddo

        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Result did not converge'
          write(luout,*)
        endif
      enddo


 9000 format(/,18x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))
 9010 format(/,18x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))

      end subroutine
