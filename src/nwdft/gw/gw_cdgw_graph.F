      subroutine gw_cdgw_graph(pars,Sigmax,Sigma,Eref,Eold,Enew,Eia,iW,
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
      double precision,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged, bracket, warning, fixed(pars%nmo), done
      integer llqp, ulqp, ipoints, icluster, qplimit, npoints
      integer isp, totalqp, llmo, iqp, iter, logs, dbls
      double precision timer, rms3(3), const, nspace
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual,Etmp
      double precision linearstep, step, der
      double precision I,dI,R,dR,a,b,fa,fb

      integer nclusters, clusters(pars%nmo),myll,myul
      double precision upper,lower,resupper,reslower,interval
      double precision resmat(100,pars%nmo)
      double precision omega(100), output(2,pars%nmo)


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

        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.025d0)

        ulqp = 0
        do icluster=1,nclusters

          llqp = ulqp + 1
          ulqp = ulqp + clusters(icluster)

          myul = ulqp
          myll = llqp
          fixed(myll:myul) = .false.

 100      continue

          !obtain a partial solution for one state in the cluster
          if (llmo+llqp-1.le.pars%nocc(isp)) then
            iqp  = myul
          else
            iqp  = myll
          endif

          Eout = Eold(llmo+iqp-1,isp)
          SigmaOld = Sigma(llmo+iqp-1,isp)
          Sigma_X = Sigmax(llmo+iqp-1,isp)

          !use previous solution in cluster as guess
          if (myul.lt.ulqp) then
            Eout = Enew(llmo+iqp,isp)
          elseif (myll.gt.llqp) then
            Eout = Enew(llmo+iqp-2,isp)
          elseif (eviter.eq.0) then
            Eout = Eout + 0.075d0*(sigma_x - sigmaold)
c            Eout = Eout + sign(5d0*pars%thresh,Eout)
          else
            Eout = Eout + sign(1.2d0*pars%thresh,Eout)
          endif

          !Trying to enhance reproducibility by rounding-off to 4 digits
          eout = nint(eout*1d3)/1d3

          !initialize variables
          call ga_zero(pars%g_sols)
          bracket = .false.
          converged = .false.
          rms3(:) = 1d100
          residualold = 0d0

          !fixed contribution
          const = Eref(llmo+iqp-1,isp) + Sigma_X -  SigmaOld

          ! we only need a guess
          do iter=1,5
            Ein = Eout

            call ga_sync()

            ! Contribution from the imaginary axis
            timer = util_wallsec()
            call gw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                      dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                      iWidx(0,iqp,isp),Ein,I,dI,
     $                      pars%ngl,mynmo)
            call ga_sync()
            pars%timer%in = pars%timer%in + util_wallsec() - timer

            ! Contribution from the poles enclosed by the contour
            timer = util_wallsec()
            if (iter.eq.1) then
              call gw_buildrn(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            else
              call gw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            endif
            call ga_sync()
            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            ! Correlation part of the self energy
            SigmaC   = I + R
            dSigmaC  = dI + dR

            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0
            z = -1d0/dresidual

            !bracket solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,ein,residual,etmp,residualold,
     $                            a,b,fa,fb)
            endif

            !check convergence
            converged = (bracket .and. b-a.lt.pars%thresh) .or.
     $                  (abs(residual).lt.pars%thresh) .or.
     $                  (sum(rms3)/3d0.lt.pars%thresh)
            call ga_brdcst(mt_log,converged,logs,0)


            if (converged) then
              eout = ein
              exit
            endif

            call gw_donewton(ein,eout,residual,dresidual,a,b,fa,fb,
     $                       bracket,pars%thresh)

            !save information
            residualold = residual
            etmp = ein
            rms3(3) = rms3(2)
            rms3(2) = rms3(1)
            rms3(1) = abs(eout-ein)

          enddo

          ! Skip graphical solution 
          if (converged) then
            enew(llmo+iqp-1,isp) = eout
            fixed(iqp) = .true.

            if (llmo+iqp-1.le.pars%nocc(isp)) then
              myul = myul - 1
            else
              myll = myll + 1
            endif

            if (myll.gt.myul) then
              goto 200
            else
              goto 100
            endif

          endif

          !Decide number of points
          nspace = sign(0.0015d0,eout)
          interval = maxval(eold(llmo+myll-1:llmo+myul-1,isp)) -
     $               minval(eold(llmo+myll-1:llmo+myul-1,isp)) + 0.007d0
          npoints = min(int(interval/abs(nspace))+1,70)

          omega(1) = eout - sign(0.0015d0,eout)
          do ipoints=1,npoints
            if (ipoints.gt.1)
     $         omega(ipoints) = omega(ipoints-1) + nspace


            timer = util_wallsec()
            call gw_buildin_graph(iW(1,1,myll,isp),Eold(1,isp),
     $        dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),iwidx(0,myll,isp),
     $        omega(ipoints),output,pars%ngl,mynmo,myul-myll+1)
            resmat(ipoints,myll:myul) = output(1,1:myul-myll+1)
            call ga_sync()
            pars%timer%in = pars%timer%in + util_wallsec() - timer


            timer = util_wallsec()
            call gw_buildrn_graph(pars,Eold(1,isp),Eia,omega(ipoints),
     $        output,llmo+myll-1,maxpoles,nri,pars%me,pars%nmo,
     $        pars%ipol,isp,myul-myll+1)
            call ga_sync()
            resmat(ipoints,myll:myul) = resmat(ipoints,myll:myul) +
     $        output(1,1:clusters(icluster)) - omega(ipoints)
            pars%timer%rn = pars%timer%rn + util_wallsec() - timer


            done = .true.
            do iqp=myll,myul
              if (fixed(iqp)) cycle

              ! Finish up building the residual
              resmat(ipoints,iqp) = resmat(ipoints,iqp) +
     $          Eref(llmo+iqp-1,isp) + sigmax(llmo+iqp-1,isp) -
     $          sigma(llmo+iqp-1,isp)

              ! If residual is small, solution has been found
              if (abs(resmat(ipoints,iqp)).lt.pars%thresh) then
                fixed(iqp) = .true.
                enew(llmo+iqp-1,isp) = omega(ipoints) 
                cycle
              endif

              ! See if there is a fixed point solution
              if (ipoints.gt.1) then
                if(resmat(ipoints-1,iqp)*resmat(ipoints,iqp).lt.0d0)then
                  fixed(iqp) = .true.
                  der=(resmat(ipoints,iqp)-resmat(ipoints-1,iqp))/nspace
                  enew(llmo+iqp-1,isp) = omega(ipoints) - 
     $                                   resmat(ipoints,iqp)/der
                  cycle
                endif
              endif

              done = .false.
            enddo
            if (done) exit
          enddo

 200      continue

          if (.not.main) cycle

          do iqp=llqp,ulqp
            if (fixed(iqp)) then
              write(luout,9020) llmo+iqp-1,
     $                          (enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev
            else
              enew(llmo+iqp-1,isp) = omega( 
     $                minloc(abs(resmat(1:npoints,iqp)),dim=1) )
              warning = .true.
              write(luout,9030) llmo+iqp-1,
     $                          (enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev
            endif
          enddo

        enddo

        !print warning if a QP had no graphical solution
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Graphical solution not found.'
          write(luout,*) '   Energy with smallest residual is shown.'
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

