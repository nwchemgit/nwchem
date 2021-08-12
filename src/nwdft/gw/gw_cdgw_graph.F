#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

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
      type(gw_params_t)            :: pars
      integer,intent(in)           :: nri,mynmo,maxqp,maxpoles,EViter
      integer,intent(in)           :: iWidx(0:mynmo,maxqp,pars%ipol)
      double precision,intent(in)  :: Efermi(pars%ipol)
      double precision,intent(in)  :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in)  :: Eref(pars%nmo,pars%ipol)
      double precision,intent(in)  :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in)  :: Eia(maxpoles,pars%ipol)
      GWTYP,intent(in)  :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(in)  :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical,target  :: lvars(3)
      logical,pointer :: converged,skipped,bracket
      logical main, warning, fixed(pars%nmo), done
      integer llqp, ulqp, ipoints, icluster, qplimit, npoints
      integer isp, totalqp, llmo, iqp, iter, logs, dbls
      double precision timer, timer2, rms3(3), const, nspace
      double precision Ein, Eout, SigmaOld, Sigma_X
      double precision residual, dresidual, der
      double precision values(20),errors(20),esterror(pars%nmo),radius
      GWTYP :: I,dI,R,dR,SigmaC,dSigmaC
      GWTYP :: resmat(100,pars%nmo)

      integer nclusters, clusters(pars%nmo),myll,myul
      double precision upper,lower,resupper,reslower,interval
      double precision omega(100), output(2,pars%nmo)


      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

      converged => lvars(1)
      bracket   => lvars(2)
      skipped   => lvars(3)

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
     $                       0.020d0)

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
          endif

          !Trying to enhance reproducibility by rounding-off to 8 digits
          eout = nint(eout*1d8)/1d8

          !initialize variables
          call ga_zero(pars%g_sols)
          bracket = .false.
          converged = .false.
          rms3(:) = 1d100

          !fixed contribution
          const = Eref(llmo+iqp-1,isp) + Sigma_X -  SigmaOld

          ! we only need a guess
          do iter=1,6
            timer2 = util_wallsec()

            Ein = Eout
            skipped = .false.

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
            call gw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            call ga_sync()
            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            ! Correlation part of the self energy
            SigmaC   = I + R
            dSigmaC  = dI + dR

            residual = const + dble(SigmaC) - Ein
            dresidual = (dble(dSigmaC) - 1d0)

            values(iter) = Ein
            errors(iter) = residual

            radius = min(0.1d0,0.1d0/residual**2)

            ! Find out if we might have skipped a solution
            if ((iter.gt.1) .and. (.not.bracket) .and. 
     $          (dble(dsigmac).lt.0d0) .and.
     $          (errors(iter)*errors(iter-1).gt.0d0)) then
              if (errors(iter).lt.errors(iter-1)) then

                bracket = .true.
                if (values(iter).lt.values(iter-1)) then
                  lower = values(iter)
                  upper = values(iter-1)
                  reslower = errors(iter)
                  resupper = errors(iter-1)
                else
                  lower = values(iter-1)
                  upper = values(iter)
                  reslower = errors(iter-1)
                  resupper = errors(iter)
                endif

                eout = values(iter) - 
     $                 0.6180d0*(values(iter)-values(iter-1))
                values(iter) = values(iter-1)
                errors(iter) = errors(iter-1)
                skipped = .true.
              endif
            endif

            !bracket solution
            if (iter.gt.1 .and. (.not.skipped) ) then
              call gw_findbracket(bracket,values(iter),errors(iter),
     $                            values(iter-1),errors(iter-1),
     $                            lower,upper,reslower,resupper)
            endif

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  sum(rms3)/3d0.lt.pars%thresh .or.
     $                  (bracket .and. abs(upper-lower).lt.pars%thresh)

            ! Exit if converged
            call ga_brdcst(mt_log,lvars,3*logs,0)

            if (converged) then
              eout = ein
              exit
            endif

            if (.not.skipped) then
              call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,pars%thresh,
     $                       radius,iter)
            endif

#ifdef GWDEBUG
            if (main) write(luout,9901) iter,iqp,util_wallsec()-timer2,
     $                                  Ein,Eout
 9901       format('  ## Finished iter',I3,' of qp ',I3,' in ',F12.2,'s'
     $             ,/,'     Ein: ',F12.6,' Eout: ',F12.6)
            if (main) write(luout,9902) residual,eout-ein,lower,upper
 9902       format('    res ',F12.6,' step ',F12.6,' lower ',F12.6,
     $             ' upper ', F12.6)
            call util_flush(luout)
#endif
            !save information
            rms3(3) = rms3(2)
            rms3(2) = rms3(1)
            rms3(1) = abs(eout-ein)

          enddo

          ! Skip graphical solution 
          if (converged) then
            enew(llmo+iqp-1,isp) = eout
            fixed(iqp) = .true.
            if (bracket) then
              esterror(llmo+iqp-1) = min(upper-lower,abs(residual))
            else
              esterror(llmo+iqp-1) = abs(residual)
            endif

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
          nspace = sign(0.00037d0,eout)
          interval = maxval(eold(llmo+myll-1:llmo+myul-1,isp)) -
     $               minval(eold(llmo+myll-1:llmo+myul-1,isp)) + 0.004d0
          npoints = min(int(interval/abs(nspace))+1,100)
          nspace = interval/dble(npoints)

          omega(1) = eout - sign(0.00037d0,eout)
          do ipoints=1,npoints
            timer2 = util_wallsec()
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
                esterror(llmo+iqp-1) = abs(resmat(ipoints,iqp))
                cycle
              endif

              ! See if there is a fixed point solution
              if (ipoints.gt.1) then
                if(resmat(ipoints-1,iqp)*resmat(ipoints,iqp).lt.0d0)then
                  fixed(iqp) = .true.
                  der=(resmat(ipoints,iqp)-resmat(ipoints-1,iqp))/nspace
                  enew(llmo+iqp-1,isp) = omega(ipoints) - 
     $                                   resmat(ipoints,iqp)/der
                  esterror(llmo+iqp-1) = nspace/2d0
                  cycle
                endif
              endif

              done = .false.
            enddo
#ifdef GWDEBUG
            if (main) write(luout,9903) ipoints,util_wallsec()-timer2
 9903       format('  ## Finished point ',I3,' in ',F12.2,' s')
#endif
            if (done) exit
          enddo

 200      continue

          if (.not.main) cycle

          do iqp=llqp,ulqp
            if (fixed(iqp)) then
              write(luout,9020) llmo+iqp-1,
     $                      (enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                       esterror(llmo+iqp-1)*ha2ev
            else
              enew(llmo+iqp-1,isp) = omega( 
     $                minloc(abs(resmat(1:npoints,iqp)),dim=1) )
              warning = .true.
              write(luout,9030) llmo+iqp-1,
     $                        (enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                         nspace*ha2ev
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
     $       5x,' Error (eV)',/,10x,36('-'))
 9010 format(/,15x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,' Error (eV)',/,10x,26('-'))
 9020 format(10x,I3,3x,F12.3,4x,F12.3)
 9030 format(10x,I3,3x,F12.3,4x,F12.3,' *')

      end subroutine

