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

      logical,target  :: lvars(2)
      logical,pointer :: converged,bracket
      logical main, warning, fixed(pars%nmo), done
      integer llqp, ulqp, ipoints, icluster, qplimit, npoints
      integer isp, totalqp, llmo, iqp, iter, logs, dbls
      double precision timer, timer2, rms3(3), const, nspace
      double precision Ein, Eout, SigmaOld, Sigma_X
      double precision residual, dresidual, der
      double precision values(20),errors(20),esterror(pars%nmo)
      double precision radius,maxradius
      GWTYP :: I,dI,R,dR,SigmaC,dSigmaC
      GWTYP :: resmat(500,pars%nmo)

      integer nclusters, clusters(pars%nmo),myll,myul
      double precision upper,lower,resupper,reslower,interval
      double precision omega(500), output(2,pars%nmo)


      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

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

        if ((isp.eq.1).and.main) then
          write(luout,9000)
        elseif(main) then
          write(luout,9010)
        endif

        llmo = pars%llmo(isp)

        !
        ! Find eigenvalues tightly clustered. The solution will be
        ! attempted starting from the value closer to the Fermi energy
        ! and using its result as guess for the subsequent orbitals.
        !
        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.020d0)

        !
        ! Loop over all clusters of eigenvalues
        ! 
        ulqp = 0
        do icluster=1,nclusters

          llqp = ulqp + 1
          ulqp = ulqp + clusters(icluster)

          myul = ulqp
          myll = llqp
          fixed(myll:myul) = .false.


 100      continue

          !
          ! First, we will do a few Newton steps in order to get an
          ! approximated solution. This will define the region were the
          ! omega grid is laid
          !
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

          !initialize variables
          call ga_zero(pars%g_sols)
          bracket = .false.
          converged = .false.
          upper = 1d100
          lower = -1d100
          rms3(:) = 1d100

          !fixed contribution
          const = Eref(llmo+iqp-1,isp) + Sigma_X -  SigmaOld

          ! we only need a guess
          do iter=1,8
            timer2 = util_wallsec()

            Ein = Eout

            call ga_sync()

            ! Contribution from the imaginary axis
            timer = util_wallsec()
            call gw_cdgw_buildin(iW(1,1,iqp,isp),Eold(1,isp),
     $                      dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),
     $                      iWidx(0,iqp,isp),Ein,pars%eta,I,dI,
     $                      pars%ngl,mynmo)
            call ga_sync()
            pars%timer%in = pars%timer%in + util_wallsec() - timer

            ! Contribution from the poles enclosed by the contour
            timer = util_wallsec()
            call gw_cdgw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
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

            radius = min(maxradius,0.1d0*maxradius/residual**2)

            !bracket solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,values(iter),errors(iter),
     $                            values(iter-1),errors(iter-1),
     $                            lower,upper,reslower,resupper)
            endif

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  sum(rms3)/3d0.lt.pars%thresh .or.
     $                  (bracket .and. abs(upper-lower).lt.pars%thresh)

            ! Exit if converged
            call ga_brdcst(1038,lvars,2*logs,0)

            if (converged) then
              eout = ein
              exit
            endif

            if (main) then
              call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,pars%thresh,
     $                       radius,iter)

#ifdef GWDEBUG
              call gw_newton_debug(residual,eout,ein,lower,
     $                  upper,efermi(isp),util_wallsec()-timer2,
     $                  bracket,llmo+iqp-1,iter)
#endif
            endif
            call ga_brdcst(1038,eout,ma_sizeof(mt_dbl,1,mt_byte),0)

            !save information
            rms3(3) = rms3(2)
            rms3(2) = rms3(1)
            rms3(1) = abs(eout-ein)

          enddo

          !
          ! If all Quasi-particle equations already converged in the
          ! Newton steps, then we need to do nothing
          !
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
          nspace = sign(0.0001d0,eout)
          interval = maxval(eold(llmo+myll-1:llmo+myul-1,isp)) -
     $               minval(eold(llmo+myll-1:llmo+myul-1,isp)) + 0.01d0
          npoints = min(int(interval/abs(nspace))+1,500)
          nspace = interval/dble(npoints)

          !
          ! Loop over omega grid
          !
          omega(1) = eout - sign(0.0001d0,eout)
          do ipoints=1,npoints
            timer2 = util_wallsec()
            if (ipoints.gt.1)
     $         omega(ipoints) = omega(ipoints-1) + sign(nspace,eout)


            timer = util_wallsec()
            call gw_cdgw_buildin_graph(iW(1,1,myll,isp),Eold(1,isp),
     $        dbl_mb(pars%k_glx),dbl_mb(pars%k_glw),iwidx(0,myll,isp),
     $        omega(ipoints),output,pars%eta,pars%ngl,mynmo,myul-myll+1)
            resmat(ipoints,myll:myul) = output(1,1:myul-myll+1)
            call ga_sync()
            pars%timer%in = pars%timer%in + util_wallsec() - timer


            timer = util_wallsec()
            call gw_cdgw_buildrn_graph(pars,Eold(1,isp),Eia,
     $        omega(ipoints),output,llmo+myll-1,maxpoles,nri,pars%me,
     $        pars%nmo,pars%ipol,isp,myul-myll+1)
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
              endif
              call ga_brdcst(1038,fixed(iqp),ma_sizeof(mt_log,1,
     $                       mt_byte),0)

              ! See if there is a fixed point solution
              if (ipoints.gt.1 .and. .not.fixed(iqp) ) then
                if(resmat(ipoints-1,iqp)*resmat(ipoints,iqp).lt.0d0)then
                  fixed(iqp) = .true.
                  der=(resmat(ipoints,iqp)-resmat(ipoints-1,iqp))/nspace
                  enew(llmo+iqp-1,isp) = omega(ipoints) - 
     $                                   resmat(ipoints,iqp)/der
                  esterror(llmo+iqp-1) = nspace/2d0
                endif
                call ga_brdcst(1038,fixed(iqp),ma_sizeof(mt_log,1,
     $                         mt_byte),0)
              endif

              if (fixed(iqp)) then
                call ga_brdcst(1038,enew(llmo+iqp-1,isp),
     $                       ma_sizeof(mt_dbl,1,mt_byte),0)
                cycle
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
            if (.not.fixed(iqp)) then
              enew(llmo+iqp-1,isp) = omega( minloc(abs(
     $               resmat(1:npoints,iqp) ),dim=1))
              esterror(llmo+iqp-1) = min(minval(abs(
     $               resmat(1:npoints,iqp))),2*nspace)
            endif
          enddo

          call gw_print_energy(llmo,llqp,ulqp,enew(llmo+llqp-1,isp),
     $           efermi(isp),esterror(llmo+llqp-1),fixed(llqp),warning)

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
     $       5x,' Error (eV)',/,10x,36('-'))

      end subroutine

