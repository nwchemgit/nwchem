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
      double precision,intent(in) :: iW(mynmo,pars%ngl,maxqp,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical main, converged, bracket, warning, fixed(maxqp)
      integer isp, totalqp, llmo, iqp, iter, logs, dbls
      double precision timer, rms3(3), m
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC, Z
      double precision d2SigmaC, dSigmaC, Linear
      double precision residual, residualOld, dresidual, d2residual,Etmp
      double precision quadratic, dquadratic, d2quadratic, dquadraticold
      double precision quadraticstep, linearstep, step, dsigmacold
      double precision I,dI,d2I,R,dR,d2R,disgmacold

      integer idiis, nclusters, clusters(pars%nmo),myll,myul
      integer llqp,ulqp,icluster
      double precision upper,lower,resupper,reslower,const


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
            Eout = Enew(llmo+iqp,isp)
          elseif (myll.gt.llqp) then
            Eout = Enew(llmo+iqp-2,isp) 
          elseif (eviter.eq.0) then
            Eout = Eout + 0.10d0*(sigma_x - sigmaold)
          else
            Eout = Eout + sign(1.2d0*pars%thresh,Eout)
          endif


          call ga_zero(pars%g_sols)
          bracket = .false.
          rms3 = 1d100
          residualOld = 0d0
          dsigmacold = 0d0
          const = Eref(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = 0d0
          upper = 0d0
c        
          do iter=1,pars%qpiter+5

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

            if (iter.eq.1) then
              call gw_buildrn(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            else
              call gw_buildrn_minres(pars,Eold(1,isp),Eia,Ein,R,dR,
     $                      llmo+iqp-1,maxpoles,nri,pars%me,pars%nmo,
     $                      pars%ipol,isp)
            endif

            pars%timer%rn = pars%timer%rn + util_wallsec() - timer

            SigmaC   = I + R
            dSigmaC  = dI + dR

            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0

            ! Find out if we have bracketed the solution
            if (iter.gt.1) then
              call gw_findbracket(bracket,ein,residual,etmp,residualold,
     $                            lower,upper,reslower,resupper)
            endif

 1234       format(5(F15.8))

            ! Check convergence of QP equation
            converged = abs(residual).lt.pars%thresh .or. 
     $                  sum(rms3)/3d0.lt.pars%thresh .or.
     $                  (bracket .and. abs(upper-lower).lt.pars%thresh)


            ! Exit if converged
            call ga_brdcst(mt_log,converged,logs,0)
            if (converged) then
              Eout = Ein
              exit
            endif

            ! Take next newton step
            call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,pars%thresh)


            ! Save information
            etmp = ein
            residualold = residual
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




      subroutine gw_findclusters(evals,clusters,nclusters,nqp,window)
        implicit none
        integer,intent(in)  :: nqp
        integer,intent(out) :: nclusters,clusters(nqp)
        double precision,intent(in) :: evals(nqp),window

        integer iqp, nevals, ll, icluster, ul
        double precision :: targ, average, delta, stdev

        ! initialization
        clusters(:) = 0
        nclusters = 1
        average = evals(1)
        ll = 1

 100    continue

        targ = average
        average = 0d0
        nevals = clusters(nclusters)
        icluster = 0
        do iqp=ll,nqp
          delta = evals(iqp) - targ
          if (iqp.gt.ll) then
            if (evals(iqp)*evals(iqp-1).lt.0d0) exit
          endif
          if (abs(delta).gt.window) then
            if (iqp.eq.ll) goto 200
            exit
          endif
          icluster = icluster + 1
          average = average + evals(iqp)
        enddo

        average = average/dble(icluster)
        clusters(nclusters) = icluster

 200    continue       

        if (nevals.eq.clusters(nclusters)) then
          ll = ll + clusters(nclusters)
          if (ll.gt.nqp) goto 300
          nclusters = nclusters + 1
          average = evals(ll) + window - 0.001d0
        endif
        
        goto 100

 300    continue

        ul = 0
        do icluster=1,nclusters
          ll = ul + 1
          ul = ul + clusters(icluster)
          average=sum(evals(ll:ul))/dble(clusters(icluster))
          stdev = maxval(evals(ll:ul))-minval(evals(ll:ul))
        enddo

      end subroutine


