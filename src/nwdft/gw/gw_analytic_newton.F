      subroutine gw_analytic_newton(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                          Efermi,nri,maxpoles,EViter)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"      
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri,maxpoles,EViter
      double precision,intent(in) :: Efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles,pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical,target  :: lvars(3)
      logical,pointer :: converged, bracket, skipped
      logical main, warning, fixed(pars%nmo)
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw,dbls,logs

      double precision timer, timer2, thresh
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC
      double precision d2SigmaC, dSigmaC
      double precision residual, dresidual, d2residual
      double precision esterror(pars%nmo)

      integer nclusters, clusters(pars%nmo), myll, myul
      integer llqp, ulqp, icluster
      double precision upper, lower, resupper, reslower, const
      double precision values(20), errors(20), radius, maxradius

      pars%ecgm = 0d0
      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

      thresh = pars%thresh
      if (eviter.gt.2) thresh = thresh/10d0

      converged => lvars(1)
      bracket   => lvars(2)
      skipped   => lvars(3)

      if (eviter.gt.4) then 
        maxradius = 0.010d0
      else
        maxradius = 0.1d0
      endif

      do isp=1,pars%ipol
        warning = .false.

        totalqp = pars%noqp(isp) + pars%nvqp(isp)
        if (totalqp.lt.1) cycle

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

        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.05d0)
        call ga_brdcst(mt_int,nclusters,ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(mt_int,clusters,
     $                 ma_sizeof(mt_int,nclusters,mt_byte),0)

        if (.not.ga_create(mt_dbl,pars%nmo,pars%npoles(isp),'wmn',
     $                     pars%nmo,0,g_wmn))
     $    call errquit('can''t create array',0,GA_ERR)

        do imo=1,pars%nmo
          if (Eold(imo,isp).gt.0d0) then
            plus = imo - 1
            exit
          endif
        enddo


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
              Eout = Enew(llmo+iqp,isp)
            elseif (myll.gt.llqp) then
              Eout = Enew(llmo+iqp-2,isp)
            endif
          endif

          bracket = .false.
          const = Eold(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = 0d0
          upper = 0d0

          timer = util_wallsec()
          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)
          pars%timer%wm = pars%timer%wm - timer + util_wallsec()

          if (llmo+iqp-1.le.pars%nocc(isp)) then
            call gw_analytic_ecgm(pars%ecgm,Eold(llmo+iqp-1,isp),Eold,
     $                          omega(1,isp),g_wmn,pars%eta**2,
     $                          pars%npoles(isp),pars%nmo,plus,pars%me)
          endif
c        
          do iter=1,pars%qpiter+5
            timer2 = util_wallsec()

            Ein = Eout
            skipped = .false.

            call ga_sync()
            timer = util_wallsec()
            call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,Eold,
     $                              omega(1,isp),g_wmn,
     $                              pars%eta**2,pars%npoles(isp),
     $                              pars%nmo,plus,pars%me)
            pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC

            values(iter) = Ein
            errors(iter) = residual
            radius = min(maxradius,0.1d0*maxradius/residual**2)


            ! Find out if we might have skipped a solution
#ifdef GW_DOSKIP
            if ((iter.gt.1) .and. (.not.bracket) .and. 
     $          (dble(dsigmac).lt.0d0) .and.
     $          (errors(iter)*errors(iter-1).gt.0d0)) then
              if (abs(errors(iter)).gt.abs(errors(iter-1))) then

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
     $                 0.70d0*(values(iter)-values(iter-1))
                values(iter) = values(iter-1)
                errors(iter) = errors(iter-1)
                skipped = .true.
              endif
            endif
#endif

            ! Find out if we have bracketed the solution
            if (iter.gt.1 .and. (.not.skipped) ) then
              call gw_findbracket(bracket,values(iter),errors(iter),
     $                            values(iter-1),errors(iter-1),
     $                            lower,upper,reslower,resupper)
            endif

            ! Check convergence of QP equation
            converged = abs(residual).lt.thresh .or. 
     $                  (bracket .and. abs(upper-lower).lt.thresh)

            ! Exit if converged
            call ga_brdcst(mt_log,lvars,3*logs,0)
            if (converged) then
              Eout = Ein
#ifdef GWDEBUG
             if (main) write(luout,9901) iter,llmo+iqp-1,
     $                                  util_wallsec()-timer2,
     $                                  ha2ev*(Ein+efermi(isp)),
     $                                  ha2ev*(Eout+efermi(isp))
              if (.not.bracket.and.main) then
                write(luout,9902) ha2ev*residual,
     $                           ha2ev*(eout-ein)
              elseif (main) then
                write(luout,9903) ha2ev*residual,
     $                           ha2ev*(eout-ein),
     $                           ha2ev*(lower+efermi(isp)),
     $                           ha2ev*(upper+efermi(isp))
              endif
#endif
              exit
            endif

            if (.not.skipped) then
              call gw_donewton(ein,eout,residual,dresidual,lower,upper,
     $                       reslower,resupper,bracket,thresh,
     $                       radius,iter)
            endif

#ifdef GWDEBUG
            if (main) write(luout,9901) iter,llmo+iqp-1,
     $                                  util_wallsec()-timer2,
     $                                  ha2ev*(Ein+efermi(isp)),
     $                                  ha2ev*(Eout+efermi(isp))
 9901       format('  ## Finished iter',I3,' of qp ',I3,' in ',F12.2,'s'
     $             ,/,'     Ein: ',F12.6,' Eout: ',F12.6)
            if (.not.bracket .and. main) then
              write(luout,9902) ha2ev*residual,
     $                          ha2ev*(eout-ein)
            elseif (main) then
              write(luout,9903) ha2ev*residual,
     $                          ha2ev*(eout-ein),
     $                          ha2ev*(lower+efermi(isp)),
     $                          ha2ev*(upper+efermi(isp))
            endif

 9902       format('    res ',F12.6,' step ',F12.6)
 9903       format('    res ',F12.6,' step ',F12.6,' lower ',F12.6,
     $             ' upper ', F12.6)
            call util_flush(luout)
#endif

          enddo

          Enew(llmo+iqp-1,isp)  = Eout
          Sigma(llmo+iqp-1,isp) = SigmaC + Sigma_X

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

          if (main) then
            do iqp=llqp,ulqp
              if (fixed(iqp)) then
                write(luout,9020) llmo+iqp-1,
     $                       (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                       esterror(llmo+iqp-1)*ha2ev
              else
                warning = .true.
                write(luout,9030) llmo+iqp-1,
     $                       (Enew(llmo+iqp-1,isp)+efermi(isp))*ha2ev,
     $                       esterror(llmo+iqp-1)*ha2ev
              endif
            enddo
          endif

        enddo
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Result did not converge'
          write(luout,*)
        endif
        if (.not.(ga_destroy(g_wmn)))
     $    call errquit('can''t destroy Wmn',0,GA_ERR)   
        if ( .not.pars%evgw0 ) then
          if (.not.(ga_destroy(pars%g_erim(isp))))
     $      call errquit('can''t destroy ERIm',0,GA_ERR)   
        endif
      enddo


 9000 format(/,18x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))
 9010 format(/,18x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))
 9020 format(10x,I3,3x,F12.3,4x,F12.3)
 9030 format(10x,I3,3x,F12.3,4x,F12.3,' *')

      end subroutine
