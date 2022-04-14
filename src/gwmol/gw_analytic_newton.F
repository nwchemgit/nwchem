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
      double precision,intent(in) :: Omega(maxpoles*pars%ipol)
      double precision,intent(inout) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(out) :: Enew(pars%nmo,pars%ipol)

      logical,target  :: lvars(2)
      logical,pointer :: converged, bracket
      logical main, warning, fixed(pars%nmo), ok
      integer isp, totalqp, llmo, iqp, iter, g_wmn, plus, imo
      integer l_sw,k_sw,dbls,logs
      integer npoles

      double precision timer, timer2, thresh
      double precision Ein, Eout, SigmaOld, Sigma_X, SigmaC
      double precision d2SigmaC, dSigmaC
      double precision residual, dresidual, d2residual
      double precision esterror(pars%nmo)

      integer nclusters, clusters(pars%nmo), myll, myul
      integer llqp, ulqp, icluster
      double precision upper, lower, resupper, reslower, const
      double precision values(20), errors(20), radius, maxradius
c
      !
      ! Initialize variables
      !
      pars%ecgm = 0d0
      main = pars%me .eq. 0
      dbls = ma_sizeof(mt_dbl,1,mt_byte)
      logs = ma_sizeof(mt_log,1,mt_byte)

c     !
      ! The total number of poles 
      !
      npoles = sum(pars%npoles(1:pars%ipol))

      !
      ! Tighten threshold for evGW0 and evGW iterations
      !
      thresh = pars%thresh
      if (eviter.gt.2) thresh = thresh/10d0

      converged => lvars(1)
      bracket   => lvars(2)

      if (eviter.gt.4) then 
        maxradius = 0.010d0
      else
        maxradius = 0.1d0
      endif

c     !
      ! Screened-Coulomb matrix elements array
      !
      if (.not.ga_create(mt_dbl,pars%nmo,npoles,'wmn',pars%nmo,0,g_wmn))
     $  call errquit('can''t create Wmn array',0,GA_ERR)

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

        !
        ! Find eigenvalues tightly clustered. The solution will be
        ! attempted starting from the value closer to the Fermi energy
        ! and using its result as guess for the subsequent orbitals.
        !
        call gw_findclusters(Eold(llmo,isp),clusters,nclusters,totalqp,
     $                       0.05d0)

        !
        ! Make sure all ranks have the same information
        !
        call ga_brdcst(1038,nclusters,ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(1038,clusters,
     $                 ma_sizeof(mt_int,nclusters,mt_byte),0)


        do imo=1,pars%nmo
          if (Eold(imo,isp).gt.0d0) then
            plus = imo - 1
            exit
          endif
        enddo


        !
        ! Loop over all clusters of eigenvalues
        !
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

          if (eviter.lt.10) then
            if (myul.lt.ulqp) then
              Eout = Enew(llmo+iqp,isp) + 0.002d0
            elseif (myll.gt.llqp) then
              Eout = Enew(llmo+iqp-2,isp) - 0.002d0
            endif
          endif

          bracket = .false.
          const = Eold(llmo+iqp-1,isp) + Sigma_X - sigmaold
          lower = 1d100
          upper = -1d100

c         Compute Screened-Coulomb matrix elements
          timer = util_wallsec()          
          call gw_analytic_wmn(pars,g_wmn,iqp+llmo-1,isp,nri)
          pars%timer%wm = pars%timer%wm - timer + util_wallsec()


c         Compute Galiitski-Migdal correlation energy
          if (llmo+iqp-1.le.pars%nocc(isp) ) then
            call gw_analytic_ecgm(pars%ecgm,Eold(llmo+iqp-1,isp),
     $       Eold(1,isp),omega,g_wmn,pars%eta**2,npoles,pars%nmo,plus,
     $       pars%me)
          endif

c        
          do iter=1,pars%qpiter+5
            timer2 = util_wallsec()

            Ein = Eout

            call ga_sync()
            timer = util_wallsec()

c           Comput self-energy            
            call gw_analytic_sigmac(sigmaC,dsigmaC,d2sigmaC,Ein,
     $           Eold(1,isp),omega,g_wmn,pars%eta**2,npoles,pars%nmo,
     $           plus,pars%me)

            pars%timer%sigmac = pars%timer%sigmac - timer +  
     $                          util_wallsec()


            residual = const + SigmaC - Ein
            dresidual = dSigmaC - 1d0
            d2residual = d2SigmaC

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

          if (main) call gw_print_energy(llmo,llqp,ulqp,
     $                    enew(llmo+llqp-1,isp),efermi(isp),
     $                    esterror(llmo+llqp-1),fixed(llqp),
     $                    warning)

        enddo
        if (warning.and.main) then
          write(luout,*)
          write(luout,*) ' * Result did not converge'
          write(luout,*)
        endif
      enddo

      if (.not.ga_destroy(g_wmn))
     $  call errquit('can''t destroy Wmn',0,GA_ERR)   
      if (.not.pars%evgw0) then
        if (.not.ga_destroy(pars%g_erim))
     $    call errquit('can''t destroy transformed ERIs',0,GA_ERR)   
      endif

 9000 format(/,18x,'Alpha Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))
 9010 format(/,18x,' Beta Orbitals',/,10x,'State',5x,'Energy (eV)',
     $       5x,'Error (eV)',/,10x,35('-'))

      end subroutine
