      subroutine gw_findclusters(evals,clusters,nclusters,nqp,window)
        implicit none
        integer,intent(in)  :: nqp
        integer,intent(out) :: nclusters,clusters(nqp)
        double precision,intent(in) :: evals(nqp),window

        integer iqp, nevals, ll, icluster, ul
        double precision :: targ, average, delta, stdev

        ! initialization
        nclusters   = 0
        clusters(:) = 0
        ll = 1

 100    continue

        nclusters = nclusters + 1
        targ = evals(ll) + window
        icluster = 1
        
        do iqp=ll+1,nqp
          if (evals(iqp)*evals(iqp-1).lt.0d0) exit
          if (evals(iqp).le.targ) then
            icluster = icluster + 1
          else
            exit
          endif
        enddo
        clusters(nclusters) = icluster
        ll = ll + icluster

        if (ll.le.nqp) goto 100

      end subroutine


      subroutine gw_diis(values,errors,old,new,ndiis,niter,nqp)
      implicit none
      integer,intent(inout) :: ndiis,nqp,niter
      double precision,intent(inout) :: errors(nqp,*),values(nqp,*)
      double precision,intent(out)  :: new(nqp)
      double precision,intent(in)  :: old(nqp)
      double precision,allocatable :: deltax(:,:),deltaf(:,:)
      double precision,allocatable :: work(:),f(:,:),b(:)
      double precision,allocatable :: x(:,:),step(:)
      integer, allocatable :: ipiv(:)

      integer idiis,jdiis,info,j,iqp

      double precision,parameter :: beta = 0.1d0
      double precision,parameter :: alpha = 0.3d0
      integer,parameter :: m = 15

 100  continue

      if (niter.eq.0) values(:,niter+1) = old(:)

      do iqp=1,nqp
        values(iqp,niter+2) = new(iqp)
        errors(iqp,niter+1) = new(iqp) - old(iqp)
      enddo

      ndiis = ndiis + 1
      if (niter.eq.0) ndiis = 0

      if (ndiis.le.1) then
        return
      endif

      allocate(f(ndiis+1,ndiis+1),b(ndiis+1),
     $         work(10*ndiis),ipiv(10*ndiis))

      do idiis=1,ndiis
        do jdiis=1,ndiis
          f(jdiis,idiis) = dot_product(errors(:,niter+1-jdiis+1),
     $                                 errors(:,niter+1-idiis+1))
        enddo
      enddo
      f(ndiis+1,:) = -1d0
      f(:,ndiis+1) = -1d0
      f(ndiis+1,ndiis+1) = 0d0

      b(1:ndiis) =  0d0
      b(ndiis+1) = -1d0

      call dsysv('u',ndiis+1,1,F,ndiis+1,ipiv,b,ndiis+1,work,10*ndiis,
     $            info)

      if (info.ne.0) then
        write(*,*) 'info 0'
        ndiis = ndiis - 1
        deallocate(f,b,work,ipiv)
        return
      endif

      new(1:nqp) = 0d0
      do idiis=1,ndiis
        new(1:nqp) = new(1:nqp) + b(idiis)*values(:,niter+1-idiis+1)
      enddo

      deallocate(f,b,work,ipiv)

      end subroutine

      subroutine gw_findbracket(found,enew,resnew,eold,resold,
     $                          xlower,xupper,flower,fupper)
      implicit none
      logical found
      double precision enew,resnew,eold,resold
      double precision xlower,xupper,flower,fupper

      if (.not.found) then
        
        if ( (resnew*resold.lt.0d0) ) then 
          found = .true.
          xlower = min(enew,eold)
          xupper = max(enew,eold)          
          if (enew.gt.eold) then
            fupper = resnew
            flower = resold
          else
            fupper = resold
            flower = resnew
          endif
        endif

      else

        if (fupper*flower.gt.0d0 .and. resnew*fupper.gt.0d0) then

          if (fupper.lt.0d0) then
            if (resnew.gt.fupper) then
              xupper = enew
              fupper = resnew
            else
              xlower = enew
              flower = resnew
            endif
          else
            if (resnew.lt.flower) then
              xlower = enew
              flower = resnew
            else
              xupper = enew
              fupper = resnew
            endif
          endif

        elseif (abs(fupper).lt.abs(flower)) then

          if (resnew*fupper.lt.0d0) then
            xlower = enew
            flower = resnew
          elseif (resnew*flower.lt.0d0) then
            xupper = enew
            fupper = resnew
          endif

        else

          if (resnew*flower.lt.0d0) then
            xupper = enew
            fupper = resnew
          elseif (resnew*fupper.lt.0d0) then
            xlower = enew
            flower = resnew
          endif
        endif

      endif


      end subroutine


      subroutine gw_cdgw_ecrpa(pi,ldpi,nri,ecrpa)
        implicit none
        integer,intent(in) :: ldpi,nri
        double precision,intent(in) :: pi(ldpi,nri)
        double precision,intent(out) :: ecrpa

        integer iri
        double precision det


        if ( mod(nri,2) .eq. 0 ) then
          det = 1d0
          do iri=1,nri/2
            det = det*(pi(iri,iri)*pi(iri+1,iri))**2
          enddo
        else
          det = pi(1,1)**2
          do iri=2,nri/2+1
            det = det*(pi(iri-1,iri)*pi(iri,iri))**2
          enddo
        endif

        ecrpa = log(det)


      end subroutine


      subroutine gw_analytic_ecrpa(rpa,eia,mynpoles,lpole,ecrpa)
      implicit none
      integer,intent(in) :: mynpoles,lpole
      double precision,intent(in) :: eia(mynpoles),rpa(mynpoles,*)
      double precision,intent(out) :: ecrpa

      integer ipole

      ecrpa = 0d0
!$omp parallel do reduction(+: ecrpa)
      do ipole=lpole,lpole+mynpoles-1
        ecrpa = ecrpa + eia(ipole-lpole+1) + 
     $    (rpa(ipole-lpole+1,ipole)/eia(ipole-lpole+1) -
     $    eia(ipole-lpole+1))/2d0
      enddo
!$omp end parallel do

      ecrpa = -0.5d0*ecrpa

      call ga_dgop(1,ecrpa,1,'+')

      end subroutine

c     -----------------------------------------------------------------

      subroutine gw_analytic_ecgm(ecgm,omegain,Eold,Omegam,wmn,eta,
     $                            npoles,nmo,plus,me)

      implicit none
#include "gw.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      integer,intent(in) :: wmn,npoles,nmo,plus,me
      double precision,intent(in) :: omegain,Eold(nmo)
      double precision,intent(in) :: eta,Omegam(npoles)
      double precision,intent(out) :: ecgm

      integer :: ilo, ihi, jlo, jhi, ipole, imo, ld, adrWMN
      double precision :: omega,w,denom,tecgm

      call ga_distribution(wmn, me, ilo, ihi, jlo , jhi)
      call ga_access(wmn, ilo, ihi ,jlo, jhi, adrWMN, ld)

      tecgm = 0d0

!$omp parallel private(omega,w,denom) reduction(+: tecgm)
!$omp do collapse(2)
      do imo=plus+1,nmo
        do ipole=jlo,jhi
          omega = omegain - Eold(imo) - omegam(ipole)
          w = dbl_mb(adrWMN + (ipole-jlo)*ld + imo - 1)
          denom = omega**2 + eta
          tecgm = tecgm + w*omega/denom
        enddo
      enddo
!$omp end  do      
!$omp end parallel
c
      call ga_dgop(1038,tecgm,1,'+')
      ecgm = ecgm + 2d0*tecgm

      call ga_release(wmn,ilo,ihi,jlo,jhi)

      end subroutine

c     -----------------------------------------------------------------

      subroutine gw_cdgw_ecgm(pi,vx,nri,ldpi,ecgm,ecrpa,diag)

c     Compute the Galitskii-Migdal correlation energy using the
c     information available in CDGW.

        implicit none
        logical diag
        integer nri,ldpi
        double precision pi(ldpi,*),vx(ldpi,*),ecgm,ecrpa
        integer info,iri,jri
        double precision evals(nri),temp
        double precision work(nri,nri),work2(nri,nri)

        ecgm = 0d0
        ecrpa = 0d0
        iri = (nri*(nri+1))/2/ldpi

        if (.not.diag) then
          call dpftri('n','l',nri,vx,info)
        endif

        if ( mod(nri,2) .eq. 0 ) then
          do iri=1,nri/2
            do jri=1,iri-1
              ecgm = ecgm + 2d0*pi(jri,iri)*vx(jri,iri)
            enddo
            do jri=iri+2,nri+1
              ecgm = ecgm + 2d0*pi(jri,iri)*vx(jri,iri)
            enddo
            ecgm = ecgm + (pi(iri,iri)-1d0)*(vx(iri,iri)-1d0)
            ecgm = ecgm + (pi(iri+1,iri)-1d0)*(vx(iri+1,iri)-1d0)
            ecrpa = ecrpa - pi(iri,iri) - pi(iri+1,iri)
          enddo
        else
          do iri=1,nri/2+1
            ecgm = ecgm+2d0*dot_product(pi(:,iri),vx(:,iri))
            ecgm = ecgm-pi(iri,iri)-vx(iri,iri)-pi(iri,iri)*vx(iri,iri)
            ecrpa = ecrpa - pi(iri,iri)
            if (iri.gt.1) then
              ecgm=ecgm-pi(iri-1,iri)-vx(iri-1,iri)-
     $                  pi(iri-1,iri)*vx(iri-1,iri)
              ecrpa = ecrpa - pi(iri-1,iri)
            endif
          enddo
          ecgm = ecgm + dble(nri)
        endif

        if (diag) then
          ecrpa = 0d0
        else
          ecrpa = ecrpa + dble(nri)
        endif

      end subroutine

c     -----------------------------------------------------------------

      subroutine gw_scissor(enew,eold,nmo,nocc,nvir,noqp,nvqp,ipol,me)

c     Applies a scissor operator to energies not corrected by GW.
c     Separate shifts for the occupied and virtual spaces.

      implicit none
      integer,intent(in) :: ipol, nmo, me
      integer,intent(in),dimension(ipol) :: nocc,nvir,noqp,nvqp
      double precision,intent(in),dimension(nmo,ipol) :: eold
      double precision,intent(out),dimension(nmo,ipol) :: enew

      integer isp,imo
      double precision shift

      do isp=1,ipol

        ! shift occupied states as HOMO
        if (noqp(isp).lt.nocc(isp) .and.noqp(isp).gt.0) then
          shift = 0d0
          do imo=nocc(isp),nocc(isp)-noqp(isp)+1,-1
            shift = shift + enew(nocc(isp),isp) - eold(nocc(isp),isp)
          enddo
          shift = shift/noqp(isp)
          if (me.eq.0) write(*,9000) shift*27.2114d0
          do imo=nocc(isp)-noqp(isp),1,-1
            enew(imo,isp) = eold(imo,isp) + shift
          enddo
        endif

        ! shift unoccupied states as LUMO
        if (nvqp(isp).lt.nvir(isp) .and. nvqp(isp).gt.0) then
          shift = 0d0
          do imo=nocc(isp)+1,nocc(isp)+nvqp(isp)
            shift = shift + enew(nocc(isp)+1,isp)-eold(nocc(isp)+1,isp)
          enddo
          shift = shift/nvqp(isp)
          if (me.eq.0) write(*,9010) shift*27.2114d0
          do imo=nocc(isp)+nvqp(isp)+1,nmo
            enew(imo,isp) = eold(imo,isp) + shift
          enddo
        endif
      enddo

 9000 format(/,3X,'Applying ',F8.4,' eV scissor-shift to rest of ',
     $         'particle states',/)
 9010 format(/,3X,'Applying ',F8.4,' eV scissor-shift to rest of ',
     $         'hole states',/)
      end subroutine

c     *****************************************************************

      logical function is_near(test,targ,tol)
      implicit none
      double precision, intent(in) :: test, targ, tol

      if (abs(test-targ).le.tol) then
        is_near = .true.
      else
        is_near = .false.
      endif

      return
      end function

c     *****************************************************************

      subroutine gw_sigmax(pars,sigmax,nri,nmo)
      implicit none
#include "mafdecls.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nri, nmo
      double precision,intent(out),dimension(nmo,*) :: sigmax

      double precision :: temp
      double precision,external :: ddot

      integer idx, idum, isp
      integer imo, jmo, kmo, ipol
      integer addrOO, addrOV, addrVV
      integer,external :: omp_get_num_threads,omp_get_thread_num

      pars%ex = 0d0

      do ipol=1,pars%ipol
        if (pars%nqp(ipol).eq.0) cycle
        if (pars%oolo(ipol).gt.pars%oohi(ipol)) then
          addrOO = 1
        else
          call ga_access(pars%g_erioo(ipol),1,nri,pars%oolo(ipol),
     $                   pars%oohi(ipol),addrOO,idum)
        endif
        if (pars%ovlo(ipol).gt.pars%ovhi(ipol)) then
          addrOV = 1
        else
          call ga_access(pars%g_eriov(ipol),1,nri,pars%ovlo(ipol),
     $                   pars%ovhi(ipol),addrOV,idum)
        endif


        do imo=min(1,pars%llmo(ipol)),
     $         max(pars%llmo(ipol)+pars%nqp(ipol)-1,pars%nocc(ipol))
          jmo = imo - pars%llmo(ipol) + 1
          temp = 0d0

          if (imo.le.pars%nocc(ipol)) then
c
            do kmo=1,pars%nocc(ipol)
              idx = (imo-1)*pars%nocc(ipol)+kmo-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOO
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
c
          else
c
            do kmo=1,pars%nocc(ipol)
              idx = (kmo-1)*pars%nvir(ipol)+imo-pars%nocc(ipol)-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOV
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
c
          endif
          sigmax(imo,ipol) = temp

        enddo

      enddo

      call ga_dgop(1038,sigmax(1,1),nmo*pars%ipol,'+')

      do ipol=1,pars%ipol
        pars%ex = pars%ex+sum(sigmax(1:pars%nocc(ipol),ipol))/pars%ipol
      enddo

      end subroutine

c     *****************************************************************

      double precision function gw_efermi(eig,nocc)
      implicit none
      integer,intent(in) :: nocc
      double precision,intent(in),dimension(*) :: eig
      gw_efermi = (eig(nocc) + eig(nocc+1))/2d0
      return
      end function

c     *****************************************************************

      subroutine gw_shift(eig,shifted,efermi,nbf,nmo,npol)
      implicit none
      integer,intent(in) :: nbf, nmo, npol
      double precision,intent(in) :: efermi(npol),eig(nbf,npol)
      double precision,intent(out) :: shifted(nmo,npol)
      integer imo,ipol
!$omp parallel do collapse(2)      
      do ipol=1,npol
        do imo=1,nmo
          shifted(imo,ipol) = eig(imo,ipol) - efermi(ipol)
        enddo
      enddo
!$omp end parallel do      
      return
      end subroutine

c     *****************************************************************

      subroutine gw_get_eia(eig,eia,nocc,nvir)
      implicit none
#include "global.fh"      
      integer,intent(in) :: nocc,nvir
      double precision,intent(in),dimension(nocc+nvir) :: eig
      double precision,intent(out),dimension(*) :: eia
      integer ipole,imo,jmo,me,nprocs,jpole

      me = ga_nodeid()
      nprocs = ga_nnodes()
c
      ipole = 0
      jpole = 0
      do imo=1,nocc
        do jmo=nocc+1,nocc+nvir
          jpole = jpole + 1
          if (me.ne.mod(jpole-1,nprocs)) cycle
          ipole = ipole + 1
          eia(ipole) = eig(jmo) - eig(imo)
        enddo
      enddo
c      
      return
      end subroutine

c     *****************************************************************

      subroutine gw_get_iwidx(idx,nmo,nocc,nvir,mynmo,nqp,llmo)
C
C     Array of MO pairs whose integrals are stored in the
C     current MPI rank. 
C
C     This info will be used to compute and store matrix elements of
C     the screened Coulomb matrix for the same MO pairs.
C
      implicit none
      integer,intent(in) :: nmo, nocc, nvir, mynmo, nqp, llmo
      integer,intent(out),dimension(0:mynmo,*) :: idx

      integer me, imo, jmo, kmo, lmo, nprocs
      integer,external :: ga_nodeid, ga_nnodes

      me = ga_nodeid()
      nprocs = ga_nnodes()

      do imo=llmo,llmo+nqp-1
        jmo = imo - llmo + 1
        lmo = 0

        if (imo.le.nocc) then
          do kmo=1,nocc
            if (me.ne.mod((imo-1)*nocc+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        else
          do kmo=1,nocc
            if (me.ne.mod((kmo-1)*nvir+imo-nocc-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-nocc-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        endif
        idx(0,jmo) = lmo

      enddo

      return
      end subroutine

c     *****************************************************************

      subroutine gw_symmat(matrix,n,lo,hi)
      implicit none
      integer,intent(in) :: n, lo, hi
      double precision,intent(inout),dimension(n,n) :: matrix
      integer i,ilo

      do i=1,n-1
        matrix(i,i+1:) = matrix(i+1:,i)
      enddo

      return
      end subroutine

c     *****************************************************************

      subroutine gw_rfpmm(rfp,vecsb,vecsc,ld,n,nvec,packed)

      ! Matrix-matrix multiply using one rectangular full-packed matrix
      !
      ! Note: Here we are using the right RFP format, not the shifted
      ! one used in other places. LD must have the right dimensions: 
      ! LD = n+(1-mod(n/2))
      !

      implicit none
      logical,intent(in) :: packed
      integer,intent(in) :: ld, n, nvec
      double precision,intent(in) :: rfp(ld,*),vecsb(n,nvec)
      double precision,intent(out) :: vecsc(n,nvec)


      if (packed .and. mod(n,2).eq.0) then
        call dsymm('l','u',n/2,nvec,1d0,rfp,n+1,vecsb(n/2+1,1),n,0d0,
     $              vecsc(n/2+1,1),n)
        call dsymm('l','l',n/2,nvec,1d0,rfp(2,1),n+1,vecsb,n,0d0,
     $              vecsc,n)
        call dgemm('n','n',n/2,nvec,n/2,1d0,rfp(n/2+2,1),n+1,vecsb,n,
     $              1d0,vecsc(n/2+1,1),n)
        call dgemm('t','n',n/2,nvec,n/2,1d0,rfp(n/2+2,1),n+1,
     $              vecsb(n/2+1,1),n,1d0,vecsc,n)
      else if (packed) then
        call dsymm('l','l',n/2+1,nvec,1d0,rfp,n,vecsb,n,0d0,
     $              vecsc,n)
        call dsymm('l','u',n/2,nvec,1d0,rfp(1,2),n,vecsb(n/2+2,1),n,0d0,
     $              vecsc(n/2+2,1),n)
        call dgemm('n','n',n/2,nvec,n/2+1,1d0,rfp(n/2+2,1),n,vecsb,n,
     $              1d0,vecsc(n/2+2,1),n)
        call dgemm('t','n',n/2+1,nvec,n/2,1d0,rfp(n/2+2,1),n,
     $              vecsb(n/2+2,1),n,1d0,vecsc,n)
      else
        call dgemm('n','n',n,nvec,n,1d0,rfp,n,vecsb,n,0d0,vecsc,n)
      endif

      end subroutine

c     *****************************************************************

      subroutine gw_newton_debug(residual,eout,ein,lower,upper,efermi,
     $                           time,bracket,iqp,iter)
      implicit none

#include "gw.fh"
#include "stdio.fh"

      integer iqp,iter
      logical bracket
      double precision residual, eout, ein, lower, upper, efermi, time

      write(luout,9901) iter,iqp,time,
     $                  ha2ev*(Ein+efermi),
     $                  ha2ev*(Eout+efermi)
      if (.not.bracket) then
        write(luout,9902) ha2ev*residual,ha2ev*(eout-ein)
      else
        write(luout,9903) ha2ev*residual,ha2ev*(eout-ein),
     $                    ha2ev*(lower+efermi),
     $                    ha2ev*(upper+efermi)
      endif

      call util_flush(luout)

 9901 format('  ## Finished iter',I3,' of qp ',I3,' in ',F12.2,'s'
     $             ,/,'     Ein: ',F12.6,' Eout: ',F12.6)
 9902 format('    res ',F12.6,' step ',F12.6)
 9903 format('    res ',F12.6,' step ',F12.6,' lower ',F12.6,
     $             ' upper ', F12.6)
      end subroutine

c     *****************************************************************

      subroutine gw_findfixed(npoints, resmat, omegamat, sigmamat,
     $             sigma, enew, sigma_x, nspace, fixed, thresh)

      implicit none

      logical fixed
      integer ipoints, npoints
      double precision sigma_x, nspace, thresh
      double precision sigma, enew 
      double precision omegamat(npoints),sigmamat(npoints)
      double precision resmat(npoints)

      do ipoints=1,npoints

        ! If residual is already small, the solution has been found
        if (abs(resmat(ipoints)).lt.thresh) then
          enew = omegamat(ipoints)
          sigma = sigmamat(ipoints) + sigma_x
          fixed = .true.
          exit
        endif

        ! if residual just changed sign, bisect solution
        if (ipoints.gt.1) then
          if (resmat(ipoints-1)*resmat(ipoints).lt.0d0) then
            enew = omegamat(ipoints) - nspace/2d0
            sigma = sigma_x +   
     $         (sigmamat(ipoints)+sigmamat(ipoints-1))/2d0
            fixed = .true.
            exit
          endif
        endif

      enddo

      end subroutine

c     *****************************************************************

      subroutine gw_print_energy(llmo,llqp,ulqp,enew,efermi,error,fixed,
     $                           warning)
      implicit none
#include "gw.fh"
#include "stdio.fh"
      integer llmo, llqp, ulqp
      logical fixed(*), warning
      double precision enew(*), efermi, error(*)

      integer iqp, jqp

      do iqp=llqp,ulqp
        jqp = iqp - llqp + 1
        if (fixed(jqp)) then
          write(luout,9020) llmo+iqp-1,
     $                     (enew(jqp)+efermi)*ha2ev,error(jqp)*ha2ev
        else
          warning = .true.
          write(luout,9030) llmo+iqp-1,
     $                      (enew(jqp)+efermi)*ha2ev,error(jqp)*ha2ev
        endif
      enddo

 9020 format('@GW',7x,I3,3x,F12.3,4x,F12.3)
 9030 format('@GW',7x,I3,3x,F12.3,4x,F12.3,' *')

      end subroutine
