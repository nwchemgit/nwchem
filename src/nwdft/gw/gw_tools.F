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
      subroutine gw_findbracket(found,enew,resnew,eold,resold,
     $                          xlower,xupper,flower,fupper)
      implicit none
      logical found
      double precision enew,resnew,eold,resold
      double precision xlower,xupper,flower,fupper

      if (.not.found) then
        
        if (resnew*resold.lt.0d0) then
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

        !update intervals
        if (resnew*fupper.lt.0d0) then
          xlower = enew
          flower = resnew
        elseif (resnew*flower.lt.0d0) then
          xupper = enew
          fupper = resnew
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

      call dgop(1,ecrpa,1,'+')

      end subroutine


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
#ifdef USE_OPENMP
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif      

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
      call ga_dgop(mt_dbl,tecgm,1,'+')
      ecgm = ecgm + tecgm

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      

      end subroutine


      subroutine gw_cdgw_ecgm(pi,vx,nri,ldpi,ecgm,ecrpa)
        implicit none
        integer nri,ldpi
        double precision pi(ldpi,*),vx(ldpi,*),ecgm,ecrpa
        integer info,iri,jri
        double precision evals(nri),temp
        double precision work(nri,nri),work2(nri,nri)

        ecgm = 0d0
        ecrpa = 0d0
        iri = (nri*(nri+1))/2/ldpi

        call dpftri('n','l',nri,vx,info)

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

        ecrpa = ecrpa + dble(nri)

        end subroutine
