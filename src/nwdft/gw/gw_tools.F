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
