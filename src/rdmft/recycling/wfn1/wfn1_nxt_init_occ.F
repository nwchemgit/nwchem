C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Initialize the occupation functions
C>
C> The occupation functions are initialized so that the electrons
C> are evenly distributed across the natural orbitals. No occupation
C> function coefficients should be hard zeroes. Of course the 
C> occupation functions should also be orthogonal. The approach to 
C> get to this is to
C>
C> 1. Set all coefficients to 1
C>
C> 2. Use sine functions to choose the signs of the coefficients
C>
C> 3. If needed use Modified Gramm-Schmidt to orthonormalize the set
C> 
      subroutine wfn1_nxt_init_occ(nbf,ncore,nvirt,oa,ob)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      integer ncore(2) !< [Input] The number of core occupation 
                       !< functions
      integer nvirt(2) !< [Input] The number of virtual occupation 
                       !< functions
c
      double precision oa(nbf,nbf) !< [Output] The alpha occupation
                                   !< functions
      double precision ob(nbf,nbf) !< [Output] The beta occupation
                                   !< functions
c
      integer k              !< Counter over frequencies
      integer i              !< Counter over coefficients
      integer p              !< Counter over occupation functions
      integer q              !< Counter over occupation functions
      double precision delta !< The width of an interval
      double precision x     !< The position associated with an
                             !< occupation function coefficient
      double precision dd    !< The overlap
      double precision pi    !< Pi
c
      double precision ddot
      external         ddot
c
      pi = acos(-1.0d0)
c
c     Alpha: Fill the occupation functions with unit coefficients
c
      call dfill(nbf*nbf,     1.0d0,oa,1)
c
c     Alpha: Zero the columns for the core and virtual occupation
c            functions
c
      call dfill(ncore(1)*nbf,0.0d0,oa,1)
      call dfill(nvirt(1)*nbf,0.0d0,oa(1,nbf+1-nvirt(1)),1)
c
c     Alpha: Zero the rows for the core and virtual occupation
c            functions
c
      do i = 1, ncore(1)
        call dfill(nbf,0.0d0,oa(i,1),nbf)
      enddo
      do i = nbf+1-nvirt(1), nbf
        call dfill(nbf,0.0d0,oa(i,1),nbf)
      enddo
c
c     Alpha: Set the diagonal elements
c
      do i = 1, ncore(1)
        oa(i,i) = 1.0d0
      enddo
      do i = nbf+1-nvirt(1), nbf
        oa(i,i) = 1.0d0
      enddo
c
c     Beta: Fill the occupation functions with unit coefficients
c
      call dfill(nbf*nbf,     1.0d0,ob,1)
c
c     Beta: Zero the columns for the core and virtual occupation
c           functions
c
      call dfill(ncore(2)*nbf,0.0d0,ob,1)
      call dfill(nvirt(2)*nbf,0.0d0,ob(1,nbf+1-nvirt(2)),1)
c
c     Beta: Zero the rows for the core and virtual occupation
c           functions
c
      do i = 1, ncore(2)
        call dfill(nbf,0.0d0,ob(i,1),nbf)
      enddo
      do i = nbf+1-nvirt(2), nbf
        call dfill(nbf,0.0d0,ob(i,1),nbf)
      enddo
c
c     Beta: Set the diagonal elements
c
      do i = 1, ncore(2)
        ob(i,i) = 1.0d0
      enddo
      do i = nbf+1-nvirt(2), nbf
        ob(i,i) = 1.0d0
      enddo
c
c     Divide the domain [0,Pi] into nbf pieces. Each occupation 
c     coefficient corresponds to the mid-point of one interval. We
c     will use the sign of sin(x_i) to decide the sign the coefficient.
c
      delta = pi/max(1,nbf-ncore(1)-nvirt(1))
      do k = 1, nbf-ncore(1)-nvirt(1)
        do i = 1, nbf-ncore(1)-nvirt(1)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            oa(ncore(1)+i,ncore(1)+k) = -oa(ncore(1)+i,ncore(1)+k)
          endif
        enddo
      enddo
c
      delta = pi/max(1,nbf-ncore(2)-nvirt(2))
      do k = 1, nbf-ncore(2)-nvirt(2)
        do i = 1, nbf-ncore(2)-nvirt(2)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            ob(ncore(2)+i,ncore(2)+k) = -ob(ncore(2)+i,ncore(2)+k)
          endif
        enddo
      enddo
c
c     Now do Modified Gramm-Schmidt
c
c     Alpha occupation functions
c
      do p = 1, nbf
        i = 0
 10     i = i + 1
        dd = ddot(nbf,oa(1,p),1,oa(1,p),1)
        dd = 1.0d0/sqrt(dd)
        call dscal(nbf,dd,oa(1,p),1)
        do q = 1, p-1
          dd = ddot(nbf,oa(1,q),1,oa(1,p),1)
          call daxpy(nbf,-dd,oa(1,q),1,oa(1,p),1)
        enddo
        dd = ddot(nbf,oa(1,p),1,oa(1,p),1)
        if (dd.lt.0.5d0) then
          if (i.lt.2) then
            go to 10
          else
            write(*,*)'This is not good: p iter dda',p,i,dd
          endif
        endif
        dd = 1.0d0/sqrt(dd)
        call dscal(nbf,dd,oa(1,p),1)
      enddo
c
c     Beta occupation functions
c
      do p = 1, nbf
        i = 0
 20     i = i + 1
        dd = ddot(nbf,ob(1,p),1,ob(1,p),1)
        dd = 1.0d0/sqrt(dd)
        call dscal(nbf,dd,ob(1,p),1)
        do q = 1, p-1
          dd = ddot(nbf,ob(1,q),1,ob(1,p),1)
          call daxpy(nbf,-dd,ob(1,q),1,ob(1,p),1)
        enddo
        dd = ddot(nbf,ob(1,p),1,ob(1,p),1)
        if (dd.lt.0.5d0) then
          if (i.lt.2) then
            go to 20
          else
            write(*,*)'This is not good: p iter ddb',p,i,dd
          endif
        endif
        dd = 1.0d0/sqrt(dd)
        call dscal(nbf,dd,ob(1,p),1)
      enddo
cDEBUG
c     call dfill(nbf*nbf,0.0d0,oa,1)
c     call dfill(nbf*nbf,0.0d0,ob,1)
c     do i = 1, nbf
c       oa(i,i) = 1.0d0
c       ob(i,i) = 1.0d0
c     enddo
cDEBUG
c
      end
C>
C> @}
