C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Initialize the occupation functions
C>
C> The occupation functions are initialized so that the electrons
C> are distributed across the natural orbitals according to Fermi.
C> The temperature used is sufficiently high to avoid zero valued
C> occupation numbers.
C>
C> 1. Work the occupation numbers out according to Fermi
C>    a. Keep increasing the temperature until the highest orbital
C>       has at least 1.0d-7 electrons and the lowest energy orbital
C>       has at most 1-1.0d-7 electron.
C>
C> 2. Compute the occupation numbers of the virtual orbitals as
C>    1-occ(ii) of the occupatied orbitals and copy the occupation
C>    numbers of the occupied orbitals
C>
C> 3. Scale the occupation numbers down by the number of electrons
C>
C> 4. Take the square root of the occupation numbers as a guess for 
C>    the coefficients.
C>
C> 5. Use sine functions to choose the signs of the coefficients
C>
C> 6. Modified Gramm-Schmidt to orthonormalize the set
C> 
      subroutine wfn1_nxt2_init_occ_fermi(nbf,nea,neb,ncore,nvirt,
     &                                    ehfa,ehfb,oa,ob)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      integer ncore(2) !< [Input] The number of core occupation
                       !< functions
      integer nvirt(2) !< [Input] The number of virtual occupation
                       !< functions
c
      double precision ehfa(nbf)   !< [Input] The alpha orbital 
                                   !< energies
      double precision ehfb(nbf)   !< [Input] The beta orbital 
                                   !< energies
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
      integer it,itt         !< Scratch variables
      double precision delta !< The width of an interval
      double precision x     !< The position associated with an
                             !< occupation function coefficient
      double precision dd    !< The overlap
      double precision pi    !< Pi
      double precision t     !< The temperature
      double precision Em    !< The Mermin energy
c
      double precision ddot
      external         ddot
c
      pi = acos(-1.0d0)
c
c     1. Work out the Fermi occupation numbers
c
      oa(1,1)   = 1.0d0
      oa(nbf,1) = 0.0d0
      t         = 0.02d0
      do while(oa(1,1).gt.(1.0d0-1.0d-5).or.oa(nbf,1).lt.1.0d-5)
        call wfn1_nxt_fermi(ehfa,oa(1,1),Em,it,itt,t,nea,nbf,6,.false.)
        t=t*1.05d0
      enddo
c
      ob(1,1)   = 1.0d0
      ob(nbf,1) = 0.0d0
      t         = 0.02d0
      do while(ob(1,1).gt.(1.0d0-1.0d-5).or.ob(nbf,1).lt.1.0d-5)
        call wfn1_nxt_fermi(ehfb,ob(1,1),Em,it,itt,t,neb,nbf,6,.false.)
        t=t*1.05d0
      enddo
c
c     2. Work out the virtual occupation numbers and copy the others
c
      do k = nea+1,nbf
        do i = 1, nbf
          oa(i,k) = 1.0d0-oa(i,1)
        enddo
      enddo
      do k = 2,nea
        do i = 1, nbf
          oa(i,k) = oa(i,1)
        enddo
      enddo
c
      do k = neb+1,nbf
        do i = 1, nbf
          ob(i,k) = 1.0d0-ob(i,1)
        enddo
      enddo
      do k = 2,neb
        do i = 1, nbf
          ob(i,k) = ob(i,1)
        enddo
      enddo
c
c     3.-4. Scale and take the square root
c
      do k = 1, nbf
        do i = 1, nbf
          oa(i,k) = sqrt(oa(i,k)/dble(nea))
          ob(i,k) = sqrt(ob(i,k)/dble(neb))
        enddo
      enddo
c
c     Setup the core and virtual occupation functions
c
      call dfill(ncore(1)*nbf,0.0d0,oa,1)
      call dfill(nvirt(1)*nbf,0.0d0,oa(1,nbf+1-nvirt(1)),1)
      do i = 1, ncore(1)
        call dfill(nbf,0.0d0,oa(i,1),nbf)
      enddo
      do i = nbf+1-nvirt(1), nbf
        call dfill(nbf,0.0d0,oa(i,1),nbf)
      enddo
      do i = 1, ncore(1)
        oa(i,i) = 1.0d0
      enddo
      do i = nbf+1-nvirt(1), nbf
        oa(i,i) = 1.0d0
      enddo
c
      call dfill(ncore(2)*nbf,0.0d0,ob,1)
      call dfill(nvirt(2)*nbf,0.0d0,ob(1,nbf+1-nvirt(2)),1)
      do i = 1, ncore(2)
        call dfill(nbf,0.0d0,ob(i,1),nbf)
      enddo
      do i = nbf+1-nvirt(2), nbf
        call dfill(nbf,0.0d0,ob(i,1),nbf)
      enddo
      do i = 1, ncore(2)
        ob(i,i) = 1.0d0
      enddo
      do i = nbf+1-nvirt(2), nbf
        ob(i,i) = 1.0d0
      enddo
cDEBUG
      write(*,*)'*** fermi: ',ncore(1),nvirt(1)
      call hess_hssout(oa,nbf,nbf,nbf)
      write(*,*)'*** fermi: ',ncore(2),nvirt(2)
      call hess_hssout(ob,nbf,nbf,nbf)
cDEBUG
c
c     Divide the domain [0,Pi] into nea pieces. Each occupation 
c     coefficient corresponds to the mid-point of one interval. We
c     will use the sign of sin(x_i) to decide the sign the coefficient.
c
      delta = pi/max(1,nea-ncore(1))
      do k = 1, nea-ncore(1)
        do i = 1, nea-ncore(1)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            oa(ncore(1)+i,ncore(1)+k) = -oa(ncore(1)+i,ncore(1)+k)
          endif
        enddo
      enddo
c
      do k = 1, nbf-nea-nvirt(1)
        do i = 1, nea-ncore(1)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            oa(ncore(1)+i,nea+k) = -oa(ncore(1)+i,nea+k)
          endif
        enddo
      enddo
c
      delta = pi/max(1,nbf-nea-nvirt(1))
      do k = 1, nea-ncore(1)
        do i = 1, nbf-nea-nvirt(1)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            oa(nea+i,ncore(1)+k) = -oa(nea+i,ncore(1)+k)
          endif
        enddo
      enddo
      do k = 1, nbf-nea-nvirt(1)
        do i = 1, nbf-nea-nvirt(1)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            oa(nea+i,nea+k) = -oa(nea+i,nea+k)
          endif
        enddo
      enddo
c
      delta = pi/max(1,neb-ncore(2))
      do k = 1, neb-ncore(2)
        do i = 1, neb-ncore(2)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            ob(ncore(2)+i,ncore(2)+k) = -ob(ncore(2)+i,ncore(2)+k)
          endif
        enddo
      enddo
      do k = 1, nbf-neb-nvirt(2)
        do i = 1, neb-ncore(2)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            ob(ncore(2)+i,neb+k) = -ob(ncore(2)+i,neb+k)
          endif
        enddo
      enddo
c
      delta = pi/max(1,nbf-neb-nvirt(2))
      do k = 1, neb-ncore(2)
        do i = 1, nbf-neb-nvirt(2)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            ob(neb+i,ncore(2)+k) = -ob(neb+i,ncore(2)+k)
          endif
        enddo
      enddo
      do k = 1, nbf-neb-nvirt(2)
        do i = 1, nbf-neb-nvirt(2)
          x = (i-0.5d0)*k*delta
          if (sin(x).lt.0.0d0) then
            ob(neb+i,neb+k) = -ob(neb+i,neb+k)
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
