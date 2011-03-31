*
*     $Id$
*

*     ***************************
*     *                         *
*     *     zero_cell_contraint *
*     *                         *
*     ***************************
*
*   This routine adjusts the cell gradients so
* that specified constraints are satisfied.
*
*   strategy==1 - a,b,c,alpha,beta,gamma constraints
*   strategy==2 - lattice vector constraints
*
      subroutine zero_cell_constraint(stress,lstress)
      implicit none
      real*8 stress(3,3)
      real*8 lstress(6)

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      integer i,u,v

*     **** external functions ****
      integer  control_optimize_cell_strategy
      integer  control_optimize_lattice
      integer  control_optimize_lattice_vectors
      external control_optimize_cell_strategy
      external control_optimize_lattice
      external control_optimize_lattice_vectors

      call Parallel_taskid(taskid)

*      **** apply lattice constraints ****
       if (control_optimize_cell_strategy().eq.1) then

*         *** dE/da ==0,dE/db==0 .... dE/dgamma==0 ***
          do i=1,6
            if (control_optimize_lattice(i).eq.0) then
              lstress(i) = 0.0d0
              if (taskid.eq.MASTER) write(*,1800) i
            end if
          end do


*      **** apply lattice_vector constraints ****
       else if (control_optimize_cell_strategy().eq.2) then
          do v=1,3
          do u=1,3
            if (control_optimize_lattice_vectors(u,v).eq.0)  then
              stress(u,v) = 0.0d0
              if (taskid.eq.MASTER) write(*,1801) u,v
            end if
          end do
          end do
       end if

      return
 1800 FORMAT('Lattice Constraint: gradient(',i1,')-->0')
 1801 FORMAT('Lattice Constraint: gradient(',i1,',',i1,')-->0')
      end

*     *******************************
*     *                             *
*     *     cell_lattice_gradient   *
*     *                             *
*     *******************************
*
*  Given the gradient, stress, wrt to lattice vectors
*  this routine calculates the gradient,lstress, wrt a,b,c,alpha,beta,gamma.
*  This is a vector projection
*
      subroutine cell_lattice_gradient(stressin,lstress)
      implicit none
! Calling args
      real*8 stressin(3,3) ! Input
      real*8 lstress(6)    ! Output (a,b,c,alpha,beta,gamma)
! Variables
      real*8 amat(3,3)     ! amat({x,y,z},{a,b,c})
      real*8 stress(3,3)   ! stress({x,y,z},{a,b,c})
      real*8 vec(3,3),vecinv(3,3)
      real*8 length(3)
      integer ii,jj,kk,ll,i
      integer iangle
      real*8 dot1, dot2
      real*8 pi, alpha, beta, gamma, a, b, c
      real*8 s1,s2,s3,c1,c2,c3,vol,dum,siner(4:6),coser(4:6)
      real*8 dum1,dum2,dum3
! External functions
      real*8 lattice_unita
      external lattice_unita

!     Initialize

      do ii = 1, 6
        lstress(ii) = 0.0d0
      enddo

      do ii = 1, 3
        do jj = 1, 3
          stress(jj,ii) = stressin(jj,ii)
          if(abs(stress(jj,ii)) .lt. 1.0e-10) stress(jj,ii)=0
        enddo
      enddo

      do jj = 1, 3
       do ii = 1, 3
         amat(ii,jj) = lattice_unita(ii,jj)
         if(abs(amat(ii,jj)) .lt. 1.0e-10) amat(ii,jj)=0
       end do
      enddo

      call lattice_abc_abg(a,b,c,alpha,beta,gamma)
      pi = 4.0d0*datan(1.0d0)
      alpha = alpha*pi/180.0d0
      beta  = beta *pi/180.0d0
      gamma = gamma*pi/180.0d0
      length(1)=a
      length(2)=b
      length(3)=c

      c1=dcos(alpha)
      c2=dcos(beta)
      c3=dcos(gamma)
      s1=dsin(alpha)
      s2=dsin(beta)
      s3=dsin(gamma)
      coser(4) = c1
      coser(5) = c2
      coser(6) = c3
      siner(4) = s1
      siner(5) = s2
      siner(6) = s3

      vol = (a*b*c)*sqrt(1-(c1**2+c2**2+c3**2)+(2*c1*c2*c3))

      do ii = 1, 3 ! a,b,c
        do jj = 1, 3 ! X, Y, Z
          dum = amat(jj,ii)*(stress(jj,ii)/length(ii)) ! Dot product
          lstress(ii) = lstress(ii) + dum
        end do
      end do

      do ii = 1, 3 ! a,b,c
        do iangle = 4, 6
          if (ii .ne. iangle-3) then ! d(a)/d(alpha) = 0
            ! For finding db(x,y,z)/dalpha, we have 3 equations and 3 unknowns
            ! -bc*sin(alpha)=cx*dbx + cy*dby + cz*dbz ! From derivative of alpha
            !              0=bx*dbx + by*dby + bz*dbz ! From derivative of b
            !              0=ax*dbx + ay*dby + az*dbz ! From derivative of gamma
            kk = 9 - (iangle+ii)  ! other atom in this angle
            jj = 6 - (kk+ii)      ! other atom not in this angle
            do ll = 1, 3
              vec(1,ll) = amat(ll,kk)
              vec(2,ll) = amat(ll,ii)
              vec(3,ll) = amat(ll,jj)
            enddo
            ! This is Ax = b problem
            ! A = "vec"
            ! b = #, 0, 0
            ! x = the solution we want for dbx,dby,dbz
            call invert3(vec,vecinv)
            ! A^-1 * Ax = A^-1*b
            !         x = A^-1*b
            do i=1, 3
             vec(i,1)=vecinv(i,1)*(-length(ii)*length(kk)*siner(iangle))
            enddo
            do i = 1, 3
             lstress(iangle) = lstress(iangle) + vec(i,1)*stress(i,ii)/2
            enddo
          endif
        enddo
      enddo
      return
      end


      subroutine invert3(vector,vecinv)
      implicit none
      real*8 vector(9)
      real*8 vecinv(9)
      real*8 deter

      deter=(vector(5)*vector(9)-vector(6)*vector(8))*vector(1)+
     >      (vector(6)*vector(7)-vector(4)*vector(9))*vector(2)+
     >      (vector(4)*vector(8)-vector(5)*vector(7))*vector(3)
      vecinv(1)=(vector(5)*vector(9)-vector(6)*vector(8))/deter
      vecinv(2)=(vector(8)*vector(3)-vector(9)*vector(2))/deter
      vecinv(3)=(vector(2)*vector(6)-vector(3)*vector(5))/deter
      vecinv(4)=(vector(6)*vector(7)-vector(4)*vector(9))/deter
      vecinv(5)=(vector(9)*vector(1)-vector(7)*vector(3))/deter
      vecinv(6)=(vector(3)*vector(4)-vector(1)*vector(6))/deter
      vecinv(7)=(vector(4)*vector(8)-vector(5)*vector(7))/deter
      vecinv(8)=(vector(7)*vector(2)-vector(8)*vector(1))/deter
      vecinv(9)=(vector(1)*vector(5)-vector(2)*vector(4))/deter
      return
      end

*     *******************************
*     *                             *
*     *     cell_unita_abc_abg      *
*     *                             *
*     *******************************
*
*  This routine calculates a,b,c,alpha,beta,gamma from unita
*
      subroutine cell_unita_abc_abg(unita,lattice)
      implicit none
      real*8 unita(3,3),lattice(6)

      real*8 a,b,c
      real*8 alpha,beta,gamma

*     *** local variables ****
      real*8 d2

*     **** determine a,b,c,alpha,beta,gmma ***
      !pi = 4.0d0*datan(1.0d0)
      a = dsqrt(unita(1,1)**2
     >        + unita(2,1)**2
     >        + unita(3,1)**2)
      b = dsqrt(unita(1,2)**2
     >        + unita(2,2)**2
     >        + unita(3,2)**2)
      c = dsqrt(unita(1,3)**2
     >        + unita(2,3)**2
     >        + unita(3,3)**2)

      d2 = (unita(1,2)-unita(1,3))**2
     >   + (unita(2,2)-unita(2,3))**2
     >   + (unita(3,2)-unita(3,3))**2
      alpha = (b*b + c*c - d2)/(2.0d0*b*c)
      alpha = dacos(alpha)

      d2 = (unita(1,3)-unita(1,1))**2
     >   + (unita(2,3)-unita(2,1))**2
     >   + (unita(3,3)-unita(3,1))**2
      beta = (c*c + a*a - d2)/(2.0d0*c*a)
      beta = dacos(beta)

      d2 = (unita(1,1)-unita(1,2))**2
     >   + (unita(2,1)-unita(2,2))**2
     >   + (unita(3,1)-unita(3,2))**2
      gamma = (a*a + b*b - d2)/(2.0d0*a*b)
      gamma = dacos(gamma)

      lattice(1) = a
      lattice(2) = b
      lattice(3) = c
      lattice(4) = alpha
      lattice(5) = beta
      lattice(6) = gamma
      return
      end


*     *******************************
*     *                             *
*     *     cell_abc_abg_unita      *
*     *                             *
*     *******************************
*
*  This routine calculates unita from a,b,c,alpha,beta,gamma 
*
*  Adapted from geom_3d code.
c
c-----------> build a-matrix
c
c Chosen convention of two axis systems.
c    Orient oblique crystal system to external rectilinear laboratory
c    coordinate system such that
c    b-axis is along j (unit rectilinear vector), ab-plane contains i (rect),
c    i.e., j along b, i in  ab-plane
c
c  Note this matrix is implicitly built to convert the input coordinates
c  from angstroms to bohr; i.e., amat is in atomic units.
c
c  the Jth column of this matrix contains the jth lattice vector
c aij is the ith Cartesian component of the jth lattice vector.
c
c     Entry - lattice, lattice constants in au and lattice angles in radians
c
c     Exit  - lattice_a, lattice vectors
c
      subroutine cell_abc_abg_unita(lattice,lattice_a)
      double precision lattice(6),lattice_a(3,3)

*     **** local variables ****
      integer i,j,iang
      double precision c1,c2,c3,s3
      double precision cdist(3),cang(3)
      double precision gmat(3,3)

*     **** external functions ****
      double precision deter3
      external deter3

      cdist(1) = lattice(1)
      cdist(2) = lattice(2)
      cdist(3) = lattice(3)
      cang(1)  = lattice(4)
      cang(2)  = lattice(5)
      cang(3)  = lattice(6)

c    **** build the metrical matrix (atomic units) ****
      do i=1,3
        gmat(i,i)=cdist(i)**2
      end do
      iang=3
      do i=1,3
        do j=i+1,3
          gmat(i,j)=cdist(i)*cdist(j)*dcos(cang(iang))
          gmat(j,i)=gmat(i,j)
          iang=iang-1
        end do  
      end do

*     **** get volume ***
      vol=dsqrt(deter3(gmat))

*     **** generate lattice_a ***
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
      lattice_a(1,1) = cdist(1)*s3
      lattice_a(1,2) = 0.0d0
      lattice_a(1,3) = (cdist(3)*(c2-c1*c3)/s3)
      lattice_a(2,1) = cdist(1)*c3
      lattice_a(2,2) = cdist(2)
      lattice_a(2,3) = cdist(3)*c1
      lattice_a(3,1) = 0.0d0
      lattice_a(3,2) = 0.0d0
      lattice_a(3,3) = (vol/(cdist(1)*cdist(2)*s3))
      return
      end

