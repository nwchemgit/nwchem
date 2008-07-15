      function leps_gradient(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
      logical leps_gradient
      integer geom,natom, nelem 
      integer h_c,i_c,h_g,i_g 
      integer h_q,i_q
      integer h_t,i_t

c      double precision gradient(3,3) ! can do this because LEPS is only for triatomic collinier sys
c
c
      write(*,*) "doing leps_gradient"
      call leps_energy(rtdb)
c     create geometry object
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('leps_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
C     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('leps_getgeom: no geometry',
     & 0,RTDB_ERR)
C   get the number of centers (atoms)
      if (.not. geom_ncent(geom,natom))
     & call errquit('leps_getgeom: ncent error',
     & 0,GEOM_ERR)
c     print tgeometry
      if (.not. geom_print(geom))
     & call errquit('leps_getgeom: print error',
     & 0,GEOM_ERR)
C ---- allocate local copies for coordinates and gradients
      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
     & h_c,i_c))
     & call errquit('leps_getgeom: coordinate not allocated',
     & 3*natom,MA_ERR)

      nelem=3*natom

      if (.not. ma_alloc_get(mt_dbl, nelem, 'gradient' ,
     & h_g,i_g))
     & call errquit('leps_getgradient: gradiant not allocated',
     & 3*natom,MA_ERR)

C  ---- load coordinates from the database
      if (.not. ma_alloc_get(mt_dbl, natom, 'charge', h_q,i_q))
     & call errquit('leps_getgeom: charge not allocated',
     & natom,MA_ERR)
c      if (.not.ma_alloc_get(mt_dbl, natom, 'mass', h_m,i_m))
c     & call errquit('leps_getgeom: mass not allocated',
c     & natom,MA_ERR)
      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
     & call errquit('leps_getgeom: tags not allocated',
     & natom,MA_ERR)

C-----Finished allocating memory for info we get from database

C ----- Fill up the local copies with data from geometry object
C------- accessed by the geom handle
      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_c), dbl_mb(i_q)))
     & call errquit('leps_getgeom: coord error',0,GEOM_ERR)

C ------ Calculate the gradient by giving the coordinates as input and 
      call leps_gradient1(natom,nelem,dbl_mb(i_c),dbl_mb(i_g))
      call leps_gradient0(natom,nelem,dbl_mb(i_g))  

C-------- Store the gradient into the database
      if (.not. rtdb_put(rtdb,'leps:gradient',mt_dbl,nelem,
     &     dbl_mb(i_g)))
     & call errquit('leps_gradient:  no grad ', 0,
     & RTDB_ERR)


      write(*,*) "looks like no error in doing leps_gradient"

      return
      end 

      function leps_gradient0(n,nel,grad)
C THis function just prints the coordinates, masses and tags(names for the atoms)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"

      logical leps_gradient0
      integer n,nel
      double precision grad(3,n)
c
      integer i
      write(*,*) 'Inside leps_gradien0:', nel

      do i=1,n
       write(*,*) i,grad(1,i),grad(2,i),grad(3,i)
      end do

      return
      END

      function leps_gradient1(n,nel,c0,grad)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
      logical leps_gradient1
      integer n,nel ! number of atoms
      integer i,j,jat
      double precision c0(3,n),grad(3,n)


c------from leps energy
      double precision c(3,n),dr(3,3,3),dvk(3,3,3),dvj(3,3,3)
      double precision r(3),vj(3),sqrj
      double precision de(3),alp(3),re(3),del(3)
      double precision tobohr, tohartree
      parameter (tobohr=0.529177249d+00)
      parameter (tohartree=627.509451d+00)
      data de     /109.4d+00,109.4d+00,109.4d+00/
      data alp    /1.943d+00,1.943d+00,1.943d+00/
      data re     /0.741d+00,0.741d+00,0.741d+00/
      data del    /0.135d+00,0.135d+00,0.135d+00/

C----- convert the coordinates to angstrom (because alp and re are in angstroms)
      do i=1,n
        do j=1,n
          c(j,i)=c0(j,i)*tobohr
        enddo
      end do
      r(1)=sqrt( (c(1,1)-c(1,2))**2+(c(2,1)-c(2,2))**2
     1                             +(c(3,1)-c(3,2))**2 )
      r(2)=sqrt( (c(1,2)-c(1,3))**2+(c(2,2)-c(2,3))**2
     1                             +(c(3,2)-c(3,3))**2 )
      r(3)=sqrt( (c(1,1)-c(1,3))**2+(c(2,1)-c(2,3))**2
     1                             +(c(3,1)-c(3,3))**2 )

c ------- Just cause I like all my arrays to be set to zero before use  
      do jat=1,n
         vj(jat)=0.0
         do j=1,3
            do i=1,3
               dr(i,j,jat)=0.0
               dvk(i,j,jat)=0.0
               dvj(i,j,jat)=0.0
            enddo
         enddo
      enddo
c      write(6,*) 'r 0 = ',r
c      write(6,*) 'dr 0 = ',dr
C derivatives of the coordinates dr(a,b,c) is der(r_a)/der(b_c)
C r1
      dr(1,1,1)=(c(1,1)-c(1,2))/r(1)
      dr(1,2,1)=(c(2,1)-c(2,2))/r(1)
      dr(1,3,1)=(c(3,1)-c(3,2))/r(1)
      dr(1,1,2)=-dr(1,1,1)
      dr(1,2,2)=-dr(1,2,1)
      dr(1,3,2)=-dr(1,3,1)
c r2
      dr(2,1,2)=(c(1,2)-c(1,3))/r(2)
      dr(2,2,2)=(c(2,2)-c(2,3))/r(2)
      dr(2,3,2)=(c(3,2)-c(3,3))/r(2)
      dr(2,1,3)=-dr(2,1,2)
      dr(2,2,3)=-dr(2,2,2)
      dr(2,3,3)=-dr(2,3,2)
c r3
      dr(3,1,1)=(c(1,1)-c(1,3))/r(3)
      dr(3,2,1)=(c(2,1)-c(2,3))/r(3)
      dr(3,3,1)=(c(3,1)-c(3,3))/r(3)
      dr(3,1,3)=-dr(3,1,1)
      dr(3,2,3)=-dr(3,2,1)
      dr(3,3,3)=-dr(3,3,1)

c      write(6,*) 'dr  = ',dr
      do i=1,3
            vj(i)=  (de(i)/(4.0*(1.0+del(i))))*
     1     (     (1.0  +3.0*del(i))*exp(-2.0*alp(i)*(r(i)-re(i)))
     2      -2.0*(3.0+      del(i))*exp(-    alp(i)*(r(i)-re(i))) )
      enddo


      do jat=1,n
        do j=1,3
          do i=1,3
            dvk(i,j,jat)=  (de(i)/(4.0*(1.0+del(i))))*
     1     (     (3.0+      del(i))*exp(-2.0*alp(i)*(r(i)-re(i)))
     2                               *   (-2.0*alp(i)* dr(i,j,jat))
     3      -2.0*(1.0  +3.0*del(i))*exp(-    alp(i)*(r(i)-re(i)))
     4                               *   (-    alp(i)* dr(i,j,jat)) )
            dvj(i,j,jat)=  (de(i)/(4.0*(1.0+del(i))))*
     1     (     (1.0  +3.0*del(i))*exp(-2.0*alp(i)*(r(i)-re(i)))
     2                               *   (-2.0*alp(i)* dr(i,j,jat))
     3      -2.0*(3.0+      del(i))*exp(-    alp(i)*(r(i)-re(i)))
     4                               *   (-    alp(i)* dr(i,j,jat)) )
          enddo
        enddo
      enddo
         sqrj=  sqrt(((vj(1)-vj(2))**2 + (vj(2)-vj(3))**2 +
     1                (vj(1)-vj(3))**2                      )/2.0 )
c         write(*,*) 'sqrj :', sqrj
c           do jat=1,n
c             do j=1,3
c                do i=1,3
c                write(*,*) 'jat,j,i,dvk',jat,j,i,dvk(i,j,jat)
c                write(*,*) 'jat,j,i,dvJ',jat,j,i,dvj(i,j,jat)  
c                enddo
c             enddo
c           enddo
         do jat=1,n
            do j=1,3
               grad(j,jat)=(   dvk(1,j,jat)+dvk(2,j,jat)+dvk(3,j,jat)
     1                   -( dvj(1,j,jat)*(vj(1)-vj(2))
     2                     -dvj(2,j,jat)*(vj(1)-vj(2))
     3                     +dvj(2,j,jat)*(vj(2)-vj(3))
     4                     -dvj(3,j,jat)*(vj(2)-vj(3))
     5                     +dvj(1,j,jat)*(vj(1)-vj(3))
     6                     -dvj(3,j,jat)*(vj(1)-vj(3)) )/sqrj/2.0)
     7                   /tohartree*tobohr

            enddo
         enddo
              
         write (*,*) 'grad 0: ',grad
         do jat=1,n
            do j=1,3
               grad(j,jat)=(   dvk(1,j,jat)+dvk(2,j,jat)+dvk(3,j,jat)
     1                   -( (dvj(1,j,jat)
     2                     -dvj(2,j,jat))*(vj(1)-vj(2))
     3                     +(dvj(2,j,jat)
     4                     -dvj(3,j,jat))*(vj(2)-vj(3))
     5                     +(dvj(1,j,jat)
     6                     -dvj(3,j,jat))*(vj(1)-vj(3)) )/sqrj/2.0)

c     1                   -( dvj(1,j,jat)*(vj(1)-vj(2))
c     2                     -dvj(2,j,jat)*(vj(1)-vj(2))
c     3                     +dvj(2,j,jat)*(vj(2)-vj(3))
c     4                     -dvj(3,j,jat)*(vj(2)-vj(3))
c     5                     +dvj(1,j,jat)*(vj(1)-vj(3))
c     6                     -dvj(3,j,jat)*(vj(1)-vj(3)) )/sqrj/2.0)
     7                   /tohartree*tobohr
            enddo
         enddo
      write (*,*) 'grad 1: ',grad
c      write (*,*) 'test - grad', 0.0+0.0+0.0 
      return
      end

      subroutine leps_gradient_test(n,nel,c0,gradient)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c      logical leps_gradient1
      integer n,nel ! number of atoms
      integer i,j
      double precision c0(3,n),gradient(3,n),cnt
    

c------from leps energy

      write(*,*) 'Inside leps_gradien1:'
C----- Testing with rangom values
      cnt=0.0
      do i=1,n
        do j=1,n
             cnt=cnt+1.0
             gradient(i,j)=cnt
        enddo
      enddo
C-end of test
      return 
      end 
