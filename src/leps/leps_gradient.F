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
      logical leps_energy
      external leps_energy
      integer geom,natom, nelem 
      integer h_c,i_c,h_g,i_g 
      integer h_q,i_q
      integer h_t,i_t
      double precision eleps
c      integer h_g1,i_g1
c      double precision gg(3,3),c(3,3),eleps
c      double precision gradient(3,3) ! can do this because LEPS is only for triatomic collinier sys
c
c
c      write(*,*) "doing leps_gradient"
c      natom=3
c      eleps=0.0
      if (.not. leps_energy(rtdb))
     & call errquit('leps_gradien: error',
     & 0,GEOM_ERR)



c       call leps_energy(rtdb)
c      call leps_energy1(natom,c,eleps)
c      if (.not. rtdb_put(rtdb,'leps:energy',mt_dbl,1,eleps))
c     & call errquit('leps_energy:  no energy ', 0,
c     & RTDB_ERR)
c      write(*,*) "LEPS = ",eleps 
c----   create geometry object
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('leps_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
c---    load geometry into the object
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
c      call leps_gradient0(natom,nelem,dbl_mb(i_g))
c      call leps_gradient1(natom,nelem,c,gg)  
c      call leps_gradient1(natom,nelem,dbl_mb(i_c),gg)
c      call leps_gradient0(natom,nelem,gg)
C-------- Store the gradient into the database
      if (.not. rtdb_put(rtdb,'leps:gradient',mt_dbl,nelem,
     &     dbl_mb(i_g)))
     & call errquit('leps_gradient:  no grad ', 0,
     & RTDB_ERR)

      if(ga_nodeid().eq.0) then
         call leps_gradient0(natom,nelem,dbl_mb(i_g))
      endif

      if(.not.ma_free_heap(h_c))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_q))
     &     call errquit('memory deallocation:charge',0,RTDB_ERR)
      if(.not.ma_free_heap(h_g))
     &     call errquit('memory deallocation:mass',0,RTDB_ERR)
      if(.not.ma_free_heap(h_t))
     &     call errquit('memory deallocation:tags',0,RTDB_ERR)

      if (.not. geom_destroy(geom))
     & call errquit('leps_getgeom:  geom_destroy failed',
     & 911,GEOM_ERR)

c      if (.not. rtdb_put(rtdb,'leps:gradient',mt_dbl,nelem,
c     &     gg))
c     & call errquit('leps_gradient:  no grad ', 0,
c     & RTDB_ERR)

    
c      write(*,*) "looks like no error in doing leps_gradient"
C ---- Let's see if I can get the gradient from the database
c      if (.not. rtdb_get(rtdb,'leps:gradient',mt_dbl,nelem,
c     &     dbl_mb(i_g1)))
c     & call errquit('leps_gradient1:  no grad ', 0,
c     & RTDB_ERR)
c       write(*,*) "looks like properly stored in database"  
c       call leps_gradient0(natom,nelem,dbl_mb(i_g1))

C ---
      leps_gradient = .true.
      return
      end 

      subroutine leps_gradient0(n,nel,grad)
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

c      logical leps_gradient0
      integer n,nel
      double precision grad(3,n)
c
      integer i
      write(*,*) 'Leps_gradient:', nel

      do i=1,n
       write(6,9000) i,grad(1,i),grad(2,i),grad(3,i)
      end do
 9000 format (i4, F18.14,F18.14,F18.14)
      return
      END

      subroutine leps_gradient1(n,nel,c0,grad)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
      logical out1 
      integer n,nel ! number of atoms
      integer i,j,jat
      double precision small,dot,cosa,c0(3,3),grad(3,3)
      logical coll

c------from leps energy
      double precision c(3,3),dr(3,3,3),dvk(3,3,3),dvj(3,3,3)
      double precision r(3),e1(3),e3(3),vj(3),sqrj
      double precision d1(3),be1(3),re1(3),d3(3),be3(3),re3(3)
      double precision tobohr, tohartree
      parameter (tobohr=0.529177249d+00)
      parameter (tohartree=627.509451d+00)
      data d1     /234.524674,234.524674,64.925971/
      data be1    /0.929968,0.929968,0.432955/
      data re1    /1.776382,1.776382,2.094857/
      data d3     /220.244820,220.244820,284.999867/
      data be3    /4.822681,4.822681,1.016811/
      data re3    /1.785014,1.785014,2.186060/
      out1=.false.
      dot=0.0
      cosa=0.0 
      small = 1d-10
      coll=.false.
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

C check collinearity
      dot =      (c(1,1)-c(1,2))*(c(1,1)-c(1,3))+
     1           (c(2,1)-c(2,2))*(c(2,1)-c(2,3))+
     2           (c(3,1)-c(3,2))*(c(3,1)-c(3,3))
      cosa=dot/(r(1)*r(3))

      cosa=abs(cosa)
      coll=abs(1.0-cosa).lt.small


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
c      write(6,*) 'C =',c
c      write(6,*)'c0=',c0
C derivatives of the coordinates dr(a,b,c) is der(r_a)/der(b_c)
C r1
      dr(1,1,1)=(c(1,1)-c(1,2))/r(1)
      dr(1,2,1)=(c(2,1)-c(2,2))/r(1)
      dr(1,3,1)=(c(3,1)-c(3,2))/r(1)
      dr(1,1,2)=-dr(1,1,1)
      dr(1,2,2)=-dr(1,2,1)
      dr(1,3,2)=-dr(1,3,1)
c      write (*,*) 'dr(1,3,2):',dr(1,3,2)
c r2
      dr(2,1,2)=(c(1,2)-c(1,3))/r(2)
      dr(2,2,2)=(c(2,2)-c(2,3))/r(2)
      dr(2,3,2)=(c(3,2)-c(3,3))/r(2)
      dr(2,1,3)=-dr(2,1,2)
      dr(2,2,3)=-dr(2,2,2)
      dr(2,3,3)=-dr(2,3,2)
c      write (*,*) 'dr(2,3,3):',dr(2,3,3)
c r3
      dr(3,1,1)=(c(1,1)-c(1,3))/r(3)
      dr(3,2,1)=(c(2,1)-c(2,3))/r(3)
      dr(3,3,1)=(c(3,1)-c(3,3))/r(3)
      dr(3,1,3)=-dr(3,1,1)
      dr(3,2,3)=-dr(3,2,1)
      dr(3,3,3)=-dr(3,3,1)
c      write (*,*) 'dr(3,3,3):',dr(3,3,3)

c      write(6,*) 'dr  = ',dr
      if(coll) then
      do i=1,3
            e1(i) = d1(i)*(1.0 - exp(-be1(i)*(r(i)-re1(i))))**2-d1(i)
            e3(i) = d3(i)*(1.0 + exp(-be3(i)*(r(i)-re3(i))))**2-d3(i)
C            vk(i) = (e1(i) + e3(i))/2.0
            vj(i) = (e1(i) - e3(i))/2.0
      enddo


      do jat=1,n
        do j=1,3
          do i=1,3
            dvk(i,j,jat)= d1(i)*(1.0 - exp(-be1(i)*(r(i)-re1(i))))*
     1       (-1)*exp(-be1(i)*(r(i)-re1(i)))*(-be1(i))* dr(i,j,jat)
     2      +             (d3(i)*(1.0 - exp(-be3(i)*(r(i)-re3(i))))*
     3      (-exp(-be3(i)*(r(i)-re3(i))))*(-be3(i))* dr(i,j,jat))

            dvj(i,j,jat)= d1(i)*(1.0 - exp(-be1(i)*(r(i)-re1(i))))*
     1      (-exp(-be1(i)*(r(i)-re1(i))))*(-be1(i))* dr(i,j,jat)
     2      -             (d3(i)*(1.0 - exp(-be3(i)*(r(i)-re3(i))))*
     3      (-exp(-be3(i)*(r(i)-re3(i))))*(-be3(i))* dr(i,j,jat))
          enddo
        enddo
      enddo
         sqrj=  sqrt(((vj(1)-vj(2))**2 + (vj(2)-vj(3))**2 +
     1                (vj(1)-vj(3))**2                      )/2.0 )

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
              
c         write (*,*) 'grad 0: ',grad
c         do jat=1,n
c            do j=1,3
c               grad(j,jat)=(   dvk(1,j,jat)+dvk(2,j,jat)+dvk(3,j,jat)
c     1                   -( (dvj(1,j,jat)
c     2                     -dvj(2,j,jat))*(vj(1)-vj(2))
c     3                     +(dvj(2,j,jat)
c     4                     -dvj(3,j,jat))*(vj(2)-vj(3))
c     5                     +(dvj(1,j,jat)
c     6                     -dvj(3,j,jat))*(vj(1)-vj(3)) )/sqrj/2.0)

c     1                   -( dvj(1,j,jat)*(vj(1)-vj(2))
c     2                     -dvj(2,j,jat)*(vj(1)-vj(2))
c     3                     +dvj(2,j,jat)*(vj(2)-vj(3))
c     4                     -dvj(3,j,jat)*(vj(2)-vj(3))
c     5                     +dvj(1,j,jat)*(vj(1)-vj(3))
c     6                     -dvj(3,j,jat)*(vj(1)-vj(3)) )/sqrj/2.0)
c     7                   /tohartree*tobohr
c            enddo
c         enddo
        else
         do jat=1,n
            do j=1,3
               grad(j,jat)=0.0
            enddo
         enddo   
        endif
         if (out1) then
          write(*,*) 'sqrj :', sqrj
           write(*,*) 'vj :', vj
           write(*,*) 'r: ',r(1),r(2),r(3)
           write(*,*) 'dr :', dr
           write(*,*) 'dvk :',dvk
           write(*,*) 'dvj :', dvj
           write(*,*) 'c :',c
           write(*,*) 'c0 :',c0
           write (*,*) 'grad :',grad
         endif
c      write (*,*) 'grad 1: ',grad
c      write (*,*) 'test - grad', 0.0+0.0+0.0 
      return
      end

c $Id$
