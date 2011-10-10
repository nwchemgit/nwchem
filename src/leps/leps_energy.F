      function leps_energy(rtdb)
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
      logical leps_energy
c
      integer geom
      integer natom
      integer h_c,i_c
      integer h_m,i_m
      integer h_q,i_q
      integer h_t,i_t
      double precision eleps
c      write(*,*) "doing leps energy"
      eleps=0.0 
c     create geometry object 
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('leps_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
C      write(*,*) "doing leps energy 1"
C     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('leps_getgeom: no geometry',
     & 0,RTDB_ERR)
C      write(*,*) "doing leps energy 2"

C     get the number of centers 
      if (.not. geom_ncent(geom,natom))
     & call errquit('leps_getgeom: ncent error',
     & 0,GEOM_ERR)      
c      write(*,*) "doing leps energy 3"
c     print tgeometry
      if (.not. geom_print(geom))
     & call errquit('leps_getgeom: print error',
     & 0,GEOM_ERR)      
C      write(*,*) "doing leps energy 4"
c    
c     allocate memory for coordinates
      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
     & h_c,i_c))
     & call errquit('leps_getgeom: coordinate not allocated',
     & 3*natom,MA_ERR)
c      write(*,*) "doing leps energy 5"
      if (.not. ma_alloc_get(mt_dbl, natom, 'charge', h_q,i_q))
     & call errquit('leps_getgeom: charge not allocated',
     & natom,MA_ERR)
C      write(*,*) "doing leps energy 6"
      if (.not.ma_alloc_get(mt_dbl, natom, 'mass', h_m,i_m))
     & call errquit('leps_getgeom: mass not allocated',
     & natom,MA_ERR)
C      write(*,*) "doing leps energy 7"
      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
     & call errquit('leps_getgeom: tags not allocated',
     & natom,MA_ERR)
C      write(*,*) "doing leps energy 8"

C-----Finished allocating memory for info we get from database

C ----- Fill up the local copies with data from geometry object
C------- accessed by the geom handle
C      write(*,*) "doing leps energy 9"
      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_c), dbl_mb(i_q)))
     & call errquit('leps_getgeom: coord error',0,GEOM_ERR)
C      write(*,*) "doing leps energy 10"
      if (.not. geom_masses_get(geom, natom, dbl_mb(i_m)))
     & call errquit('leps_getgeom:  mass error',0,GEOM_ERR)
C--- now compute the leps energy
C      write(*,*) 'Before calling leps_energy 1:',eleps
      call leps_energy1(natom,dbl_mb(i_c),eleps)
C      write(*,*) 'After calling leps_energy 1:',eleps
C--- end of leps energy computation

C-----Print what you get out of the database
c      call leps_energy0(natom,dbl_mb(i_c),dbl_mb(i_m),byte_mb(i_t))
C-----Calculate the energy
      if (.not. rtdb_put(rtdb,'leps:energy',mt_dbl,1,eleps))
     & call errquit('leps_energy:  no energy ', 0,
     & RTDB_ERR)

       write(*,*) 'Leps energy: ',eleps
      if(.not.ma_free_heap(h_c))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_q))
     &     call errquit('memory deallocation:charge',0,RTDB_ERR)
      if(.not.ma_free_heap(h_m))
     &     call errquit('memory deallocation:mass',0,RTDB_ERR)
      if(.not.ma_free_heap(h_t))
     &     call errquit('memory deallocation:tags',0,RTDB_ERR)





      if (.not. geom_destroy(geom))
     & call errquit('leps_getgeom:  geom_destroy failed',
     & 911,GEOM_ERR)

      leps_energy = .true.

      return
      END

      subroutine leps_energy0(n,c,m,t)
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

c      logical leps_energy0
      integer n
      double precision c(3,3)
      double precision m(3)
      character*16 t(3)
c
      integer i

      do i=1,n
       write(*,*) i,c(1,i),c(2,i),c(3,i),m(i),t(i)
      end do    

      return
      END

      subroutine leps_energy_test(n,c0,eleps)
            implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c      double precision leps_energy_test
      integer n
      double precision c0(3,3),eleps
      eleps = 0.15 

      end

      subroutine leps_energy1(n,c0,eleps)
c      subroutine leps_energy1(n,c0,eleps)
C THis function calculates the leps energy and returns teh value
C n is the number of atoms, c is array holding the coordinates 
C coef for SN2 reaction Cl-Ch3-Cl-
C ref. J.Chem.Phys Vol86, No3, 1356. Bergsma et.al.
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"

c      double precision leps_energy1
      integer n
      double precision c(3,3),c0(3,3),eleps
      double precision r(3),e1(3),e3(3),vq(3),vj(3),vleps
      double precision dot,cosa,small
      double precision d1(3),be1(3),re1(3),d3(3),be3(3),re3(3)
      double precision tobohr, tohartree
      integer i,j 
      logical coll,out1
      parameter (tobohr=0.529177249d+00)
      parameter (tohartree=627.509451d+00)
      data d1     /234.524674,234.524674,64.925971/
      data be1    /0.929968,0.929968,0.432955/
      data re1    /1.776382,1.776382,2.094857/
      data d3     /220.244820,220.244820,284.999867/
      data be3    /4.822681,4.822681,1.016811/
      data re3    /1.785014,1.785014,2.186060/

      small = 1d-10
      out1=.false.
      coll=.false.
c
c      write (*,*) 'entering leps_energy1'
      eleps = 0.0
      if (n.ne.3) then
         write(*,*) "Number of atoms is not 3. Leps only works with 3"
      endif 
C----- convert the coordinates to angstrom (because alp and re are in angstroms)
      do i=1,3
        do j=1,3
          c(j,i)=c0(j,i)*tobohr
        enddo
      end do

C      do i=1,n
C       write(*,*) i,c(1,i),c(2,i),c(3,i)
C       write(*,*) i,c0(1,i),c0(2,i),c0(3,i)
C      end do



c   must make sure that r and c are in th right units
C distance should be in angstroms, de (kcal/mol), alp (1/angs), re(angs), del(no unit)
      r(1)=sqrt( (c(1,1)-c(1,2))**2+(c(2,1)-c(2,2))**2
     1                             +(c(3,1)-c(3,2))**2 )
      r(2)=sqrt( (c(1,2)-c(1,3))**2+(c(2,2)-c(2,3))**2
     1                             +(c(3,2)-c(3,3))**2 )
      r(3)=sqrt( (c(1,1)-c(1,3))**2+(c(2,1)-c(2,3))**2
     1                             +(c(3,1)-c(3,3))**2 )

C the dot product gives the colinearity chack, if dot is zero (i.e. small)
C  dot=r1*r3 cosa
C if a=180 or 0 degrees, cosa== -1 or 1
      dot =      (c(1,1)-c(1,2))*(c(1,1)-c(1,3))+
     1           (c(2,1)-c(2,2))*(c(2,1)-c(2,3))+
     2           (c(3,1)-c(3,2))*(c(3,1)-c(3,3))
      cosa=dot/(r(1)*r(3))

      cosa=abs(cosa)
      coll=abs(1.0-cosa).lt.small
              
      if(out1) then
         write(*,*) 'coll=',coll 
         write(*,*) 'd     = ',d1,d3
         write(*,*) 'be    = ',be1,be3
         write(*,*)  're     = ',re1,re3
      endif
      if(coll) then
         do i=1,3
            e1(i) = d1(i)*(1.0 - exp(-be1(i)*(r(i)-re1(i))))**2-d1(i)
            e3(i) = d3(i)*(1.0 + exp(-be3(i)*(r(i)-re3(i))))**2-d3(i)
            vq(i) = (e1(i) + e3(i))/2.0
            vj(i) = (e1(i) - e3(i))/2.0
         enddo
         vleps = vq(1)+vq(2)+vq(3)  
     1         -sqrt(vj(1)**2+vj(2)**2+vj(3)**2
     2              -vj(1)*vj(2)-vj(2)*vj(3)-vj(3)*vj(1))

C------- convert the energy from kcal_mol (because de is inkcal_mol) to hartree
         eleps=vleps/tohartree
         if(out1) then
             write(6,*)'r      = ',r
             write(6,*)'vq     = ',vq
             write(6,*)'vj     = ',vj
             write(6,*)'eleps (hartree) = ',eleps
c            call fiflsh(iw)
         endif
      else
         eleps=0.0
      endif
c      leps_energy1=eleps

      return
      END

c $Id$
