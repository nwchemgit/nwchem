c $Id: cons.F,v 1.1 2004-01-28 01:30:59 marat Exp $
      subroutine cons_add_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      double precision energy
      integer geom
      integer nat
      logical status
      integer i
      integer inb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e
c
      call cons_init(rtdb)
c
      e=0.0d0
      do inb=1,n_bond
       iat=dbl_mb(i_bond+4*(inb-1))
       jat=dbl_mb(i_bond+4*(inb-1)+1)
       r0 =dbl_mb(i_bond+4*(inb-1)+2)
       k  =dbl_mb(i_bond+4*(inb-1)+3) 
       write(6,*) "adding constraints to atoms",iat,jat
c
       r=0.0d0
       do i=1,3
         ciat(i)=dbl_mb(i_c+(iat-1)*3+i-1)
         cjat(i)=dbl_mb(i_c+(jat-1)*3+i-1)
         r=r+(ciat(i)-cjat(i))**2
       end do
       r=sqrt(r)
c
       e=e+k*(r-r0)**2
      end do
      energy = energy + e
c
      call cons_end(rtdb)
      return
      end
c
c
      subroutine cons_add_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical status
      integer i
      integer inb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e,f(3)
c
      write(6,*) "cons_add enter"
c
      call cons_init(rtdb)
c
      do inb=1,n_bond
       iat=dbl_mb(i_bond+4*(inb-1))
       jat=dbl_mb(i_bond+4*(inb-1)+1)
       r0 =dbl_mb(i_bond+4*(inb-1)+2)
       k  =dbl_mb(i_bond+4*(inb-1)+3) 
       write(6,*) "adding constraints to atoms",iat,jat
c
       r=0.0d0
       do i=1,3
         ciat(i)=dbl_mb(i_c+(iat-1)*3+i-1)
         cjat(i)=dbl_mb(i_c+(jat-1)*3+i-1)
         r=r+(ciat(i)-cjat(i))**2
       end do
       r=sqrt(r)
c
       e=k*(r-r0)**2
       do i=1,3
         f(i)=k*(r-r0)*
     &     (ciat(i)-cjat(i))/r
       end do

       write(6,'("cons00  ",5F12.6)') r,e,f(1),f(2),f(3)     
       write(6,'("cons11 ",5F12.6)') r,energy,
     & gx((iat-1)*3+1),gx((iat-1)*3+2),gx((iat-1)*3+3)

       do i=1,3
         gx((iat-1)*3+i)=gx((iat-1)*3+i)+f(i)
         gx((jat-1)*3+i)=gx((jat-1)*3+i)-f(i)
         if(i.eq.1) write(6,*) "cons gradient",f
       end do
       energy = energy+e

       write(6,'("cons22 ",5F12.6)') r,energy,
     & gx((iat-1)*3+1),gx((iat-1)*3+2),gx((iat-1)*3+3)


      end do
c      write(6,*) "cons_add force"
c      do iat=1,nat
c        write(6,*) iat,i,gx(1,iat),
c     &    gx(2,iat),gx(3,iat)
c      end do 

      call cons_end(rtdb)

      return
      end  
c
c
c
      subroutine cons_init(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_data.fh"
c     
      integer rtdb
c
      integer geom
      integer nat
      integer i
      integer ma_type
      logical status
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_init: no geometry ',72, RTDB_ERR)
c
c     load harmonic constraints
c
      status = rtdb_ma_get(rtdb, 'cons:harm_bonds', ma_type, n_bond,
     &      h_bond)
      if (.not. status) call errquit(
     &    'cons_init: failed harmonic bonds ', 0, RTDB_ERR)
      n_bond=n_bond/4
      write(6,*) "number of harm bonds",n_bond
      status = ma_get_index(h_bond,i_bond)
      if (.not. status) call errquit(
     &    'cons_init: failed retrieving index ', 0, RTDB_ERR)

      do i=1,n_bond
       write(*,*) "harmonic bond",i
       write(*,*) "i,j,r0,k",dbl_mb(i_bond+4*(i-1)),
     &  dbl_mb(i_bond+4*(i-1)+1),dbl_mb(i_bond+4*(i-1)+2),
     &  dbl_mb(i_bond+4*(i-1)+3)
      end do
c
c     get cart coordinates
c     --------------------
      status=geom_ncent(geom,nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if ( .not. ma_alloc_get( MT_DBL, nat*3, 'cons coord',
     &      h_c, i_c) ) call errquit(
     &      'input_set: unable to allocate cons coord',
     &      nat*3, MA_ERR)
      status=geom_cart_coords_get(geom, dbl_mb(i_c))
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
 
      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)

      return
      end
c
c
c
      subroutine cons_end()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      logical status

      status = ma_free_heap(h_c)
      if (.not. status) call errquit(
     &    'cons_end: failed coord ', 0, MA_ERR)
      status = ma_free_heap(h_bond)
      if (.not. status) call errquit(
     &    'cons_end: failed harmonic bonds ', 0, MA_ERR)
      return
      end  

