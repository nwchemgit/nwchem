c
c $Id: brillouin_zone_input.F,v 1.5 2004-01-31 01:46:30 bylaska Exp $
c

*     *************************************************
*     *                                               *
*     *             brillouin_zone_input              *
*     *                                               *
*     *************************************************
      subroutine brillouin_zone_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
      character*22 dirs(num_dirs)
      character*255 test

      character*30 rtdb_name
      character*25 zone_name
      integer          num_kvectors,k,l
      double precision kvector(4)
      integer kvs(2),kvs_new(2)

      data dirs / 'zone_name:','zone_name',
     >            'kvector:','kvector',
     >            'end'/

c     
c     **** intialize stuff ****
      num_kvectors = 0
      zone_name    = 'zone_default'
c
c
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'brillouine_zone_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'brillouine_zone_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'brillouine_zone_input: unknown directive', 0,
     &       INPUT_ERR)


      goto ( 100,100,200,200, 
     >      9999) ind
      call errquit(
     >      'brillouine_zone_input: unimplemented directive', ind,
     &       INPUT_ERR)

c     
c  zone_name
c     
 100  if (.not. inp_a(zone_name))
     >     call errquit(
     >          'brillouine_zone_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c kvector
c
 200  if (.not. inp_f(kvector(1)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)
      if (.not. inp_f(kvector(2)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)
      if (.not. inp_f(kvector(3)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)

      if (.not. inp_f(kvector(4))) kvector(4) = -1.0d0
     

      num_kvectors = num_kvectors + 1
      if (.not. MA_alloc_get(mt_dbl,(4*num_kvectors), 
     >        'kvs_new',kvs_new(2),kvs_new(1)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 1', 0,
     &       INPUT_ERR)

      if (num_kvectors.gt.1) then
        do k=1,(num_kvectors-1)
         dbl_mb(kvs_new(1) + 4*(k-1))   = dbl_mb(kvs(1) + 4*(k-1))
         dbl_mb(kvs_new(1) + 4*(k-1)+1) = dbl_mb(kvs(1) + 4*(k-1)+1)
         dbl_mb(kvs_new(1) + 4*(k-1)+2) = dbl_mb(kvs(1) + 4*(k-1)+2) 
         dbl_mb(kvs_new(1) + 4*(k-1)+3) = dbl_mb(kvs(1) + 4*(k-1)+3) 
        end do
        if (.not.MA_free_heap(kvs(2)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 2', 0, MA_ERR)
      end if
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1))   = kvector(1) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+1) = kvector(2) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+2) = kvector(3) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+3) = kvector(4) 
      kvs(1) = kvs_new(1)
      kvs(2) = kvs_new(2)

      goto 10


*     ***** add brillouin_zone to rtdb ****
 9999 continue

      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)


      if (.not.MA_free_heap(kvs(2)))
     >  call errquit(
     >       'brillouine_zone_input: heap failed 3', 0, MA_ERR)
      return
      end


*     *************************************************
*     *                                               *
*     *             monkhost_pack_set                  *
*     *                                               *
*     *************************************************

      subroutine monkhorst_pack_set(rtdb,zone_name,nx,ny,nz)
      implicit none
      integer rtdb
      character*(*) zone_name
      integer nx,ny,nz

#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

      !**** local variables ****
      character*30 rtdb_name
      integer i1,i2,i3,k,l,nx,ny,nz,num_kvectors
      integer kvs(2)
      real*8  kx,ky,kz,xx,yy,zz,xxx,yyy,zzz


      num_kvectors = nx*ny*nz

      if(.not.MA_alloc_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >     call errquit(
     >          'monkhorst_pack_set: heap failed 1', 0,
     >       MA_ERR)

      xxx = 1.0d0/(2.0d0*nx)
      yyy = 1.0d0/(2.0d0*ny)
      zzz = 1.0d0/(2.0d0*nz)
      k = 0
      do i3=0,(nz-1)
      do i2=0,(ny-1)
      do i1=0,(nx-1)
         xx = 1.0 + 2*i1 - nx
         yy = 1.0 + 2*i2 - ny
         zz = 1.0 + 2*i3 - nz
 
         kx = xx*xxx
         ky = yy*yyy
         kz = zz*zzz
 
         dbl_mb(kvs(1) + 4*(k))   = kx
         dbl_mb(kvs(1) + 4*(k)+1) = ky
         dbl_mb(kvs(1) + 4*(k)+2) = kz
         dbl_mb(kvs(1) + 4*(k)+3) = -1.0d0

         k = k + 1
      end do
      end do
      end do


*     ***** add brillouin_zone to rtdb ****
      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'monkhorst_pack_set: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'monkhorst_pack_set: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)


      if (.not.MA_free_heap(kvs(2)))
     >  call errquit(
     >       'monkhorst_pack_set: heap failed 3', 0, MA_ERR)

      return 
      end
