c
c $Id: brillouin_zone_input.F,v 1.1 2001-11-26 21:34:21 bylaska Exp $
c
      subroutine brillouin_zone_input(rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 3)
      character*22 dirs(num_dirs)
      character*255 test

      character*30 rtdb_name
      character*25 zone_name
      integer          num_kvectors,k,l
      double precision kvector(3)
      integer kvs(2),kvs_new(2)

      data dirs / 'zone_name:',
     >            'kvector:',
     >            'end'/

c     
c     **** intialize stuff ****
      num_kvectors = 0
      zone_name    = 'zone_default'
c
c
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'brillouine_zone_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'brillouine_zone_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'brillouine_zone_input: unknown directive', 0)


      goto ( 100,  200, 
     >      9999) ind
      call errquit(
     >      'brillouine_zone_input: unimplemented directive', ind)

c     
c  zone_name
c     
 100  if (.not. inp_a(zone_name))
     >     call errquit(
     >          'brillouine_zone_input: failed to read keyword', 0)
      goto 10


c
c kvector
c
 200  if (.not. inp_f(kvector(1)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0)
      if (.not. inp_f(kvector(2)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0)
      if (.not. inp_f(kvector(3)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0)

      num_kvectors = num_kvectors + 1
      if (.not. MA_alloc_get(mt_dbl,(3*num_kvectors), 
     >        'kvs_new',kvs_new(2),kvs_new(1)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 1', 0)

      if (num_kvectors.gt.1) then
        do k=1,(num_kvectors-1)
         dbl_mb(kvs_new(1) + 3*(k-1))   = dbl_mb(kvs(1) + 3*(k-1))
         dbl_mb(kvs_new(1) + 3*(k-1)+1) = dbl_mb(kvs(1) + 3*(k-1)+1)
         dbl_mb(kvs_new(1) + 3*(k-1)+2) = dbl_mb(kvs(1) + 3*(k-1)+2) 
        end do
        if (.not.MA_free_heap(kvs(2)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 2', 0)
      end if
      dbl_mb(kvs_new(1) + 3*(num_kvectors-1))   = kvector(1) 
      dbl_mb(kvs_new(1) + 3*(num_kvectors-1)+1) = kvector(2) 
      dbl_mb(kvs_new(1) + 3*(num_kvectors-1)+2) = kvector(3) 
      kvs(1) = kvs_new(1)
      kvs(2) = kvs_new(2)

      goto 10


*     ***** add brillouin_zone to rtdb ****
 9999 continue

      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0)

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (3*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0)


      if (.not.MA_free_heap(kvs(2)))
     >  call errquit(
     >       'brillouine_zone_input: heap failed 3', 0)
      return
      end

