c
c $Id$
c

      subroutine nwpw_set_library(rtdb)
      implicit none
      integer rtdb


#include "inp.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
      character*30 dirs(num_dirs)
      data dirs / 'pspw_library',
     >            'library',
     >            'paw_library',
     >            'cpi',
     >            'teter' /

      character*4  atom
      character*255 rtdb_name,library_type,library_name
      integer idum


 10   if (.not. inp_read())
     >     call errquit(
     >     'nwpw_set_library: inp_read failed', 0, INPUT_ERR)

      if (.not. inp_a(atom))
     >     call errquit(
     >     'nwpw_set_library: inp_a failed', 0, INPUT_ERR)

      if (inp_compare(.false.,atom,'end')) goto 100
      
      if (.not. inp_a(library_type))
     >     call errquit(
     >     'nwpw_set_library: inp_a failed', 0, INPUT_ERR)

      if (.not. inp_a(library_name))
     >     call errquit(
     >     'nwpw_set_library: inp_a failed', 0, INPUT_ERR)


      if(.not.inp_match(num_dirs,.false.,library_type,dirs,ind))
     > call errquit('nwpw_set_library: unknown library type', 
     >              0,INPUT_ERR)
c
c
      goto ( 20, 20, 30, 40, 50) ind
      call errquit('nwpw_set_library: unimplemented directive', ind,
     &       INPUT_ERR)


c
c pspw_library
c
  20  ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':pspw_library'
      ind = index(library_name,' ') - 1
      if (.not. btdb_cput(rtdb,rtdb_name,
     >                  1,library_name(1:ind)))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)

      rtdb_name = atom(1:ind)//':pspw_library_type'
      if (.not. btdb_put(rtdb,rtdb_name,mt_int,1,0))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)
      goto 10

c
c  paw_library
c
  30  ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':paw_library'
      ind = index(library_name,' ') - 1
      if (.not. btdb_cput(rtdb,rtdb_name,
     >                  1,library_name(1:ind)))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)
      goto 10

c
c CPI filename
c
  40  ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':pspw_library'
      ind = index(library_name,' ') - 1
      if (.not. btdb_cput(rtdb,rtdb_name,
     >                  1,library_name(1:ind)))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)

      idum = 2
      ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':pspw_library_type'
      if (.not. btdb_put(rtdb,rtdb_name,mt_int,1,idum))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)
      goto 10

c
c TETER filename
c
  50  ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':pspw_library'
      ind = index(library_name,' ') - 1
      if (.not. btdb_cput(rtdb,rtdb_name,
     >                  1,library_name(1:ind)))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)

      idum = 3
      ind = index(atom,' ') - 1
      rtdb_name = atom(1:ind)//':pspw_library_type'
      if (.not. btdb_put(rtdb,rtdb_name,mt_int,1,idum))
     > call errquit('nwpw_set_library: rtdb write error', 
     >              0,INPUT_ERR)
      goto 10



 100  continue

      return
      end
