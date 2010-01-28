      subroutine metadynamcis_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "errquit.fh"
c
c

      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 4)
      character*30 dirs(num_dirs)
      data dirs / 'bond',
     >            'angle',
     >            'dihedral',
     >            'end'/

c
      logical       value
      character*80  rtdb_name
      character*255 test,dir_name,sdir_name
      integer ind               ! Index of matched directive
      integer nka,nkb,nkd
      integer sindx(50)
      real*8  sdist(50)

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      nka   = 0
      nkb   = 0
      nkd   = 0

 10   if (.not. inp_read()) 
     >  call errquit(
     >  'metadyanmics_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'pspw_qmmm_input: failed to read keyword',0,INPUT_ERR)
      if (.not. inp_match(num_dirs,.false.,test,dirs,ind))
     >  call errquit(
     >  'metadynamics_input: unknown directive',0,INPUT_ERR)

      goto ( 100, 200, 300,
     >      9999) ind
      call errquit(
     >     'metadynamics_input: unimplemented directive',ind,INPUT_ERR)

c
c  bond
c
 100  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         nkb = nkb + 1
         rtdb_name = 'metadynamics_nbond'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkb)
         rtdb_name = 
     >       'metadynamics_bond_index:'//c_index_name(nkb)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,2,sindx)
         rtdb_name = 
     >       'metadynamics_bond_ksigma:'//c_index_name(nkb)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
         if (.not. value) call errquit(
     > 'metadynamics_input:failed writing bond',0, RTDB_ERR)
      end if
      end if
      end if
      end if

      goto 10

c
c  angle
c
 200  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         nka = nka + 1
         rtdb_name = 'metadynamics_nangle'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,nka)
         rtdb_name = 
     >       'metadynamics_angle_index:'//c_index_name(nka)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,3,sindx)
         rtdb_name = 
     >       'metadynamics_angle_ksigma:'//c_index_name(nka)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
         if (.not. value) call errquit(
     > 'metadynamics_input:failed writing angle',0, RTDB_ERR)
      end if
      end if
      end if
      end if
      end if

      goto 10

c
c  dihedral
c
 300  if (inp_i(sindx(1))) then
      if (inp_i(sindx(2))) then
      if (inp_i(sindx(3))) then
      if (inp_i(sindx(4))) then
      if (inp_f(sdist(1))) then
      if (inp_f(sdist(2))) then
         nkd = nkd + 1
         rtdb_name = 'metadynamics_ndihedral'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkd)
         rtdb_name =
     >       'metadynamics_dihedral_index:'//c_index_name(nkd)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,4,sindx)
         rtdb_name = 
     >       'metadynamics_dihedral_ksigma:'//c_index_name(nkd)
         value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,sdist)
         if (.not. value) call errquit(
     > 'metadynamics_input:failed writing dihedral',0, RTDB_ERR)
      end if
      end if
      end if
      end if
      end if
      end if

      goto 10


 9999 continue

      return
      end
