      subroutine pspw_qmmm_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'lj_ion_parameters:',
     >            'lj_mm_parameters:',
     >            'polarization_mm_parameter:',
     >            'sr_basis:',
     >            'end'/
c
      logical       value
      character*4   element
      character*80  rtdb_name
      integer       element_length
      real*8        se(2),alpha


 10   if (.not. inp_read()) 
     >  call errquit(
     >  'pspw_qmmm_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >  call errquit(
     >  'pspw_qmmm_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >  call errquit(
     >  'pspw_qmmm_input: unknown directive', 0)


      goto ( 100, 200, 300, 400,
     >      9999) ind
      call errquit(
     >      'pspw_qmmm_input: unimplemented directive', ind)

c
c  LJ Ion parameters
c
 100  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(se(1))
      value = value.and.inp_f(se(2))
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lj_ion_parameters:',0)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_LJ_param_ion:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,se)
      
      goto 10

c
c  LJ MM parameters
c
 200  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(se(1))
      value = value.and.inp_f(se(2))
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lj_mm_parameters:',0)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_LJ_param_mm:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,se)
      
      goto 10

c
c  Polarization parameter
c
 300  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading polarization_parameter:',0)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_polarization_mm:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,alpha)
      
      goto 10

c
c  SR Basis
c
 400  continue
      goto 10

 9999 continue

      return
      end
