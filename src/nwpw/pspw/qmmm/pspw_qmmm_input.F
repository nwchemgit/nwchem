#
# $Id: pspw_qmmm_input.F,v 1.3 2001-10-26 04:44:16 bylaska Exp $
#

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
      parameter (num_dirs = 6)
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'lj_ion_parameters:',
     >            'lj_mm_parameters:',
     >            'polarization_mm_parameter:',
     >            'sr_basis:',
     >            'lr_cutoff:',
     >            'end'/
c
      logical       value
      character*4   element
      character*80  rtdb_name
      integer       element_length,basis_length,tmp(2)
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


      goto ( 100, 200, 300, 400, 500,
     >      9999) ind
      call errquit(
     >      'pspw_qmmm_input: unimplemented directive', ind)

c
c  LJ Ion parameters - used by model potential and model charge
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
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading LJ ion parameters:',0)
      
      goto 10

c
c  LJ MM parameters - used by model potential and model charge
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
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading LJ MM parameters:',0)
      
      goto 10

c
c  Polarization parameter - used by model potential and model charge
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
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading polarization_parameter:',0)
      
      goto 10

c
c  SR Basis - used by model potential
c
 400  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(se(1))
      value = value.and.inp_f(se(2))

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_SR_basis_length:'//element(1:element_length)
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,basis_length)) then
         basis_length = 0
      end if
      basis_length = basis_length+1
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,1,basis_length)


      value = value.and.
     >         MA_push_get(mt_dbl,(2*basis_length),'tmp',tmp(2),tmp(1))

      rtdb_name = 'pspw_SR_basis:'//element(1:element_length)
      if (basis_length.gt.1) then
        value = value.and.
     >          rtdb_get(rtdb,rtdb_name,
     >                 mt_dbl,(2*basis_length-2),dbl_mb(tmp(1)))
      end if
      dbl_mb(tmp(1)+2*basis_length-2) = se(1)
      dbl_mb(tmp(1)+2*basis_length-1) = se(2)
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,
     >                 mt_dbl,(2*basis_length),dbl_mb(tmp(1)))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading SR_Basis:',0)
      
      goto 10

c
c  lr_cutoff - used by model potential
c
 500  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lr_cutoff:',0)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_lr_cutoff:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lr_cutoff:',0)
      
      goto 10

 9999 continue

      return
      end
