      subroutine pspw_qmmm_input(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
c

      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 7)
      character*30 dirs(num_dirs)
      data dirs / 'lj_ion_parameters:',
     >            'lj_mm_parameters:',
     >            'polarization_mm_parameter:',
     >            'sr_basis:',
     >            'lr_cutoff:',
     >            'fragment',
     >            'end'/

      integer num_fdirs          ! No. of known directives
      parameter (num_fdirs = 3)
      character*30 fdirs(num_fdirs)
      data fdirs / 'size',
     >            'index_start',
     >            'end'/
c
      logical       value
      character*4   element
      character*80  rtdb_name
      integer       element_length,basis_length,tmp(2)
      integer       ni,nkf,nfrag,frag_size,index_start(30)
      real*8        se(2),alpha
      character*255 test
      integer ind               ! Index of matched directive

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      nkf   = 0
      nfrag = 0

 10   if (.not. inp_read()) 
     >  call errquit(
     >  'pspw_qmmm_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'pspw_qmmm_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >  call errquit(
     >  'pspw_qmmm_input: unknown directive', 0, INPUT_ERR)


      goto ( 100, 200, 300, 400, 500, 600,
     >      9999) ind
      call errquit(
     >      'pspw_qmmm_input: unimplemented directive', ind, INPUT_ERR)

c
c  LJ Ion parameters - used by model potential and model charge
c
 100  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(se(1))
      value = value.and.inp_f(se(2))
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lj_ion_parameters:',0, INPUT_ERR)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_LJ_param_ion:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,se)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading LJ ion parameters:',0, INPUT_ERR)
      
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
     > 'pspw_qmmm_input:failed reading lj_mm_parameters:',0, INPUT_ERR)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_LJ_param_mm:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,2,se)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading LJ MM parameters:',0, RTDB_ERR)
      
      goto 10

c
c  Polarization parameter - used by model potential and model charge
c
 300  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading polarization_parameter:',0,
     &       INPUT_ERR)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_polarization_mm:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading polarization_parameter:',0,
     &       RTDB_ERR)
      
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
     > 'pspw_qmmm_input:failed reading SR_Basis:',0, MA_ERR)
      
      goto 10

c
c  lr_cutoff - used by model potential
c
 500  element = '    '
      value = inp_a(element)
      value = value.and.inp_f(alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lr_cutoff:',0, INPUT_ERR)

      element_length = index(element,' ') - 1
      rtdb_name = 'pspw_lr_cutoff:'//element(1:element_length)
      value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,alpha)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed reading lr_cutoff:',0, RTDB_ERR)
      
      goto 10

c
c  fragment
c
 600  nkf = nkf + 1
      rtdb_name = 'pspw_qmmm_nkfrag'
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,nkf)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing nkfrag',0, RTDB_ERR)

 601  if (.not. inp_read()) 
     >  call errquit('pspw_qmmm_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >call errquit('pspw_qmmm_input:failed reading keyword',0,INPUT_ERR)
      if (.not. inp_match(num_fdirs, .false.,test,fdirs,ind))
     >  call errquit('pspw_qmmm_input: unknown directive', 0, INPUT_ERR)


      goto ( 610, 620,
     >       699) ind
      call errquit('pspw_qmmm_input:directive not found',ind,INPUT_ERR)

c
c  size
c
 610  if (.not.inp_i(frag_size)) frag_size = 1
      rtdb_name = 'pspw_qmmm_frag_size:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,frag_size)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing frag_size',0, RTDB_ERR)

      goto 601
c
c index_start
c
 620  ni = 0
      do while (inp_i(index_start(ni+1)))
        ni = ni + 1
      end do
      rtdb_name = 'pspw_qmmm_frag_index_start:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,ni,index_start)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing index_start',0, RTDB_ERR)

      rtdb_name = 'pspw_qmmm_frag_nindex_start:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,ni)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing nindex_start',0, RTDB_ERR)

      nfrag = nfrag + ni
      rtdb_name = 'pspw_qmmm_nfrag'
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,nfrag)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing nfrag',0, RTDB_ERR)

      goto 601


 699  continue
      goto 10

 9999 continue

      return
      end
