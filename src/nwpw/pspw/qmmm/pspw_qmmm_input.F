      subroutine pspw_qmmm_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "errquit.fh"
c
c

      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
      character*30 dirs(num_dirs)
      data dirs / 'lj_ion_parameters:',
     >            'lj_ion_parameters',
     >            'lj_mm_parameters:',
     >            'polarization_mm_parameter:',
     >            'sr_basis:',
     >            'lr_cutoff:',
     >            'fragment',
     >            'end'/

      integer num_fdirs          ! No. of known directives
      parameter (num_fdirs = 5)
      character*30 fdirs(num_fdirs)
      data fdirs / 'size',
     >            'index_start',
     >            'shake',
     >            'switching',
     >            'end'/
c
      logical       value
      character*4   element
      character*80  rtdb_name
      integer       element_length,basis_length,tmp(2)
      integer       ni,nkf,nfrag,frag_size,index_start(2)
      real*8        se(2),alpha,Rin,Rout,sdist2(30),scal
      character*255 test
      integer ind               ! Index of matched directive
      integer j,jstart,jlast,jstride

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


      goto ( 100,100, 200, 300, 400, 500, 600,
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


      goto ( 610, 620,630,640,
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
 620  value = MA_push_get(mt_int,nw_max_atom,'indx_start',
     >                    index_start(2),index_start(1))
      if (.not. value) call errquit(
     > 'pspw_qmmm_input:failed allocating index_start',0,MA_ERR)

      ni = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start(1)+ni) = j
          ni = ni+1
         end do
      end do
c      do while (inp_i(index_start(ni+1)))
c        ni = ni + 1
c      end do
      rtdb_name = 'pspw_qmmm_frag_index_start:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,ni,int_mb(index_start(1)))
      if (.not. value) call errquit(
     > 'pspw_qmmm_input:failed writing index_start',0,RTDB_ERR)

      value = MA_pop_stack(index_start(2))
      if (.not. value) call errquit(
     > 'pspw_qmmm_input:failed deallocating index_start',0,MA_ERR)

      rtdb_name = 'pspw_qmmm_frag_nindex_start:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,ni)
      if (.not. value) call errquit(
     > 'pspw_qmmm_input:failed writing nindex_start',0, RTDB_ERR)

      nfrag = nfrag + ni
      rtdb_name = 'pspw_qmmm_nfrag'
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,nfrag)
      if (.not. value) call errquit(
     > 'pspw_qmmm_input:failed writing nfrag',0, RTDB_ERR)

      goto 601

c
c shake
c
 630  ni = 0
      if (inp_f(alpha)) then
         scal = 1.0d0
         sdist2(ni+1) = (scal*alpha)**2
         ni = ni + 1
      else
         call get_scalefrominput(scal)
      end if
      do while (inp_f(alpha))
        sdist2(ni+1) = (scal*alpha)**2
        ni = ni + 1
      end do
      rtdb_name = 'pspw_qmmm_frag_sdist:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_dbl,ni,sdist2)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing sdist',0, RTDB_ERR)

      rtdb_name = 'pspw_qmmm_frag_n_sdist:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_int,1,ni)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing n_sdist',0, RTDB_ERR)

      goto 601
c
c switching
c
 640  if (.not.inp_f(Rin))  Rin = 0.0d0
      if (.not.inp_f(Rout)) Rout = 0.0d0

      rtdb_name = 'pspw_qmmm_frag_switch_Rin:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,Rin)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing switch_Rin',0, RTDB_ERR)

      rtdb_name = 'pspw_qmmm_frag_switch_Rout:'//c_index_name(nkf)
      value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,Rout)
      if (.not. value)
     >  call errquit(
     > 'pspw_qmmm_input:failed writing switch_Rout',0, RTDB_ERR)

      goto 601


 699  continue
      goto 10

 9999 continue

      return
      end
