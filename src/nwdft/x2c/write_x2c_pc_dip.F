c
      subroutine write_x2c_pc_dip(icomp,g_pcdip)
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "dra.fh"
#include "inp.fh"
c
      integer icomp
      integer g_pcdip

      character*(nw_max_path_len) fn_pcdip
      character*3 ch_comp
c
      logical dmat_to_file
      external dmat_to_file
c
c     prepare file name
      call util_file_name('pcdip',.false.,.false.,fn_pcdip)
      if (icomp.eq.1) ch_comp='.1'
      if (icomp.eq.2) ch_comp='.2'
      if (icomp.eq.3) ch_comp='.3'
      fn_pcdip = fn_pcdip(1:inp_strlen(fn_pcdip))//ch_comp   ! append component
c
c     resolve path name
      call util_file_name_resolve(fn_pcdip, .false.)
c
c     write ga to file
      if (.not. dmat_to_file(g_pcdip, fn_pcdip))
     &   call errquit('write_x2c_pc_dip: dmat_to_file error',0,
     &             UNKNOWN_ERR)
c
      return
      end
