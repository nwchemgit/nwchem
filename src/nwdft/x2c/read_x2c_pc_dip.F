c
      subroutine read_x2c_pc_dip(icomp,g_pcdip)
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
#include "msgtypesf.h"
c
      integer icomp
      integer g_pcdip
c
      character*(nw_max_path_len) fn_pcdip
      character*3 ch_comp
c
      logical dmat_from_file
      external dmat_from_file
c
      integer inntsize,ok
c
c     prepare file name
      call util_file_name('pcdip',.false.,.false.,fn_pcdip)
      if (icomp.eq.1) ch_comp='.1'
      if (icomp.eq.2) ch_comp='.2'
      if (icomp.eq.3) ch_comp='.3'
      fn_pcdip = fn_pcdip(1:inp_strlen(fn_pcdip))//ch_comp   ! append component
c
c     resolve path
      call util_file_name_resolve(fn_pcdip, .false.)
c
c     read ga from file
      if (.not. dmat_from_file(g_pcdip,fn_pcdip))
     &     call errquit('read_x2c_pc_dip: dmat_from_file',0,
     &             UNKNOWN_ERR)
c
c     propagate status
      ok = 1
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      return
      end
