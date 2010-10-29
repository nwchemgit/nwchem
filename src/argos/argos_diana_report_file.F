      subroutine argos_diana_report_file(opt,fil,lfn)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
c
      character*8 opt
      character*255 fil
      integer lfn
c
      integer i,j,num,numi,numj
      integer i_ii,l_ii,i_ic,l_ic
      integer i_ji,l_ji,i_jc,l_jc,i_jr,l_jr
c
      if(me.ne.0) return
c
      if(opt.eq.'local   '.and.index(fil,'.loc').gt.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1))
      rewind(lfn)
      read(lfn,1000) num
 1000 format(i7)
      read(lfn,1000) numi
      do 1 i=1,numi
      read(lfn,1000) j
    1 continue
      read(lfn,1000) numj
      if(.not.ma_push_get(mt_int,3*numi,'ii',l_ii,i_ii))
     + call md_abort('Could not allocate ii',0)
      if(.not.ma_push_get(mt_byte,12*numi,'ic',l_ic,i_ic))
     + call md_abort('Could not allocate ic',0)
      if(.not.ma_push_get(mt_int,4*numj,'ji',l_ji,i_ji))
     + call md_abort('Could not allocate ji',0)
      if(.not.ma_push_get(mt_byte,12*numj,'jc',l_jc,i_jc))
     + call md_abort('Could not allocate jc',0)
      if(.not.ma_push_get(mt_dbl,numj,'jr',l_jr,i_jr))
     + call md_abort('Could not allocate jr',0)
      call argos_diana_report_loc(lfn,numi,numj,int_mb(i_ii),
     + int_mb(i_ji),
     + byte_mb(i_ic),byte_mb(i_jc),dbl_mb(i_jr))
      if(.not.ma_pop_stack(l_jr))
     + call md_abort('Failed to deallocate jr',0)
      if(.not.ma_pop_stack(l_jc))
     + call md_abort('Failed to deallocate jc',0)
      if(.not.ma_pop_stack(l_ji))
     + call md_abort('Failed to deallocate ji',0)
      if(.not.ma_pop_stack(l_ic))
     + call md_abort('Failed to deallocate ic',0)
      if(.not.ma_pop_stack(l_ii))
     + call md_abort('Failed to deallocate ii',0)
      close(unit=lfn)
      endif
c
      return
      end
