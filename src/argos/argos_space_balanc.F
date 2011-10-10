      subroutine argos_space_balanc(stime,syntim,waltim,lpsyn)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      logical argos_space_diffbb
      external argos_space_diffbb
c
      real*8 stime,syntim,waltim
      logical lpsyn
c
      integer l_dlb,i_dlb,i_nod,l_nod,i_nen,l_nen,i_lst,l_lst
c
c      if(loadb.eq.0) return
c
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in argos_space_balanc at 1',me)
      if(.not.ma_push_get(mt_int,np,'nod',l_nod,i_nod))
     + call md_abort('Failed to allocate nod',me)
      if(.not.ma_push_get(mt_int,2*(np+1),'lst',l_lst,i_lst))
     + call md_abort('Failed to allocate lst',me)
      if(.not.ma_push_get(mt_int,mbbl*np,'nen',l_nen,i_nen))
     + call md_abort('Failed to allocate nen',me)
      if(.not.ma_push_get(mt_dbl,np+1,'dlb',l_dlb,i_dlb))
     + call md_abort('Failed to allocate dlb',me)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in argos_space_balanc at 2',me)
c
      call argos_space_dldbal(stime,syntim,waltim,int_mb(i_nod),
     + int_mb(i_lst),
     + int_mb(i_nen),int_mb(i_bb),dbl_mb(i_dlb),lpsyn,int_mb(i_bindex),
     + int_mb(i_buren),np+1)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in argos_space_balanc at 3',me)
c
      if(.not.ma_pop_stack(l_dlb))
     + call md_abort('Failed to deallocate dlb',me)
      if(.not.ma_pop_stack(l_nen))
     + call md_abort('Failed to deallocate nen',me)
      if(.not.ma_pop_stack(l_lst))
     + call md_abort('Failed to deallocate lst',me)
      if(.not.ma_pop_stack(l_nod))
     + call md_abort('Failed to deallocate nod',me)
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in argos_space_balanc at 4',me)
c
      if(argos_space_diffbb(dbl_mb(i_boxs),int_mb(i_rng))) then
      call argos_space_numbb(int_mb(i_iown),dbl_mb(i_boxs))
      call argos_space_listbb(int_mb(i_iown),dbl_mb(i_boxs),
     + int_mb(i_bb))
      endif
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('Problems in argos_space_balanc at 5',me)
c
      return
      end
c $Id$
