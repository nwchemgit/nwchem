      subroutine argos_prop_mcti(n,filnam)
c
      implicit none
c
#include "argos_prop.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer n,i_dec,l_dec,ibl,i,igp_handle
      character*80 filnam,string
      character*5 gid
c
      if(npg.gt.1) then
      if(me.eq.0) then
      close(unit=lfngib,status='keep')
      endif
c
c     now go global
c
      igp_handle=ga_pgroup_get_default()
      call ga_pgroup_set_default(ga_pgroup_get_world())
      call ga_sync()
      if(ga_nodeid().eq.0) then
c
      ibl=index(filnam,' ')-1
      string=filnam(1:ibl)//'.gib'
      open(unit=lfngib,file=string,status='unknown',form='formatted')
      rewind(lfngib)
      do 1 i=1,npg
      write(gid,3300) i-1
 3300 format(i5.5)
      string=filnam(1:ibl)//gid//'.gib'
      open(unit=75,file=string,status='old',form='formatted')
      rewind(unit=75)
    2 continue
      read(75,3301,end=99,err=99) string
      write(lfngib,3301) string
 3301 format(a)
      goto 2
   99 continue
      close(unit=75)
    1 continue
      close(unit=lfngib,status='keep')
c
      endif
c
c     go back to groups
c
      call ga_sync()
      call ga_pgroup_set_default(igp_handle)
c
      endif
c
      if(meg.eq.npg-1.or.npg.eq.1) then
c
      if(me.eq.0) then
      ibl=index(filnam,' ')-1
      string=filnam(1:ibl)//'.gib'
      open(unit=lfngib,file=string,status='old',form='formatted')
      rewind(lfngib)
      endif
c
      npgdec=n
c
      if(npgdec.gt.0) then
      if(.not.ma_push_get(mt_dbl,6*nsa,'dec',l_dec,i_dec))
     + call md_abort('Failed to allocate dec',0)
      else
      if(.not.ma_push_get(mt_dbl,1,'dec',l_dec,i_dec))
     + call md_abort('Failed to allocate dec',0)
      endif
c
      call argos_prop_mcti_s(dbl_mb(i_dfr),dbl_mb(i_dfrm),dbl_mb(i_dec))
c
      if(.not.ma_pop_stack(l_dec))
     + call md_abort('Failed to deallocate dec',0)
c
      endif
c      endif
c
      return
      end
c $Id$
