      subroutine prp_start(ltwn,msfi,nsfi,mpei)
c
      implicit none
c
#include "prp_common.fh"
c
      integer msfi,nsfi,mpei
      logical ltwn
c
      ltwin=ltwn
      msf=msfi
      nsf=nsfi
      mpe=mpei
c
      return
      end
      subroutine prp_step(mdstep,temp,pres,eww,esw,ess)
c
      implicit none
c
#include "prp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer mdstep
      real*8 temp,pres
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2)
c
      integer i,j,k,me,np
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      if(ltwin) then
      do 1 i=1,mpe
      eww(i,1)=eww(i,1)+eww(i,2)
    1 continue
      do 2 i=1,msf
      do 3 j=1,mpe
      esw(i,j,1)=esw(i,j,1)+esw(i,j,2)
    3 continue
    2 continue
      do 4 i=1,msf
      do 5 j=1,msf
      do 6 k=1,mpe
      ess(i,j,k,1)=ess(i,j,k,1)+ess(i,j,k,2)
    6 continue
    5 continue
    4 continue
      endif
c
      if(np.gt.0) then
      call ga_dgop(mrg_d40,eww,mpe,'+')
      call ga_dgop(mrg_d41,esw,msf*mpe,'+')
      call ga_dgop(mrg_d42,ess,msf*msf*mpe,'+')
      endif
      if(me.ne.0) return
c
      write(*,'(i5,e12.5,9f12.3)') mdstep,pres,temp,(eww(i,1),i=1,8)
      do 7 i=1,nsf
      write(*,'(24x,i5,8f12.3)') i,(esw(i,j,1),j=1,8)
    7 continue
      do 8 i=1,nsf
      write(*,'(24x,i5,8f12.3)') i,(ess(i,i,j,1),j=1,8)
    8 continue
c
      return
      end
