      subroutine argos_forces(iwl,isl,xw,xwm,xs,fw,fs)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
c
      integer i,j,k
c
      call argos_fclass(iwl,isl,xw,xwm,xs,fw,fs)
c
      call timer_start(44)
      call timer_start(202)
      call timer_start(204)
      call ga_sync()
      call timer_stop(204)
      call timer_stop(202)
      call timer_stop(44)
c
      call timer_start(45)
      call argos_space_final(fw,fs,lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(45)
c
      if(ltwin) then
      do 1 j=1,3
      do 2 k=1,nwa
      do 3 i=1,nwmloc
      fw(i,j,k,1)=fw(i,j,k,1)+fw(i,j,k,2)
    3 continue
    2 continue
      do 4 i=1,nsaloc
      fs(i,j,1)=fs(i,j,1)+fs(i,j,2)
    4 continue
    1 continue
      endif
c
cx moved to argos_stdesc where the norm of forces is needed
cx      call timer_start(47)
cx      call argos_cafe_fnorm(iwl(1,lwdyn),fw,nwmloc,
cx     + isl(1,lsatt),isl(1,lsdyn),fs,nsaloc,fnorm,fmax)
cx      call timer_stop(47)
c
      return
      end
c $Id$
