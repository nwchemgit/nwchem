      subroutine argos_cafe_proper(vol,tmp,tmpw,tmps,prs,tsw,tss,psc,
     + ewkin,nwwl,nwws,nswl,nsws,nssl,nsss,nshw,nshs,ept2,ept3,
     + ep2m,ep3m,edistr,epomef,vir,prest,virt,ekct,ubiasi,pmfu)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 vol,tmp,tmpw,tmps,prs,tsw,tss,psc,ewkin,ept2,ept3
      real*8 edistr,epomef,vir,ubiasi,ep2m,ep3m
      real*8 prest(3,3),virt(3,3),ekct(3,3),pmfu(100)
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshw,nshs
c
      integer i,j
c
      vol=volume
      tmp=temp
      tmpw=tempw
      tmps=temps
      prs=pres
      tsw=scalew
      tss=scales
      psc=pscale
      vir=virial
c
      ewkin=ewk
c
      nwwl=llww
      nwws=lsww
      nswl=llsw
      nsws=lssw
      nssl=llss
      nsss=lsss
c
      nshw=niterw
      nshs=niters
c
      ept2=ep2(1)
      ept3=ep3(1)
      ep2m=ep2(3)
      ep3m=ep3(3)
c
      edistr=edrs
      epomef=epmf
c
      ubiasi=ubias
c
      do 1 i=1,3
      do 2 j=1,3
      prest(i,j)=prx(i,j)
      virt(i,j)=z(i,j)
      ekct(i,j)=ekc(i,j)
    2 continue
    1 continue
c
      if(npmfi.gt.0) then
      do 3 i=1,npmfi
      pmfu(i)=dbl_mb(i_upmf+i-1)
    3 continue
      endif
c
      return
      end
c $Id$
