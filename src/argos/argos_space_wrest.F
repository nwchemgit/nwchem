      subroutine argos_space_wrest(lfn,lbbl,ltemp,mdim,boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfn,mdim
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
      real*8 boxsiz(maxbox,3)
c
      integer i,j,k
c
      if(me.eq.0) then
      write(lfn,1000)
 1000 format('restart space')
      write(lfn,1001) np,mbbl,npx,npy,npz,nbx,nby,nbz,np
 1001 format(9i7)
      write(lfn,1002) (boxsiz(i,1),i=1,nbx)
      write(lfn,1002) (boxsiz(i,2),i=1,nby)
      write(lfn,1002) (boxsiz(i,3),i=1,nbz)
 1002 format(4e20.12)
      endif
c
      do 1 i=1,np
      if(i.eq.me+1) then
      do 2 k=1,mbb2
      ltemp(1,k)=0
      do 3 j=1,nbbl
      ltemp(j+1,k)=lbbl(j,k)
    3 continue
    2 continue
      ltemp(1,1)=nbbl
      else
      do 4 k=1,mbb2
      do 5 j=1,mdim
      ltemp(j,k)=0
    5 continue
    4 continue
      endif
c
      call ga_igop(msp_09,ltemp,mdim*mbb2,'+')
c
      if(me.eq.0) then
      write(lfn,1003) i-1,ltemp(1,1)
 1003 format(2i7)
      do 6 j=1,ltemp(1,1)
      write(lfn,1004) (ltemp(j+1,k),k=1,4)
 1004 format(8i10)
    6 continue
      endif
c
    1 continue
c
      return
      end
c $Id$
