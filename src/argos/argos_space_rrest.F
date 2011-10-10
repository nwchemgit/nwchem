      logical function argos_space_rrest(lfn,fil,lbbl,ltemp,mdim,boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfn,mdim
      character*255 fil
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
c
      integer i,j,k,node,npp,nbl,nbytes,nbxs
      integer npxp,npyp,npzp,nbxp,nbyp,nbzp,mbblx
      character*13 string
      real*8 boxsiz(maxbox,3)
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      npp=0
      read(lfn,1000,end=9997) string
 1000 format(a13)
      if(string.ne.'restart space') goto 1
      read(lfn,1001) npp,mbblx,npxp,npyp,npzp,nbxp,nbyp,nbzp,nbxs
 1001 format(9i7)
      if(mbblp.lt.mbblx) npp=0
      if(npxp.ne.npx) npp=0
      if(npyp.ne.npy) npp=0
      if(npzp.ne.npz) npp=0
      if(nbxp.ne.nbx) npp=0
      if(nbyp.ne.nby) npp=0
      if(nbzp.ne.nbz) npp=0
      if(nbxs.gt.0.and.npp.gt.0) then
      read(lfn,1002) (boxsiz(i,1),i=1,nbxp)
      read(lfn,1002) (boxsiz(i,2),i=1,nbyp)
      read(lfn,1002) (boxsiz(i,3),i=1,nbzp)
 1002 format(4e20.12)
      endif
 9997 continue
      do 2 i=1,mbb2
      ltemp(1,i)=0
    2 continue
      endif
c
      nbytes=ma_sizeof(mt_int,1,mt_byte)
      call ga_brdcst(msp_10,npp,nbytes,0)
c
      if(np.ne.npp) goto 9998
c
      nbytes=ma_sizeof(mt_dbl,3*maxbox,mt_byte)
      call ga_brdcst(msp_10,boxsiz,nbytes,0)
c
      nbytes=mdim*mbb2*ma_sizeof(mt_int,1,mt_byte)
c
      do 3 i=1,np
c
      if(me.eq.0) then
      read(lfn,1003) node,nbl
 1003 format(2i7)
      ltemp(1,1)=node
      ltemp(1,2)=nbl
      do 4 j=1,nbl
      read(lfn,1004) (ltemp(j+1,k),k=1,4)
 1004 format(8i10)
      if(mbb2.gt.4) then
      do 44 k=5,mbb2
      ltemp(j+1,k)=0
   44 continue
      endif
    4 continue
      endif
c
      call ga_brdcst(msp_11,ltemp,nbytes,0)
c
      if(ltemp(1,1).eq.me) then
      nbbl=ltemp(1,2)
      do 5 k=1,mbb2
      do 6 j=1,nbbl
      lbbl(j,k)=ltemp(j+1,k)
    6 continue
    5 continue
      endif
c
    3 continue
c
      nable=1
      argos_space_rrest=.true.
      return
 9998 continue
      nable=2
      argos_space_rrest=.false.
      return
 9999 continue
      call md_abort('Failed to open restart file',0)
      argos_space_rrest=.false.
      return
      end
c $Id$
