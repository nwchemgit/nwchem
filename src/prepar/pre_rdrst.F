      logical function pre_rrst(lfnout,lfnrst,filrst,
     + npbtyp,nbxtyp,box,
     + iwmr,xw,vw,mwm,mwa,nwm,nwa,iwmrc,xwc,vwc,mwmc,nwmc,
     + isar,xs,vs,msa,nsa,
     + inoe,dnoe,mnoe,noe,msm,nsm)
c
c $Id$
c
      implicit none
c
      integer lfnrst,lfnout
      integer mwm,mwa,nwm,nwa,msa,nsa,mnoe,noe,msm,nsm,mwmc,nwmc
      integer npbtyp,nbxtyp
      character*255 filrst
      real*8 xw(3,mwa,mwm),xs(3,msa),vw(3,mwa,mwm),vs(3,msa),box(3)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),dnoe(5,mnoe)
      integer iwmr(mwm),iwmrc(mwmc),isar(msa),inoe(2,mnoe)
c
      integer i,j,k,length,kinsep,nprev,nbxp,nbyp,nbzp
      character*10 rstdat,rsttim
      real*8 rdummy
c
#include "util.fh"
c
c     open the restart file
c     ---------------------
c
      length=index(filrst,' ')-1
      open(unit=lfnrst,file=filrst(1:length),form='formatted',
     + status='old',err=9999)
      rewind(lfnrst)
c
      call swatch(rstdat,rsttim)
c
      read(lfnrst,1001) rdummy,rstdat,rsttim
 1001 format(///,f12.6,2a10)
      read(lfnrst,1002) npbtyp,box,nbxtyp
 1002 format(i5,3f12.6,7x,i5)
      read(lfnrst,1003) kinsep
 1003 format(40x,i5)
      read(lfnrst,1004) rdummy,rdummy,rdummy
 1004 format(3f12.6)
      if(kinsep.eq.1) then
      read(lfnrst,1015) rdummy
      read(lfnrst,1015) rdummy
      read(lfnrst,1015) rdummy
 1015 format(e20.12)
      endif
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,nprev,noe
 1005 format(5i10,20x,2i5)
c
      if(nwa.gt.mwa) then
      write(lfnout,3001) nwa,mwa
 3001 format('Number of solvent atoms is ',i5,' while expecting ',i5)
      goto 9999
      endif
      if(nwm-nwmc.gt.mwm) then
      write(lfnout,3002) nwm-nwmc,mwm
 3002 format('Number of solvent molecules is ',i5,
     + ' while expecting ',i5)
      goto 9999
      endif
      if(noe.gt.mnoe) then
      write(lfnout,3003) noe,mnoe
 3003 format('Number of NOE constraints is ',i5,' while expecting ',i5)
      goto 9999
      endif
      if(nwmc.gt.mwmc) then
      write(lfnout,3004) nwmc,mwmc
 3004 format('Number of crystal solvent molecules is ',i5,' 
     + while expecting ',i5)
      goto 9999
      endif
      if(nsa.gt.msa) then
      write(lfnout,3005) nsa,msa
 3005 format('Number of solute atoms is ',i5,' while expecting ',i5)
      goto 9999
      endif
c
      if(nprev.gt.0) then
      read(lfnrst,1016) nbxp,nbyp,nbzp
 1016 format(3i5)
      read(lfnrst,1017) (rdummy,i=1,nbxp)
      read(lfnrst,1017) (rdummy,i=1,nbyp)
      read(lfnrst,1017) (rdummy,i=1,nbzp)
 1017 format(4e20.12)  
      endif
c
c     read coordinates
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      read(lfnrst,1006) (((xwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc),
     + (((xw(k,j,i),k=1,3),j=1,nwa),i=1,nwm-nwmc)
      else
      read(lfnrst,1006) (((xw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      endif
 1006 format(6f13.8)
      endif
c
      if(nsa.gt.0) then
      read(lfnrst,1007) ((xs(k,i),k=1,3),i=1,nsa)
 1007 format(6f13.8)
      endif
c
c     read velocities
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      read(lfnrst,1008) (((vwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc),
     + (((vw(k,j,i),k=1,3),j=1,nwa),i=nwmc+1,nwm)
      else
      read(lfnrst,1008) (((vw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      endif
 1008 format(6f13.8)
      endif
c
      if(nsa.gt.0) then
      read(lfnrst,1009) ((vs(k,i),k=1,3),i=1,nsa)
 1009 format(6f13.8)
      endif
c
c     read reference coordinates
c
      if(nwm.gt.0) then
      read(lfnrst,1010) ((rdummy,k=1,3),i=1,nwm)
 1010 format(6f13.8)
      endif
c
      if(nsm.gt.0) then
      read(lfnrst,1011) ((rdummy,k=1,3),i=1,nsm)
 1011 format(6f13.8)
      endif
c
c     read dynamics type
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      read(lfnrst,1012) (iwmrc(i),i=1,nwmc),(iwmr(i),i=nwmc+1,nwm)
      else
      read(lfnrst,1012) (iwmr(i),i=1,nwm)
      endif
 1012 format(80i1)
      endif
c
      if(nsa.gt.0) then
      read(lfnrst,1013) (isar(i),i=1,nsa)
 1013 format(80i1)
      endif
c
c     read solute NOE restraints
c
      if(noe.gt.0) then
      read(lfnrst,1014)
     + ((inoe(i,j),i=1,2),(dnoe(i,j),i=1,3),rdummy,
     + (dnoe(i,j),i=4,5),j=1,noe)
 1014 format(2i10,4f12.6,/,20x,2e12.5)
      endif
c
      close(unit=lfnrst)
c
      nwm=nwm-nwmc
c
      if(util_print('where',print_default)) then
      write(lfnout,3000) filrst(1:length)
 3000 format(' Read restart',t40,a,/)
      endif
c
      pre_rrst=.true.
      return
c
 9999 continue
      pre_rrst=.false.
      return
      end
      logical function pre_rdrst(lfnout,lfnrst,filrst,
     + npbtyp,nbxtyp,box,
     + iwmr,xw,vw,mwm,mwa,nwm,nwa,iwmrc,xwc,vwc,mwmc,nwmc,
     + isar,xs,vs,msa,nsa,
     + inoe,dnoe,mnoe,noe,msm,nsm,mseq,nseq,lseq,ips)
c
c $Id$
c
      implicit none
c
      integer lfnrst,lfnout
      integer mwm,mwa,nwm,nwa,msa,nsa,mnoe,noe,msm,nsm,mwmc,nwmc
      integer mseq,nseq
      integer npbtyp,nbxtyp
      character*255 filrst
      real*8 xw(3,mwa,mwm),xs(3,msa),vw(3,mwa,mwm),vs(3,msa),box(3)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),dnoe(5,mnoe)
      integer iwmr(mwm),iwmrc(mwmc),isar(msa),inoe(2,mnoe),lseq(mseq)
      integer ips(msa)
c
      integer i,j,k,length,nhist,nseqi
      character*10 rstdat,rsttim
      real*8 rdummy
      character*1 cdummy
      logical lforces
c
#include "util.fh"
c
c     open the restart file
c     ---------------------
c
      length=index(filrst,' ')-1
      open(unit=lfnrst,file=filrst(1:length),form='formatted',
     + status='old',err=9999)
      rewind(lfnrst)
c
      call swatch(rstdat,rsttim)
c
      read(lfnrst,1001) rdummy,rstdat,rsttim,nhist,lforces
 1001 format(///,f12.6,2a10,i5,4x,l1)
      if(nhist.gt.0) then
      do 1 i=1,nhist
      read(lfnrst,1000) cdummy
 1000 format(a1)
    1 continue
      endif
      read(lfnrst,1002) npbtyp,nbxtyp,box(1),rdummy,rdummy,rdummy,
     + box(2),rdummy,rdummy,rdummy,box(3)
 1002 format(2i5,/,(3f12.6))
      read(lfnrst,1003) rdummy
 1003 format(e12.5)
      read(lfnrst,1004) rdummy,rdummy,rdummy
 1004 format(3f12.6)
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,i,nseqi,i,noe
 1005 format(7i10,2i5)
c
      if(nwa.gt.mwa) then
      write(lfnout,3001) nwa,mwa
 3001 format('Number of solvent atoms is ',i5,' while expecting ',i5)
      goto 9999
      endif
      if(nwm-nwmc.gt.mwm) then
      write(lfnout,3002) nwm-nwmc,mwm
 3002 format('Number of solvent molecules is ',i5,
     + ' while expecting ',i5)
      goto 9999
      endif
      if(noe.gt.mnoe) then
      write(lfnout,3003) noe,mnoe
 3003 format('Number of NOE constraints is ',i5,' while expecting ',i5)
      goto 9999
      endif
      if(nwmc.gt.mwmc) then
      write(lfnout,3004) nwmc,mwmc
 3004 format('Number of crystal solvent molecules is ',i5,' 
     + while expecting ',i5)
      goto 9999
      endif
      if(nsa.gt.msa) then
      write(lfnout,3005) nsa,msa
 3005 format('Number of solute atoms is ',i5,' while expecting ',i5)
      goto 9999
      endif
      if(nseqi.ne.nseq) then
      write(lfnout,3055) nseqi,nseq
 3055 format('Number of segments is ',i5,' while expecting ',i5)
      goto 9999
      endif
c
c     read coordinates
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      do 2 i=1,nwmc
      read(lfnrst,1006) ((xwc(k,j,i),k=1,3),(vwc(k,j,i),k=1,3),j=1,nwa)
 1006 format(2x,6f13.8)
      if(lforces) read(lfnrst,1109) ((rdummy,j=1,3),k=1,nwa)
 1109 format(2x,6e13.6)
      read(lfnrst,1007) iwmrc(i),(rdummy,k=1,3)
 1007 format(i1,1x,3f13.8)
    2 continue
      endif
      do 3 i=1,nwm-nwmc
      read(lfnrst,1006) ((xw(k,j,i),k=1,3),(vw(k,j,i),k=1,3),j=1,nwa)
      if(lforces) read(lfnrst,1109) ((rdummy,j=1,3),k=1,nwa)
      read(lfnrst,1007) iwmr(i),(rdummy,k=1,3)
    3 continue
      endif
      if(nsa.gt.0) then
      do 4 i=1,nsa
      read(lfnrst,1008) isar(i),(xs(k,i),k=1,3),(vs(k,i),k=1,3),ips(i)
 1008 format(i1,1x,6f13.8,i5)
      if(lforces) read(lfnrst,1111) (rdummy,j=1,3)
 1111 format(2x,6e13.6) 
    4 continue
      endif
c
      if(nsm.gt.0) then
      read(lfnrst,1011) ((rdummy,k=1,3),i=1,nsm)
 1011 format(2x,3f13.8)
      endif
c
      if(nseq.gt.0) then
      read(lfnrst,1012) (lseq(i),i=1,nseq)
 1012 format(20i3)
      endif
c
c     read solute NOE restraints
c
c      if(noe.gt.0) then
c      read(lfnrst,1014) ((inoe(i,j),i=1,2),(dnoe(i,j),i=1,3),rdummy),
c     + (dnoe(i,j),i=4,5),j=1,noe)
c 1014 format(2i10,4f12.6,/,20x,2e12.5)
c      endif
c
      close(unit=lfnrst)
c
      nwm=nwm-nwmc
c
      if(util_print('where',print_default)) then
      write(lfnout,3000) filrst(1:length)
 3000 format(' Read restart',t40,a,/)
      endif
c
      pre_rdrst=.true.
      return
c
 9999 continue
      pre_rdrst=.false.
      return
      end
