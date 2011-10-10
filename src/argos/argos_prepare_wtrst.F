      logical function argos_prepare_wtrst(lfnout,lfnrst,filrst,
     + npbtyp,nbxtyp,box,
     + iwmr,xw,vw,mwm,mwa,nwm,nwa,iwmrc,xwc,vwc,mwmc,nwmc,
     + isar,xs,vs,msa,nsa,
     + inoe,dnoe,mnoe,noe,msm,nsm)
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnrst
      integer mwm,mwa,nwm,nwa,msa,nsa,mnoe,noe,msm,nsm,mwmc,nwmc
      integer npbtyp,nbxtyp
      character*255 filrst
      real*8 xw(3,mwa,mwm),xs(3,msa),vw(3,mwa,mwm),vs(3,msa),box(3)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),dnoe(5,mnoe)
      integer iwmrc(mwmc),iwmr(mwm),isar(msa),inoe(2,mnoe)
c
      integer i,j,k,length
      character*10 rstdat,rsttim
c
c     open the restart file
c     ---------------------
c
      length=index(filrst,' ')-1
      open(unit=lfnrst,file=filrst(1:length),form='formatted',
     + status='unknown',err=9999)
c
      call swatch(rstdat,rsttim)
c
      write(lfnrst,1001) 3.1,rstdat,rsttim
 1001 format('This is an automatically generated restart file',/,
     + ' ',/,' ',/,f12.6,2a10)
      write(lfnrst,1002) npbtyp,box,nbxtyp
 1002 format(i5,3f12.6,7x,i5)
      write(lfnrst,1003) 0.0d0
 1003 format(e12.5)
      write(lfnrst,1004) 0.0d0,0.0d0,0.0d0
 1004 format(3f12.6)
      write(lfnrst,1005) nwmc+nwm,nwa,nsm,nsa,nwmc,0,0,0,noe
 1005 format(7i10,2i5)
c
c     check wich version
c
c     write coordinates
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      write(lfnrst,1006) (((xwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc),
     + (((xw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      else
      write(lfnrst,1006) (((xw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      endif
      elseif(nwmc.gt.0) then
      write(lfnrst,1006) (((xwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc)
 1006 format(6f13.8)
      endif
c
      if(nsa.gt.0) then
      write(lfnrst,1007) ((xs(k,i),k=1,3),i=1,nsa)
 1007 format(6f13.8)
      endif
c
c     write velocities
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      write(lfnrst,1008) (((vwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc),
     + (((vw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      else
      write(lfnrst,1008) (((vw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      endif
      elseif(nwmc.gt.0) then
      write(lfnrst,1008) (((vwc(k,j,i),k=1,3),j=1,nwa),i=1,nwmc)
 1008 format(6f13.8)
      endif
c
      if(nsa.gt.0) then
      write(lfnrst,1009) ((vs(k,i),k=1,3),i=1,nsa)
 1009 format(6f13.8)
      endif
c
c     write reference coordinates
c
      if(nwm+nwmc.gt.0) then
      write(lfnrst,1010) ((0.0d0,k=1,3),i=1,nwm+nwmc)
 1010 format(6f13.8)
      endif
c
      if(nsm.gt.0) then
      write(lfnrst,1011) ((0.0d0,k=1,3),i=1,nsm)
 1011 format(6f13.8)
      endif
c
c     write dynamics type
c
      if(nwm.gt.0) then
      if(nwmc.gt.0) then
      write(lfnrst,1012) (iwmrc(i),i=1,nwmc),(iwmr(i),i=1,nwm)
      else
      write(lfnrst,1012) (iwmr(i),i=1,nwm)
      endif
      elseif(nwmc.gt.0) then
      write(lfnrst,1012) (iwmrc(i),i=1,nwmc)
 1012 format(80i1)
      endif
c
      if(nsa.gt.0) then
      write(lfnrst,1013) (isar(i),i=1,nsa)
 1013 format(80i1)
      endif
c
c     write solute NOE restraints
c
      if(noe.gt.0) then
      write(lfnrst,1014)
     + ((inoe(i,j),i=1,2),(dnoe(i,j),i=1,3),1.0d0/(dnoe(2,j)**3),
     + (dnoe(i,j),i=4,5),j=1,noe)
 1014 format(2i10,4f12.6,/,20x,2e12.5)
      endif
c
      close(unit=lfnrst)
c
      if(util_print('files',print_default)) then
      write(lfnout,2000) filrst(1:length)
 2000 format(' Created old format restart',t40,a,/)
      endif
c
      argos_prepare_wtrst=.true.
      return
c
 9999 continue
      argos_prepare_wtrst=.false.
      return
      end
c $Id$
