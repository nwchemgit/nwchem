      logical function pre_wrtrst(lfnout,lfnrst,filrst,title,
     + npbtyp,nbxtyp,box,
     + iwmr,xw,vw,mwm,mwa,nwm,nwa,iwmrc,xwc,vwc,mwmc,nwmc,
     + isar,isgm,xs,vs,msa,nsa,msm,nsm,irrep,nxrep,nyrep,nzrep,
     + mseq,nseq,lseq,ips,istat,lfnmod,filmod,itopol)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      logical pre_ssize,pre_rstmod
      external pre_ssize,pre_rstmod
c
      integer lfnout,lfnrst
      integer mwm,mwa,nwm,nwa,msa,nsa,msm,nsm,mwmc,nwmc,mseq,nseq
      integer npbtyp,nbxtyp,irrep,nxrep,nyrep,nzrep
      character*255 filrst,filmod
      character*80 title(2,3)
      real*8 xw(3,mwa,mwm),xs(3,msa),vw(3,mwa,mwm),vs(3,msa),box(3)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc)
      integer iwmrc(mwmc),iwmr(mwm),isar(msa),isgm(msa),lseq(mseq)
      integer ips(msa),istat(msa),lfnmod,itopol
c
      integer i,j,k,length,ixrep,iyrep,izrep,nrep,idhop
      real*8 orep(3),rsgm
      character*10 rstdat,rsttim
c
c     process modifications
c
      if(.not.pre_rstmod(lseq,mseq,nseq,lfnmod,filmod,lfnout))
     + call md_abort('pre_rstmod failed',9999)
c
      if(itopol.eq.3) then
      do 111 i=1,nsa
      idhop=2**(lseq(isgm(i))-1)
      if(iand(istat(i),idhop).eq.0) ips(i)=-ips(i)
c      write(*,'(a,5i5)') 'wrtrst ',i,isgm(i),lseq(isgm(i)),istat(i),
c     + ips(i)
  111 continue
      else
      do 112 i=1,nsa
      ips(i)=0
  112 continue
      endif
c
      if(.not.pre_ssize(isgm,xs,msa,nsa,rsgm))
     + call md_abort('pre_ssize failed',9999)
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
      if(npbtyp.eq.0) then
      if(box(1).lt.1.0e-6) box(1)=max(box(1),box(2),box(3))
      if(box(2).lt.1.0e-6) box(2)=max(box(1),box(2),box(3))
      if(box(3).lt.1.0e-6) box(3)=max(box(1),box(2),box(3))
      endif
c
      if(irrep.eq.0) then
      orep(1)=0.0d0
      orep(2)=0.0d0
      orep(3)=0.0d0
      nxrep=1
      nyrep=1
      nzrep=1
      nrep=1
      else
      nrep=nxrep*nyrep*nzrep
      endif
c
      write(lfnrst,1001) (title(2,i),i=1,3),4.6,rstdat,rsttim,0,.false.
 1001 format(a,/,a,/,a,/,f12.6,2a10,i5,4x,l1)
      write(lfnrst,1002) npbtyp,nbxtyp,rsgm,
     + dble(nxrep)*box(1),0.0d0,0.0d0,
     + 0.0d0,dble(nyrep)*box(2),0.0d0,
     + 0.0d0,0.0d0,dble(nzrep)*box(3)
 1002 format(2i5,f12.6,/,(3f12.6))
      write(lfnrst,1003) 0.0d0
 1003 format(e12.5)
      write(lfnrst,1004) 0.0d0,0.0d0,0.0d0
 1004 format(3f12.6)
c
c     dimensions on the restart file
c
c     1 i10 nwm   number of solvent molecules
c     2 i10 nwa   number of atoms per solvent molecule
c     3 i10 nsm   number of solute molecules
c     4 i10 nsa   number of solute atoms
c     5 i10 nwmc  number of crystal solvent molecules
c     6 i10 nsf   number of solute fractions
c     7 i10 nss   number of solute segments
c     8 i5  nprev number of processors used in previous job
c     9 i5  noe   number of noe constraints
c
c
      write(lfnrst,1005) nrep*(nwmc+nwm),nwa,
     + nrep*nsm,nrep*nsa,nrep*nwmc,0,nrep*nseq,0,0
 1005 format(7i10,2i5)
c
c     check wich version
c
      if(nwm+nwmc.gt.0) then
      if(nwmc.gt.0) then
      do 101 izrep=1,nzrep
      do 102 iyrep=1,nyrep
      do 103 ixrep=1,nxrep
      if(irrep.gt.0) then
      orep(1)=0.5d0*dble(2*ixrep-nxrep-1)*box(1)
      orep(2)=0.5d0*dble(2*iyrep-nyrep-1)*box(2)
      orep(3)=0.5d0*dble(2*izrep-nzrep-1)*box(3)
      endif
      do 1 i=1,nwmc
      write(lfnrst,1006) ((xwc(k,j,i)+orep(k),k=1,3),
     + (vwc(k,j,i),k=1,3),j=1,nwa)
 1006 format(2x,6f13.8)
      write(lfnrst,1007) iwmrc(i),(0.0d0,k=1,3)
 1007 format(i1,1x,3f13.8)
    1 continue
  103 continue
  102 continue
  101 continue
      endif
      do 201 izrep=1,nzrep
      do 202 iyrep=1,nyrep
      do 203 ixrep=1,nxrep
      if(irrep.gt.0) then
      orep(1)=0.5d0*dble(2*ixrep-nxrep-1)*box(1)
      orep(2)=0.5d0*dble(2*iyrep-nyrep-1)*box(2)
      orep(3)=0.5d0*dble(2*izrep-nzrep-1)*box(3)
      endif
      do 2 i=1,nwm
      write(lfnrst,1006) ((xw(k,j,i)+orep(k),k=1,3),
     + (vw(k,j,i),k=1,3),j=1,nwa)
      write(lfnrst,1007) iwmr(i),(0.0d0,k=1,3)
    2 continue
  203 continue
  202 continue
  201 continue
      endif
      if(nsa.gt.0) then
      do 301 izrep=1,nzrep
      do 302 iyrep=1,nyrep
      do 303 ixrep=1,nxrep
      if(irrep.gt.0) then
      orep(1)=0.5d0*dble(2*ixrep-nxrep-1)*box(1)
      orep(2)=0.5d0*dble(2*iyrep-nyrep-1)*box(2)
      orep(3)=0.5d0*dble(2*izrep-nzrep-1)*box(3)
      endif
      write(lfnrst,1008)
     + (isar(i),(xs(k,i)+orep(k),k=1,3),(vs(k,i),k=1,3),ips(i),i=1,nsa)
 1008 format(i1,1x,6f13.8,i5)
  303 continue
  302 continue
  301 continue
      endif
c
      if(nsm.gt.0) then
      write(lfnrst,1011) ((0.0d0,k=1,3),i=1,nrep*nsm)
 1011 format(2x,3f13.8)
      endif
c
      if(nseq.gt.0) then
      write(lfnrst,1012) ((lseq(i),i=1,nseq),j=1,nrep)
 1012 format(20i3)
      endif
c
      close(unit=lfnrst)
c
      if(util_print('files',print_default)) then
      write(lfnout,2000) filrst(1:length)
 2000 format(/,' Created restart',t40,a,/)
      endif
c
      pre_wrtrst=.true.
      return
c
 9999 continue
      pre_wrtrst=.false.
      return
      end
      logical function pre_wtrst(lfnout,lfnrst,filrst,
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
      pre_wtrst=.true.
      return
c
 9999 continue
      pre_wtrst=.false.
      return
      end
      logical function pre_ssize(isgm,xs,msa,nsa,rsgm)
c
      implicit none
c
      integer msa,nsa
      integer isgm(msa)
      real*8 xs(3,msa)
      real*8 rsgm
c
      integer i,j
      real*8 d
c
      rsgm=0.0d0
c
      do 1 i=1,nsa-1
      do 2 j=i+1,nsa
      if(isgm(i).eq.isgm(j)) then
      d=(xs(1,i)-xs(1,j))**2+(xs(2,i)-xs(2,j))**2+(xs(3,i)-xs(3,j))**2
      rsgm=max(rsgm,d)
      endif
    2 continue
    1 continue
      if(rsgm.gt.0) rsgm=dsqrt(rsgm)
c
      pre_ssize=.true.
c
      return
      end
