      logical function argos_prepare_toprd(lfntop,filtop,lfnout,
     + num,amass,mat,cwa,qwa,mwa,nwa,iwat,isat,isgm,isml,isfr,
     + csa,qsa,msa,nsa,nsm,
     + idsb,cdsb,msb,nsb,qu,nqu,slvnam,mseq,nseq,lseq,ihop,istat)
c
c $Id$
c
      implicit none
c
      integer lfntop,lfnout,mwa,nwa,msa,nsa,nsm,msb,nsb,mat,mseq,nseq
      character*16 cwa(mwa),csa(msa)
      integer num(mat),isat(msa),isgm(msa),isml(msa),idsb(2,msb)
      integer iwat(mwa),isfr(msa),lseq(mseq),ihop(msa),istat(msa)
      real*8 amass(mat),qsa(msa),cdsb(msb),qwa(mwa)
      character*255 filtop
      character*10 slvmdl
      character*3 slvnam
      character*1 cdummy
c
      integer natl,nwb,nwh,nwd,nwo,nwt,nwn,nsh,nsd,nso,nst,nsx
      integer i,j,k,idummy,len,m,mf,nqu,npar,ival(25)
      real*8 qu(nqu),rval(25)
      character*3 string
c
      len=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:len),form='formatted',
     + status='old',err=9999)
c
      nsm=0
c
      do 1 i=1,4
      read(lfntop,1000)
 1000 format(a1)
    1 continue
      read(lfntop,1001) npar
      read(lfntop,1001) natl
      read(lfntop,1001) nqu
      read(lfntop,1001) nseq
 1001 format(i5)
      read(lfntop,1000)
      do 2 i=1,natl
      read(lfntop,1002) (ival(j),rval(j),j=1,npar)
 1002 format(5x,i5,7x,f12.6)
      amass(i)=rval(1)
      num(i)=ival(1)
    2 continue
      do 3 i=1,natl
      do 4 j=i,natl
      do 5 k=1,npar
      read(lfntop,1000) cdummy
    5 continue
    4 continue
    3 continue
      do 6 i=1,nqu
      read(lfntop,1333) (rval(j),j=1,npar)
 1333 format(5x,f12.6)
      qu(i)=rval(1)
    6 continue
      do 77 i=1,nseq
      read(lfntop,1334) lseq(i)
 1334 format(48x,i5)
      if(lseq(i).eq.0) lseq(i)=1
   77 continue
      read(lfntop,1003) nwa,nwb,nwh,nwd,nwo,nwt,nwn
      read(lfntop,1003) nsa,nsb,nsh,nsd,nso,nst,nsx
 1003 format(5i7,2i10)
      do 7 i=1,nwa
      read(lfntop,1336) cwa(i),iwat(i),idummy
 1336 format(a16,25x,2i5)
      qwa(i)=qu(idummy)
c      read(lfntop,1004) cwa(i)
c 1004 format(35x,a16)
c      read(lfntop,1010) qwa(i)
c 1010 format(f12.6)
    7 continue
      do 8 i=1,nwb*(npar+1)
      read(lfntop,1000) cdummy
    8 continue
      do 9 i=1,nwh*(npar+1)
      read(lfntop,1000) cdummy
    9 continue
      do 10 i=1,nwd*(npar+1)
      read(lfntop,1000) cdummy
   10 continue
      do 11 i=1,nwo*(npar+1)
      read(lfntop,1000) cdummy
   11 continue
      if(nwt.gt.0) then
      read(lfntop,1005) (idummy,i=1,nwt)
      read(lfntop,1005) (idummy,i=1,nwt)
 1005 format(11i7)
      endif
      if(nwn.gt.0) then
      read(lfntop,1005) (idummy,i=1,nwn)
      read(lfntop,1005) (idummy,i=1,nwn)
      endif
      if(nwa.gt.0) then
      read(lfntop,1006) slvmdl,string
 1006 format(a10,a3)
      do 14 i=1,npar
      read(lfntop,1000) cdummy
   14 continue
      endif
      if(string.ne.'   ') slvnam=string
      do 12 i=1,nsa
      read(lfntop,1339) csa(i),mf,m,isgm(i),isat(i),idummy,
     + ihop(i),istat(i)
 1339 format(a16,i3,2i7,14x,2i5,15x,i5,i3)
c 1339 format(a16,3i5,10x,2i5,15x,2i5)
      isml(i)=m
      isfr(i)=mf
      qsa(i)=qu(idummy)
c      read(lfntop,1007) isat(i),m,isgm(i),csa(i),n
c 1007 format(i7,i7,14x,i7,a16,7x,i5)
      if(nsm.lt.m) nsm=m
c      read(lfntop,1008) qsa(i)
c 1008 format(f12.6)
c      if(n.eq.3) then
c      read(lfntop,1000) cdummy
c      read(lfntop,1000) cdummy
c      endif
   12 continue
      do 13 i=1,nsb
      read(lfntop,1009) idsb(1,i),idsb(2,i)
 1009 format(2i7)
      read(lfntop,1010) (rval(j),j=1,npar)
 1010 format(f12.6)
      cdsb(i)=rval(1)
   13 continue
c
      close(unit=lfntop)
c
      argos_prepare_toprd=.true.
      return
c
 9998 continue
      close(unit=lfntop)
 9999 continue
      argos_prepare_toprd=.false.
      return
      end
