      logical function pre_toprd(lfntop,filtop,lfnout,
     + num,amass,mat,cwa,qwa,mwa,nwa,iwat,isat,isgm,isml,
     + csa,qsa,msa,nsa,nsm,
     + idsb,cdsb,msb,nsb,qu,nqu)
c
c $Id: pre_toprd.F,v 1.15 2000-06-27 18:07:44 d3j191 Exp $
c
      implicit none
c
      integer lfntop,lfnout,mwa,nwa,msa,nsa,nsm,msb,nsb,mat
      character*16 cwa(mwa),csa(msa)
      integer num(mat),isat(msa),isgm(msa),isml(msa),idsb(2,msb)
      integer iwat(mwa)
      real*8 amass(mat),qsa(msa),cdsb(msb),qwa(mwa)
      character*255 filtop
      character*10 slvnam
      character*1 cdummy
c
      integer natl,nwb,nwh,nwd,nwo,nwt,nwn,nsh,nsd,nso,nst,nsx
      integer i,j,idummy,len,m,nqu
      real*8 qu(nqu)
c
      len=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:len),form='formatted',
     + status='old',err=9999)
c
      nsm=0
c
      do 1 i=1,5
      read(lfntop,1000)
 1000 format(a1)
    1 continue
      read(lfntop,1001) natl
      read(lfntop,1001) nqu
 1001 format(i5)
      do 2 i=1,natl
      read(lfntop,1002) amass(i),num(i)
 1002 format(6x,f12.6,36x,i5)
    2 continue
      do 3 i=1,natl
      do 4 j=i,natl
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    4 continue
    3 continue
      do 335 i=1,nqu
      read(lfntop,1333) qu(i)
 1333 format(f12.6)
  335 continue
      read(lfntop,1003) nwa,nwb,nwh,nwd,nwo,nwt,nwn
      read(lfntop,1003) nsa,nsb,nsh,nsd,nso,nst,nsx
 1003 format(7i7)
      do 5 i=1,nwa
      read(lfntop,1336) cwa(i),iwat(i),idummy
 1336 format(a16,25x,2i5)
      qwa(i)=qu(idummy)
c      read(lfntop,1004) cwa(i)
c 1004 format(35x,a16)
c      read(lfntop,1010) qwa(i)
c 1010 format(f12.6)
    5 continue
      do 6 i=1,nwb
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    6 continue
      do 7 i=1,nwh
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    7 continue
      do 8 i=1,nwd
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    8 continue
      do 9 i=1,nwo
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    9 continue
      if(nwt.gt.0) then
      read(lfntop,1005) (idummy,i=1,nwt)
      read(lfntop,1005) (idummy,i=1,nwt)
 1005 format(11i7)
      endif
      if(nwn.gt.0) then
      read(lfntop,1005) (idummy,i=1,nwn)
      read(lfntop,1005) (idummy,i=1,nwn)
      endif
      if(nwa.gt.0) read(lfntop,1006) slvnam
 1006 format(a10)
      do 10 i=1,nsa
      read(lfntop,1339) csa(i),m,isgm(i),isat(i),idummy
 1339 format(a16,5x,2i5,10x,2i5)
      isml(i)=m
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
   10 continue
      do 11 i=1,nsb
      read(lfntop,1009) idsb(1,i),idsb(2,i),cdsb(i)
 1009 format(2i7,/,f12.6)
   11 continue
c
      close(unit=lfntop)
c
      pre_toprd=.true.
      return
c
 9998 continue
      close(unit=lfntop)
 9999 continue
      pre_toprd=.false.
      return
      end
