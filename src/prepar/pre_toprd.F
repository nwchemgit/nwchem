      logical function pre_toprd(lfntop,filtop,lfnout,lprint,
     + cwa,mwa,nwa,isat,isgm,csa,qsa,msa,nsa,nsm)
c
      implicit none
c
      integer lfntop,lfnout,lprint,mwa,nwa,msa,nsa,nsm
      character*16 cwa(mwa),csa(msa)
      integer isat(msa),isgm(msa)
      real*8 qsa(msa)
      character*255 filtop
      character*10 slvnam
      character*1 cdummy
c
      integer natl,nwb,nwh,nwd,nwo,nwt,nwn,nsb,nsh,nsd,nso,nst,nsx
      integer i,j,idummy,len,n,m
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
      read(lfntop,1001) natl
 1001 format(i5)
      do 2 i=1,natl
      read(lfntop,1000) cdummy
    2 continue
      do 3 i=1,natl
      do 4 j=i,natl
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    4 continue
    3 continue
      read(lfntop,1002) nwa,nwb,nwh,nwd,nwo,nwt,nwn
      read(lfntop,1002) nsa,nsb,nsh,nsd,nso,nst,nsx
 1002 format(7i7)
      do 5 i=1,nwa
      read(lfntop,1003) cwa(i)
 1003 format(35x,a16)
      read(lfntop,1000) cdummy
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
      read(lfntop,1004) (idummy,i=1,nwt)
      read(lfntop,1004) (idummy,i=1,nwt)
 1004 format(11i7)
      endif
      if(nwn.gt.0) then
      read(lfntop,1004) (idummy,i=1,nwn)
      read(lfntop,1004) (idummy,i=1,nwn)
      endif
      if(nwa.gt.0) read(lfntop,1005) slvnam
 1005 format(a10)
      do 10 i=1,nsa
      read(lfntop,1006) isat(i),m,isgm(i),csa(i),n
 1006 format(i7,i7,14x,i7,a16,10x,i5)
      if(nsm.lt.m) nsm=m
      read(lfntop,1007) qsa(i)
 1007 format(f12.6) 
      if(n.eq.3) then
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      endif
   10 continue
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
