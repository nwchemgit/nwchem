      logical function pre_links(lfnout,ffield,lfnmod,filmod,
     + xatm,catm,latm,matm,lseq,cseq,mseq,nseq,
     + llnk,clnk,mlnk,nlnk,nlnkf)
c
c $Id: pre_links.F,v 1.13 2000-03-21 23:58:04 d3j191 Exp $
c
c     function to scan database for link segments in sequence
c
c     in    lseq(1,*) = sequence segment numbers
c           mseq      = dimension of the sequence list
c           nseq      = length of sequence list
c
c     out   cseq(2,*) = sequence segment names
c
      implicit none
c
#include "util.fh"
c
      logical pre_short,pre_protein,pre_nucleic
      external pre_short,pre_protein,pre_nucleic
c
      integer lfnout,matm,mseq,nseq,mlnk,nlnk,nlnkf,lfnmod
      character*255 filmod
      character*80 card
      real*8 xatm(3,matm)
      integer latm(5,matm)
      character*6 catm(3,matm)
      character*80 ffield
      integer lseq(4,mseq),llnk(4,mlnk)
      character*5 cseq(2,mseq)
      character*4 clnk(2,mlnk)
      integer ishort,jshort
      real*8 rshort
c
      integer i,j,iu,k,l,nlnko,nlnkp
      integer link1,link2,link3,mol
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'pre_links ')
      endif
c
c     get explicit links from input
c     -----------------------------
      open(unit=lfnmod,file=filmod(1:index(filmod,' ')-1),
     + form='formatted',status='old',err=999)
      rewind(unit=lfnmod)
   26 continue
      read(lfnmod,1000,end=99) card
 1000 format(a)
c
      if(card(1:4).eq.'link') then
      nlnk=nlnk+1
      read(card(6:80),1001)
     + llnk(1,nlnk),clnk(1,nlnk),llnk(2,nlnk),clnk(2,nlnk)
 1001 format(i5,a4,i5,a4)
      write(lfnout,1002)
     +  llnk(1,nlnk),clnk(1,nlnk),llnk(2,nlnk),clnk(2,nlnk)
 1002 format(' Input link   ',i5,':',a4,i5,':',a4)
      endif
      goto 26
   99 continue
      close(unit=lfnmod)
  999 continue
c
c     get atom index for links input from pdb
c     ---------------------------------------
c
      if(nlnk.gt.0) then
      do 21 i=1,nlnk
      do 22 j=1,nseq
      if(llnk(1,i).eq.lseq(1,j)) then
      llnk(1,i)=j
      goto 121
      endif
   22 continue
      call errquit('Illegal link on pdb',i)
  121 continue
      do 212 j=1,nseq
      if(llnk(2,i).eq.lseq(1,j)) then
      llnk(2,i)=j
      goto 21
      endif
  212 continue
      call errquit('Illegal link on pdb',i)
   21 continue
      do 23 i=1,nlnk
      do 24 j=lseq(3,llnk(1,i)),lseq(3,llnk(1,i)+1)-1
      if(clnk(1,i).eq.catm(2,j)(1:4)) llnk(3,i)=j
   24 continue
      do 25 j=lseq(3,llnk(2,i)),lseq(3,llnk(2,i)+1)-1
      if(clnk(2,i).eq.catm(2,j)(1:4)) llnk(4,i)=j
   25 continue
   23 continue
      endif
c
      nlnkf=0
      nlnkp=nlnk
      mol=0
c
c     check main chain links
c     ----------------------
c
      do 1 i=1,nseq-1
      mol=mol+1
      if(i.eq.1) lseq(4,i)=mol
      lseq(4,i+1)=mol+1
      if(pre_short(xatm,latm,matm,
     + lseq(3,i),lseq(3,i+1)-1,lseq(3,i+1),lseq(3,i+2)-1,
     + ishort,jshort,rshort)) then
      lseq(4,i+1)=lseq(4,i)
      if(latm(5,ishort).eq.0.or.latm(5,ishort).eq.2) then
      latm(5,ishort)=2
      else
      call errquit('inconsistent links',9999)
      endif
      if(latm(5,jshort).eq.0.or.latm(5,jshort).eq.1) then
      latm(5,jshort)=1
      else
      call errquit('inconsistent links',9999)
      endif
      if(.not.pre_protein(ffield,cseq(1,i)).or.
     + .not.pre_protein(ffield,cseq(1,i+1))) then
      nlnk=nlnk+1
      if(nlnk.gt.mlnk) call errquit('increase mlnk',9999)
      llnk(1,nlnk)=i
      llnk(2,nlnk)=i+1
      llnk(3,nlnk)=ishort
      llnk(4,nlnk)=jshort
      clnk(1,nlnk)=catm(2,ishort)(1:4)
      clnk(2,nlnk)=catm(2,jshort)(1:4)
      endif
      endif
    1 continue
c
c     check side chain links
c     ----------------------
c
      do 2 i=1,nseq-2
      do 3 j=i+2,nseq
      if(pre_short(xatm,latm,matm,
     + lseq(3,i),lseq(3,i+1)-1,lseq(3,j),lseq(3,j+1)-1,
     + ishort,jshort,rshort)) then
      if(latm(5,ishort).eq.0) latm(5,ishort)=3
      if(latm(5,jshort).eq.0) latm(5,jshort)=3
      nlnk=nlnk+1
      if(nlnk.gt.mlnk) call errquit('increase mlnk',9999)
      llnk(1,nlnk)=i
      llnk(2,nlnk)=j
      llnk(3,nlnk)=ishort
      llnk(4,nlnk)=jshort
      clnk(1,nlnk)=catm(2,ishort)(1:4)
      clnk(2,nlnk)=catm(2,jshort)(1:4)
      l=lseq(4,j)
      do 4 k=1,nseq
      if(lseq(4,k).eq.l) lseq(4,k)=lseq(4,i)
    4 continue
      endif
    3 continue
    2 continue
c
c     remove duplicate links
c     ----------------------
c
      if(nlnk.gt.0) then
      nlnko=nlnk
      nlnk=1
      if(nlnkp.ge.1) nlnkf=1
      do 5 i=2,nlnko
      do 6 j=1,nlnk
      if(llnk(1,i).eq.llnk(1,j).and.llnk(2,i).eq.llnk(2,j).and.
     + llnk(3,i).eq.llnk(3,j).and.llnk(4,i).eq.llnk(4,j).and.
     + clnk(1,i).eq.clnk(1,j).and.clnk(2,i).eq.clnk(2,j)) goto 5
      if(llnk(1,i).eq.llnk(2,j).and.llnk(2,i).eq.llnk(1,j).and.
     + llnk(3,i).eq.llnk(4,j).and.llnk(4,i).eq.llnk(3,j).and.
     + clnk(1,i).eq.clnk(2,j).and.clnk(2,i).eq.clnk(1,j)) goto 5
    6 continue
      nlnk=nlnk+1
      llnk(1,nlnk)=llnk(1,i)
      llnk(2,nlnk)=llnk(2,i)
      llnk(3,nlnk)=llnk(3,i)
      llnk(4,nlnk)=llnk(4,i)
      clnk(1,nlnk)=clnk(1,i)
      clnk(2,nlnk)=clnk(2,i)
      if(i.le.nlnkp) nlnkf=nlnk
    5 continue
      do 55 i=1,nlnk
      if(llnk(1,i).ne.llnk(2,i)) then
      if(latm(5,llnk(3,i)).eq.0) latm(5,llnk(3,i))=1
      if(latm(5,llnk(4,i)).eq.0) latm(5,llnk(4,i))=2
      endif
   55 continue
      endif
c
      do 7 i=1,nseq
      link1=0
      link2=0
      do 8 k=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,k).eq.1) link1=k
      if(latm(5,k).eq.2) link2=k
    8 continue
      do 9 j=1,nseq
      if(cseq(1,i).eq.cseq(1,j)) then
      do 10 l=lseq(3,j),lseq(3,j+1)-1
      if(latm(5,l).ge.3) then
      if(catm(1,l).eq.catm(1,link1)) latm(5,l)=1
      if(catm(1,l).eq.catm(1,link2)) latm(5,l)=2
      endif
   10 continue
      endif
    9 continue
    7 continue   
c     
c     change amino acid segment names based on links
c     ----------------------------------------------
c
      do 11 i=1,nseq
      if(pre_protein(ffield,cseq(1,i))) goto 14
      if(pre_nucleic(ffield,cseq(1,i))) goto 16
      do 15 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,j).eq.1) latm(5,j)=3
      if(latm(5,j).eq.2) latm(5,j)=3
   15 continue
      goto 11
   14 continue
      link1=0
      link2=0
      link3=0
      do 12 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,j).eq.1) link1=link1+1
      if(latm(5,j).eq.2) link2=link2+1
      if(latm(5,j).ge.3) link3=link3+1
   12 continue
      iu=index(cseq(1,i),'_')
      if(iu.eq.0) iu=index(cseq(1,i),' ')
      if(link1.gt.0) then
      if(link2.eq.0) cseq(1,i)(iu:iu+1)='_C'
      else
      if(link2.gt.0) then
      cseq(1,i)(iu:iu+1)='_N'
      else
      cseq(1,i)(iu:iu+1)='_M'
      endif
      endif
      goto 11
   16 continue
      link1=0
      link2=0
      link3=0
      do 17 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,j).eq.1) link1=link1+1
      if(latm(5,j).eq.2) link2=link2+1
      if(latm(5,j).ge.3) link3=link3+1
   17 continue
      iu=index(cseq(1,i),'_')
      if(iu.eq.0) iu=index(cseq(1,i),' ')
      if(link1.gt.0) then
      if(link2.eq.0) cseq(1,i)(iu:iu+1)='_3'
      else
      if(link2.gt.0) then
      cseq(1,i)(iu:iu+1)='_5'
      else
      cseq(1,i)(iu:iu+1)='_M'
      endif
      endif
   11 continue
c
c     adjust molecule number based on explicit links
c     ----------------------------------------------
c
      do 18 i=1,nlnk
      if(lseq(4,llnk(1,i)).ne.lseq(4,llnk(2,i))) then
      j=lseq(4,llnk(2,i))
      do 19 k=1,nseq
      if(lseq(4,k).eq.j) then
      lseq(4,k)=lseq(4,llnk(1,i))
      elseif(lseq(4,k).gt.j) then
      lseq(4,k)=lseq(4,k)-1
      endif
   19 continue
      endif
   18 continue
c
c     assign molecule number
c     ----------------------
c
      pre_links=.true.
      return
      end
