      logical function pre_links(lfnout,ffield,xatm,catm,latm,matm,
     + lseq,cseq,mseq,nseq,
     + llnk,clnk,mlnk,nlnk)
c
c $Id: pre_links.F,v 1.6 1998-07-14 15:47:27 d3j191 Exp $
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
      logical pre_short,pre_protein
      external pre_short,pre_protein
c
      integer lfnout,matm,mseq,nseq,mlnk,nlnk
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
      integer i,j,iu,k,l
      integer link1,link2,link3,mol
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'pre_links ')
      endif
c
c     get atom index for links input from pdb
c     ---------------------------------------
c
      if(nlnk.gt.0) then
      do 21 j=1,nseq
      do 22 i=1,nlnk
      if(llnk(1,i).eq.lseq(1,j)) then
      llnk(1,i)=i
      goto 21
      endif
      if(llnk(2,i).eq.lseq(1,j)) then
      llnk(2,i)=i
      goto 21
      endif
   22 continue
      call errquit('Illegal link on pdb',9999)
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
      do 10 k=1,nseq
      if(lseq(4,k).eq.l) lseq(4,k)=lseq(4,i)
   10 continue
      endif
    3 continue
    2 continue
c
      do 4 i=1,nseq
      link1=0
      link2=0
      do 5 k=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,k).eq.1) link1=k
      if(latm(5,k).eq.2) link2=k
    5 continue
      do 6 j=1,nseq
      if(cseq(1,i).eq.cseq(1,j)) then
      do 7 l=lseq(3,j),lseq(3,j+1)-1
      if(latm(5,l).ge.3) then
      if(catm(1,l).eq.catm(1,link1)) latm(5,l)=1
      if(catm(1,l).eq.catm(1,link2)) latm(5,l)=2
      endif
    7 continue
      endif
    6 continue
    4 continue   
c     
c     change amino acid segment names based on links
c     ----------------------------------------------
c
      do 8 i=1,nseq
      if(pre_protein(ffield,cseq(1,i))) goto 11
      do 12 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,j).eq.1) latm(5,j)=3
      if(latm(5,j).eq.2) latm(5,j)=3
   12 continue
      goto 8
   11 continue
      link1=0
      link2=0
      link3=0
      do 9 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(5,j).eq.1) link1=link1+1
      if(latm(5,j).eq.2) link2=link2+1
      if(latm(5,j).ge.3) link3=link3+1
    9 continue
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
    8 continue
c
c     assign molecule number
c     ----------------------
c
      pre_links=.true.
      return
      end
