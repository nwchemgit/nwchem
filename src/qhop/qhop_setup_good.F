      subroutine qhop_setup(mparms,nparms,matm,natm,mseq,nseq,
     + catm,latm,cseq,lseq,mbnd,nbnd,lbnd,rbnd,mang,nang,lang,rang)
c
c $Id$
c
      implicit none
c
c     input variables
c     ---------------
c
c     matm            : dimension of atom arrays
c     natm            : number of atoms
c     mseq            : dimension of residue arrays
c     nseq            : number of residues
c     nparms          : number of parameter sets
c     mparms          : number of parameter sets plus 1
c     cseq            : residue names
c     catm(1,1:natm)  : atom names
c     catm(k,1:natm)  : atom types in parameter set k-1
c     latm(5,1:natm)  : residue number
c     latm(10,1:natm) : + heavy atom attached to protonatable hydrogen
c                       - protonatable hydrogen, relative index to heavy atom
c     lseq(5,1:nseq)  : number of protonation states of residue
c     lseq(6,1:nseq)  : current protonation state of residue
c
c     mbnd            : dimension of bond index array
c     nbnd            : number of bonds
c     lbnd(1,1:nbnd)  : atom index i
c     lbnd(2,1:nbnd)  : atom index j
c
      integer matm,natm,mseq,nseq,nparms,mparms
      character*10 cseq(mseq)
      character*6 catm(mparms,matm)
      integer latm(10,matm),lseq(6,mseq)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(nparms,2,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(nparms,2,mang)
c
      integer i,j,k,l,m,n
      integer mxseq,maxhv,maxpar,mxar
      parameter(mxseq=30)
      parameter(maxhv=4) 
      parameter(maxpar=10000)
      parameter(mxar=mxseq*30)
      integer numseq,numpar,mnprot,mnhv,nuseq,larat,nbnat
      integer prtin(2),flseq,flhv,nseql,eof,ptat(3)
      integer ptseq(mseq,2),numhv(mxseq),ihv,numprot(mxseq),
     + ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
      integer shbit(5),nbseq,nbhv,nbpt
      real*8 deq(mxar),aneq(mxar)
      character*10 cuseq(mxseq),cseqin(2)
      character*6 chv(mxseq,maxhv),chvin(2)  
      character*80 fileinp,fileout    
c
      real*8 e12fxy(maxpar,3),defe12(maxpar,7),zpef(maxpar,6),
     + tdsgl(maxpar,6),tunnel(maxpar,7),tstval(maxpar,5) 
c
c
      fileinp='qhop.par'
      open(unit=11,file=fileinp(1:index(fileinp,' ')-1),err=112)
      fileout='qhopcmp.par'
      open(unit=12,file=fileout(1:index(fileout,' ')-1),err=113)      
c
c     Initializing variables
c
      mnprot=0
      mnhv=0
      numpar=0
      numseq=0	
      nuseq=0
      flhv=0
      flseq=0
      ptat(1)=0
      larat=0
      do 31 i=1,mxar
      do j=1,4
      arat(i,j)=0
      enddo
      deq(i)=0
      aneq(i)=0
   31 continue
      do 12 i=1,mxseq
      ptseq(i,1)=0
      ptseq(i,2)=0
      ptarat(i)=0
   12 continue
c
c     Setting pointers to static array cuseq(numseq) and calculating max
c     heavy atom #  and prot states  
c       
      do i=1,nbnd
      write(*,101)i,(lbnd(j,i),j=1,2),rbnd(1,1,i)
  101 format(3i5,f12.6)
      enddo
c
      do i=1,nang
      write(*,102)i,(lang(j,i),j=1,3),rang(1,1,i)
  102 format(4i5,f12.6)
      enddo
c	
      do 1 i=1,natm 
c
      write(*,1000) i,catm(1,i),latm(5,i),cseq(latm(5,i)),latm(10,i),
     + lseq(5,latm(5,i)),lseq(6,latm(5,i))
c
 1000 format(i5,1x,a6,1x,i5,1xa10,3i5)
c	
      if((i.gt.1.and.latm(5,i).ne.latm(5,i-1)).or.i.eq.1)then
      if(flhv.eq.1)then
      ptat(2)=i-1
      flhv=0
      endif
      ptat(3)=i
      endif
c

      if(latm(10,i).gt.0)then
      flseq=1
      j=1
      do while (j.le.numseq.and.flseq.ne.0)
      if(cuseq(j).eq.cseq(latm(5,i)))then
      flseq=0
      ptseq(latm(5,i),1)=j      
      ptseq(latm(5,i),2)=ptat(3)
      endif
      j=j+1
      enddo	
c
      if(flseq.ne.0)then           
c
      flhv=1	
      if(lseq(5,latm(5,i)).gt.mnprot)mnprot=lseq(5,latm(5,i))
      numseq=numseq+1
      nseql=latm(5,i)
      cuseq(numseq)=cseq(latm(5,i))
      ptseq(latm(5,i),1)=numseq
      ptseq(latm(5,i),2)=ptat(3)
      numhv(numseq)=1
      if(mnhv.eq.0)mnhv=1
      chv(numseq,1)=catm(1,i)
c
      if(ptat(1).gt.0)then
      ptarat(numseq-1)=larat+1
c      print 1212,'ptat',numseq,ptat(1),ptat(2)
c 1212 format(a5,3i5) 
      ihv=0
      do 13 k=ptat(1),ptat(2)
      if(latm(10,k).gt.0)ihv=ihv+1
      nbnat=0
      larat=larat+1
      do 14 j=1,nbnd
      do 28 m=1,2
      l=2**(2-m)
      if(lbnd(m,j).eq.k.and.(latm(10,k).lt.0.or.latm(10,lbnd(l,j))
     + .ge.0))then
      nbnat=nbnat+1
      arat(larat,nbnat)=lbnd(l,j)-k
c      print 1029,'check ',numseq-1,k,lbnd(l,j),arat(larat,nbnat),
c     + larat,nbnat
c 1029 format(a6,6i6)
      if(latm(10,k).lt.0)then 
      deq(larat)=rbnd(1,1,j)
      print 121,'larat ',numseq-1,k,arat(larat,1),
     + latm(10,k+arat(larat,1))
  121 format(a6,4i5)
      arat(larat,nbnat+1)=latm(10,k+arat(larat,1)) 
      arat(larat,nbnat+2)=ihv
c
      if(latm(10,k+1).ge.0)then
      do 25 n=1,-arat(larat,nbnat)
      arat(larat-n+1,nbnat+3)=-arat(larat,nbnat)
   25 continue
      endif
c
      if(arat(larat,nbnat+1).eq.1)then
      do 18 n=1,nang
      if(k.eq.50)
     + print 222,'angle ',k,lbnd(l,j),k+arat(larat,1)
     + +arat(larat+arat(larat,1),1)
  222 format(a6,1x,3i5)
      if((lang(1,n).eq.k+arat(larat,1)+arat(larat+arat(larat,1),1)
     + .and.lang(2,n).eq.lbnd(l,j).and.lang(3,n).eq.k).or.
     + (lang(3,n).eq.k+arat(larat,1)+arat(larat+arat(larat,1),1)
     + .and.lang(2,n).eq.lbnd(l,j).and.lang(1,n).eq.k))
     + aneq(larat)=rang(1,1,n) 
   18 continue
       
      endif
      endif    
      endif
   28 continue
   14 continue
   13 continue
      endif      
      ptat(1)=ptat(3)
c
      elseif(flhv.eq.1.and.latm(5,i).eq.nseql)then
      numhv(numseq)=numhv(numseq)+1
      if(numhv(numseq).gt.mnhv)mnhv=numhv(numseq)
      chv(numseq,numhv(numseq))=catm(1,i)
      else
      flhv=0
      endif
      endif
    1 continue
c
c     Writing array arat for last unique res   
c  
      if(ptat(1).ge.0)then
      ihv=0 
      ptarat(numseq)=larat+1
      do 15 k=1,ptat(2)-ptat(1)+1
      nbnat=0
      larat=larat+1
      do 16 j=1,nbnd
      do 29 m=1,2
      l=2**(m-1)
      if(lbnd(m,j).eq.k.and.(latm(10,k).lt.0.or.latm(10,lbnd(l,j))
     + .ge.0))then
      nbnat=nbnat+1
      arat(larat,nbnat)=lbnd(l,j)-k
      if(latm(10,k).lt.0)then
      deq(larat)=rbnd(1,1,j)
      arat(larat,nbnat+1)=latm(10,k+arat(larat,nbnat))
      arat(larat,nbnat+2)=ihv
   
      if(latm(10,k+1).ge.0)then
      do 23 n=1,-arat(larat,nbnat)
      arat(larat-n+1,nbnat+3)=-arat(larat,nbnat)
   23 continue
      endif

      if(arat(larat,nbnat+1).eq.1)then
      do 19 n=1,nang
      if((lang(1,n).eq.k+arat(larat,1)+arat(larat+arat(larat,1),1)
     + .and.lang(2,n).eq.lbnd(l,j).and.lang(3,n).eq.k).or.
     + (lang(3,n).eq.k+arat(larat,1)+arat(larat+arat(larat,1),1)
     + .and.lang(2,n).eq.lbnd(l,j).and.lang(1,n).eq.k))
     + aneq(larat)=rang(1,1,n) 
   19 continue
      endif
      endif       
      endif
   29 continue
   16 continue
   15 continue
      endif
c
c     Calculating number of bits for ptrs to parameter sets
c
      do while(2**nbpt.le.mnprot)
      nbpt=nbpt+1
      enddo
      do while(2**(nbhv).le.mnhv)
      nbhv=nbhv+1
      enddo
      do while(2**(nbseq).le.numseq)
      nbseq=nbseq+1
      enddo
c     
      shbit(1)=nbpt
      shbit(2)=shbit(1)+nbhv
      shbit(3)=shbit(2)+nbseq
      shbit(4)=shbit(3)+nbpt
      shbit(5)=shbit(4)+nbhv 
c      print 1011,(shbit(i),i=1,5),mnprot,mnhv,numseq
c 1011 format('shbit ',8i5)
c
c     reading in the relevant hopping parameters
c
      eof=1
      do while (eof.ge.0) 
c
      read(11,1001,end=8)cseqin(1),chvin(1),prtin(1),cseqin(2),chvin(2)
     + ,prtin(2)     
 1001 format(a10,a6,i3,10x,a10,a6,i3)
      if(cseqin(1).eq.'          ')goto 8
      read(11,1003) (e12fxy(numpar+1,i),i=1,3)
 1003 format(3f12.6)
      read(11,1004) (defe12(numpar+1,i),i=1,7)
 1004 format(3f12.6,f8.4,2f12.6,f12.8)
      read(11,1005) (zpef(numpar+1,i),i=1,6)
 1005 format(6f12.6)
      read(11,1006) (tdsgl(numpar+1,i),i=1,6)
 1006 format(6f12.6)
      read(11,1007) (tunnel(numpar+1,i),i=1,7)
 1007 format(f8.2,f8.3,f10.5,e12.3,f10.5,2e12.3)
      read(11,1008) (tstval(numpar+1,i),i=1,5)
 1008 format(f10.2,2f10.4,f10.6,f10.4)
      do 4 i=1,numseq
      if(cuseq(i).eq.cseqin(1))then  
      do 5 j=1,numseq   
      if(cuseq(j).eq.cseqin(2))then
      do 6 k=1,numhv(i)
      if(chvin(1).eq.chv(i,k))then
      do 7 l=1,numhv(j)
      if(chvin(2).eq.chv(j,l))then
c     MANCA UN CONTROLLO SUI NUMERI CHE LEGGE SUI PROT STATES
      numpar=numpar+1
      ptpar(numpar)=lshift(i,shbit(5))+lshift(k,shbit(4))+
     + lshift(prtin(1),shbit(3))+lshift(j,shbit(2))+
     + lshift(l,shbit(1))+prtin(2) 
c      write(12,1030) cseqin(1),chvin(1),prtin(1),cseqin(2),chvin(2)
c     + ,prtin(2),ptpar(numpar)
c 1030 format(a10,a6,i3,10x,a10,a6,i3,i10)
c      write(12,1009) i,k,prtin(1),j,l,prtin(2)
c      write(12,1003)(e12fxy(numpar,m),m=1,3)     	
c      write(12,1004)(defe12(numpar,m),m=1,7)
c      write(12,1005)(zpef(numpar,m),m=1,6)
c      write(12,1006)(tdsgl(numpar,m),m=1,6)
c      write(12,1007)(tunnel(numpar,m),m=1,7)
c      write(12,1008)(tstval(numpar,m),m=1,5)
 1009 format(3i7,5x,3i7)
      else
c     ERR HEAVY ATOM2 NOT FIND 
      endif
    7 continue
      else 
c     ERR HEAVY ATOM1 NOT FIND   
      endif
    6 continue	 
      endif
    5 continue	
      endif
    4 continue
      enddo      
    8 continue
c    
      do 9 i=1,nseq
      if(ptseq(i,1).ne.0)nuseq=nuseq+1
    9 continue
c
      write(12,1014)numseq,larat,numpar
 1014 format(2i5,i10)
c
      do 11 i=1,nseq
      write(12,1010) i,cseq(i),(ptseq(i,j),j=1,2)
 1010 format(i5,1x,a10,2i5)
   11 continue
c
      do 21 i=1,numseq
      write(12,1020) i,cuseq(i),ptarat(i)
 1020 format(i5,1x,a10,1x,i5)
   21 continue
c
      do 22 i=1,larat
      write(12,1021) i,(arat(i,j),j=1,4), deq(i),aneq(i)  
 1021 format(i5,4i4,2f8.4)
   22 continue
c
      do 10 i=1,numpar 
      write(12,1015)ptpar(i)
      write(12,1003)(e12fxy(i,m),m=1,3)     	
      write(12,1004)(defe12(i,m),m=1,7)
      write(12,1005)(zpef(i,m),m=1,6)
      write(12,1006)(tdsgl(i,m),m=1,6)
      write(12,1007)(tunnel(i,m),m=1,7)
      write(12,1008)(tstval(i,m),m=1,5)
   10 continue  
 1015 format(i10)
c
      close(11)
      close(12)         
      return
c
  112 write(*,1030) 'err opening input file'
 1030 format(a40)
      return
  113 write(*,1030) 'err opening output file'
      return
      end
