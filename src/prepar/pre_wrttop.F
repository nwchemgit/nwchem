      logical function pre_wrttop(ffield,lfntop,filtop,
     + releps,q14fac,lseq,cseq,mseq,nseq,
     + latt,lats,catt,patt,ratt,matt,natt,nats,latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp,l3rd,m3rd,n3rd,lexc,mexc,nexc)
c
      implicit none
c
      integer lfntop
      character*80 ffield
      character*255 filtop
      real*8 releps,q14fac
      integer mseq,nseq
      integer lseq(4,mseq)
      character*10 cseq(mseq)
      integer matt,natt,nats
      integer latt(3,matt),lats(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(8,matm)
      character*6 catm(4,matm)
      real*8 qatm(6,matm)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(6,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(6,mang)
      integer mdih,ndih
      integer ldih(9,mdih),kdih(6,3,mdih)
      real*8 rdih(6,6,mdih)
      integer mimp,nimp
      integer limp(9,mimp)
      real*8 rimp(6,mimp)
      integer m3rd,n3rd
      integer l3rd(2,m3rd)
      integer mexc,nexc
      integer lexc(2,mexc)
c
      integer i,j,k,l,length,nd
      character*10 topdat,toptim
c
      length=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:length),form='formatted',
     + status='new',err=9999)
c
      call arg_swatch(topdat,toptim)
c
      write(lfntop,1000) filtop(1:length),3.3,topdat,toptim,ffield(1:10)
 1000 format('TOPOLOGY FILE ',a,/,
     + 'GENERATED AUTOMICALLY BY NWCHEM PREPARE MODULE',//,
     + f12.6,3a10)
c
      
      write(lfntop,1004) nats,q14fac,releps,0.0d0,0
 1004 format(i5,3f12.6,i5)
      do 1 i=1,nats
      write(lfntop,1005) (catt(1,lats(j,i)),ratt(lats(j,i)),j=1,3),
     + (latt(1,lats(j,i)),j=1,3)
 1005 format(3(a6,f12.6),3i5)
    1 continue
      do 2 i=1,nats
      do 3 j=i,nats
      write(lfntop,1006) ((patt(1,k,lats(l,i),lats(l,j)),l=1,3),k=1,2)
      write(lfntop,1006) ((patt(2,k,lats(l,i),lats(l,j)),l=1,3),k=1,2)
 1006 format(6(1pe12.5))
    3 continue
    2 continue
c
      nd=0
      do 30 i=1,ndih
      nd=nd+max(ldih(7,i),ldih(8,i),ldih(9,i))
   30 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1009) natm,nbnd,nang,nd+nimp,0,n3rd,nexc
      else
      write(lfntop,1009) natm,nbnd,nang,nd,nimp,n3rd,nexc
 1009 format(5i7,2i10)
      endif
c
      do 10 i=1,natm
      write(lfntop,1010) latm(3,i),0,latm(1,i),latm(2,i),latm(5,i),
     + cseq(latm(5,i)),catm(1,i)(1:5),i,0,1
 1010 format(5i7,a10,a5,i10,2i5)
      write(lfntop,1011) (qatm(2*j-1,i),j=1,3),(qatm(2*j,i),j=1,3)
 1011 format(6f12.6)
   10 continue
c
      do 11 i=1,nbnd
      write(lfntop,1012) lbnd(1,i),lbnd(2,i),0,i
 1012 format(4i7)
      write(lfntop,1013) (rbnd(j,i),j=1,6)
 1013 format(3(0pf12.6,1pe12.5))
   11 continue
c
      do 12 i=1,nang
      write(lfntop,1014) lang(1,i),lang(2,i),lang(3,i),0,i
 1014 format(5i7)
      write(lfntop,1015) (rang(j,i),j=1,6)
 1015 format(3(0pf12.6,1pe12.5))
   12 continue
c
      do 13 i=1,ndih
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 14 k=1,nd
      write(lfntop,1016) (ldih(j,i),j=1,4),0,i
 1016 format(6i7)
      write(lfntop,1017)
     +  (kdih(k,j,i),rdih(k,2*j-1,i),rdih(k,2*j,i),j=1,3)
 1017 format(3(i3,0pf12.6,1pe12.5))
   14 continue
   13 continue
c
      do 15 i=1,nimp
      write(lfntop,1018) (limp(j,i),j=1,4),0,i
 1018 format(6i7)
      write(lfntop,1019) (0,rimp(2*j-1,i),rimp(2*j,i),j=1,3)
 1019 format(3(i3,0pf12.6,1pe12.5))
   15 continue
c
      if(n3rd.gt.0) then
      write(lfntop,1020) (l3rd(1,i),i=1,n3rd)
      write(lfntop,1020) (l3rd(2,i),i=1,n3rd)
 1020 format(11i7)
      endif
c
      if(nexc.gt.0) then
      write(lfntop,1021) (lexc(1,i),i=1,nexc)
      write(lfntop,1021) (lexc(2,i),i=1,nexc)
 1021 format(11i7)
      endif
c
      close(unit=lfntop)
c
c     write(lfntop,1007) natmw,nbndw,nangw,ndihw,nimpw
c      write(lfntop,1008) natm,nbnd,nang,ndih,nimp
c 1008 format(5i7)
c
      pre_wrttop=.true.
      return
 9999 continue
      pre_wrttop=.false.
      return
      end
