      logical function pre_wrttop(lfnout,ffield,lfntop,filtop,
     + releps,q14fac,lseq,cseq,mseq,nseq,mqu,
     + latt,lats,catt,patt,ratt,matt,natt,nats,latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp,l3rd,m3rd,n3rd,lexc,mexc,nexc,lnon,mnon,nnon,
     + natmt,nbndt,nangt,ndiht,nimpt,n3rdt,wcorr,slvnam)
c
c $Id: pre_wrttop.F,v 1.34 2001-07-06 03:45:14 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfntop
      character*80 ffield
      character*255 filtop
      real*8 releps,q14fac
      integer mseq,nseq,mqu
      integer lseq(4,mseq)
      character*10 cseq(mseq)
      integer matt,natt,nats
      integer latt(3,matt),lats(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(9,matm)
      character*6 catm(4,matm)
      real*8 qatm(6,mqu)
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
      integer mnon,nnon
      integer lnon(2,mnon)
      integer natmt,nbndt,nangt,ndiht,nimpt,n3rdt
      real*8 wcorr(3)
      character*3 slvnam
c
      integer i,j,k,l,length,nd,ifrac,imol,iq,isgm
      character*10 topdat,toptim
      real*8 qsum(3),qtot(3),qatom(6)
      logical lprt
c
      integer nqu,iql(3)
c
      nqu=0
c
c     generate list with unique charges and polarizabilities
c
      do 331 i=1,natmt
      do 332 j=1,nqu
      if(qatm(1,i).eq.qatm(1,j).and.qatm(2,i).eq.qatm(2,j).and.
     + qatm(3,i).eq.qatm(3,j).and.qatm(4,i).eq.qatm(4,j).and.
     + qatm(5,i).eq.qatm(5,j).and.qatm(6,i).eq.qatm(6,j)) then
      latm(7,i)=j
      latm(8,i)=j
      latm(9,i)=j
      goto 331
      endif
  332 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call errquit('Increase mqu',0)
      latm(7,i)=nqu
      do 333 j=1,6
      qatm(j,nqu)=qatm(j,i)
  333 continue
  331 continue
      nqu=nqu+1
      do 334 j=1,6
      qatm(j,nqu)=0.0d0
  334 continue
      do 335 i=1,natm
      do 336 j=1,3
      qatom(2*j-1)=qatm(2*j-1,latm(7,i))
      qatom(2*j)=qatm(2*j,latm(7,i))
      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'D') then
      qatom(2*j-1)=0.0d0
      qatom(2*j)=0.0d0
      endif
  336 continue
      do 337 j=1,nqu
      if(qatom(1).eq.qatm(1,j).and.qatom(2).eq.qatm(2,j).and.
     + qatom(3).eq.qatm(3,j).and.qatom(4).eq.qatm(4,j).and.
     + qatom(5).eq.qatm(5,j).and.qatom(6).eq.qatm(6,j)) then
      latm(7,i)=j
      latm(8,i)=j
      latm(9,i)=j
      goto 335
      endif
  337 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call errquit('Increase mqu',0)
      latm(7,i)=nqu
      do 338 j=1,6
      qatm(j,nqu)=qatom(j)
  338 continue
      latm(7,i)=nqu
      latm(8,i)=nqu
      latm(9,i)=nqu
  335 continue
      do 339 i=1,natm
      do 340 j=1,3
      qatom(2*j-1)=qatm(2*j-1,latm(7,i))
      qatom(2*j)=qatm(2*j,latm(7,i))
      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'S') then
      qatom(2*j-1)=0.0d0
      qatom(2*j)=0.0d0
      endif
  340 continue
      do 341 j=1,nqu
      if(qatom(1).eq.qatm(1,j).and.qatom(2).eq.qatm(2,j).and.
     + qatom(3).eq.qatm(3,j).and.qatom(4).eq.qatm(4,j).and.
     + qatom(5).eq.qatm(5,j).and.qatom(6).eq.qatm(6,j)) then
      latm(7,i)=j
      latm(8,i)=j
      goto 339
      endif
  341 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call errquit('Increase mqu',0)
      latm(7,i)=nqu
      latm(8,i)=nqu
      do 342 j=1,6
      qatm(j,nqu)=qatom(j)
  342 continue
  339 continue
c
      length=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:length),form='formatted',
     + status='new',err=9999)
c
      call swatch(topdat,toptim)
c
      write(lfntop,1000) filtop(1:length),4.0,topdat,toptim,ffield(1:10)
 1000 format('TOPOLOGY FILE ',a,/,
     + 'GENERATED AUTOMICALLY BY NWCHEM PREPARE MODULE',//,
     + f12.6,3a10)
c     
      write(lfntop,1004) q14fac,releps,0.0d0,0
 1004 format(3f12.6,i5)
      write(lfntop,1104) nats
      write(lfntop,1104) nqu
 1104 format(i5)
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
      do 322 i=1,nqu
      write(lfntop,1304) (qatm(j,i),j=1,6),i
 1304 format(3(0pf12.6,1pe12.5),i5)
  322 continue
c
      nd=0
      do 29 i=ndih+1,ndiht
      nd=nd+max(ldih(7,i),ldih(8,i),ldih(9,i))
   29 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1008) natmt-natm,nbndt-nbnd,nangt-nang,
     + nd+nimpt-nimp,0,n3rdt-n3rd,nnon
      else
      write(lfntop,1008) natmt-natm,nbndt-nbnd,nangt-nang,nd,
     + nimpt-nimp,n3rdt-n3rd,nnon
 1008 format(5i7,2i10)
      endif
      nd=0
      do 30 i=1,ndih
      nd=nd+max(ldih(7,i),ldih(8,i),ldih(9,i))
   30 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1009) natm,nbnd,nang,nd+nimp,0,n3rd,nexc,0,0,0
      else
      write(lfntop,1009) natm,nbnd,nang,nd,nimp,n3rd,nexc,0,0,0
 1009 format(5i7,3i10,2i5)
      endif
c
      if(natmt.gt.natm) then
      do 4 i=natm+1,natmt
      write(lfntop,1007) cseq(nseq+1),catm(1,i)(1:5),latm(3,i),
     + latm(7,i),latm(7,i),latm(7,i),i-natm
 1007 format(a10,a5,1x,25x,4i5,i7)
c      write(lfntop,1007) latm(3,i),0,0,0,0,cseq(nseq+1),catm(1,i)(1:5),
c     + i-natm
c 1007 format(5i7,a10,a5,1x,i7)
c      write(lfntop,1217) (qatm(2*j-1,i),j=1,3),(qatm(2*j,i),j=1,3)
c 1217 format(6f12.6)
    4 continue
c
      if(nbndt.gt.nbnd) then
      do 20 i=nbnd+1,nbndt
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 120 j=1,3
      if(latt(1,lats(j,latm(3,lbnd(1,i)-natm))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)-natm))).ne.1) lbnd(4,i)=0
  120 continue
      endif
      write(lfntop,1112) lbnd(1,i)-natm,lbnd(2,i)-natm,lbnd(4,i),i-nbnd
 1112 format(4i7)
      write(lfntop,1113) (rbnd(j,i),j=1,6)
 1113 format(3(0pf12.6,1pe12.5))
   20 continue
      endif
c
      if(nangt.gt.nang) then
      do 21 i=nang+1,nangt
      write(lfntop,1114) lang(1,i)-natm,lang(2,i)-natm,lang(3,i)-natm,
     + lang(5,i),i-nang
 1114 format(5i7)
      write(lfntop,1115) (rang(j,i),j=1,6)
 1115 format(3(0pf10.6,1pe12.5))
   21 continue
      endif
c
      if(ndiht.gt.ndih) then
      do 22 i=ndih+1,ndiht
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 23 k=1,nd
      write(lfntop,1116) (ldih(j,i)-natm,j=1,4),ldih(6,i),i-ndih
 1116 format(6i7)
      write(lfntop,1117)
     +  (kdih(k,j,i),rdih(k,2*j-1,i),rdih(k,2*j,i),j=1,3)
 1117 format(3(i3,0pf10.6,1pe12.5))
   23 continue
   22 continue
      endif
c
      if(nimpt.gt.nimp) then
      do 24 i=nimp+1,nimpt
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1118) limp(2,i)-natm,limp(3,i)-natm,limp(1,i)-natm,
     + limp(4,i)-natm,limp(6,i),i-nimp
      write(lfntop,1119) (limp(6+j,i),rimp(2*j-1,i),rimp(2*j,i),j=1,3)
      else
      write(lfntop,1118) (limp(j,i)-natm,j=1,4),limp(6,i),i-nimp
 1118 format(6i7)
      write(lfntop,1119) (0,rimp(2*j-1,i),rimp(2*j,i),j=1,3)
 1119 format(3(i3,0pf10.6,1pe12.5))
      endif
   24 continue
      endif
c
      if(n3rdt.gt.n3rd) then
      write(lfntop,1120) (l3rd(1,i)-natm,i=n3rd+1,n3rdt)
      write(lfntop,1120) (l3rd(2,i)-natm,i=n3rd+1,n3rdt)
 1120 format(11i7)
      endif
c
      if(nnon.gt.0) then
      write(lfntop,1121) (lnon(1,i)-natm,i=1,nnon)
      write(lfntop,1121) (lnon(2,i)-natm,i=1,nnon)
 1121 format(11i7)
      endif
c
      write(lfntop,1122) cseq(nseq+1),wcorr,slvnam
 1122 format(a10,3f12.6,a3)
      endif
c
      ifrac=1
      isgm=0
      lprt=.false.
      qsum(1)=0.0d0
      qsum(2)=0.0d0
      qsum(3)=0.0d0
      qtot(1)=0.0d0
      qtot(2)=0.0d0
      qtot(3)=0.0d0
      if(natm.gt.0) then
      if(util_print('topology',print_high)) then
      write(lfnout,2000)
 2000 format(' Segment charges',t40,7x,'set 1',7x,'set 2',7x,'set 3',/)
      endif
      do 10 i=1,natm
      imol=latm(6,i)
      if(imol.lt.0) then
      imol=-imol
      ifrac=ifrac+1
      endif
      iq=0
      iql(1)=latm(7,i)
      iql(2)=latm(8,i)
      iql(3)=latm(9,i)
c
c     bits set in iq : 0 set 1 s-w is dummy
c                      1 set 2 s-w is dummy
c                      2 set 3 s-w is dummy
c                      3 set 1 s-s non-self is dummy
c                      4 set 2 s-s non-self is dummy
c                      5 set 3 s-s non-self is dummy
c                      6 set 1 s-s self is dummy
c                      7 set 2 s-s self is dummy
c                      8 set 3 s-s self is dummy
c
      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'S') iq=ior(iq,9)
      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'S') iq=ior(iq,18)
      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'S') iq=ior(iq,36)
      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'D') iq=ior(iq,73)
      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'D') iq=ior(iq,146)
      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'D') iq=ior(iq,292)
c
c      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'S') iq=iq+4
c      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'S') iq=iq+2
c      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'S') iq=iq+1
c
c      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'D') iq=iq+36
c      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'D') iq=iq+18
c      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'D') iq=iq+9
c      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'S') then
c      iq=3
c      iql(1)=nqu
c      endif
c      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'S') then
c      iq=3
c      iql(2)=nqu
c      endif
c      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'S') then
c      iq=3
c      iql(3)=nqu
c      endif
c      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'S') iq=3
c      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'S') iq=3
c      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'S') iq=3
      write(lfntop,1024) cseq(latm(5,i)),catm(1,i)(1:5),ifrac,imol,
     + lseq(1,latm(5,i)),latm(1,i),latm(2,i),latm(3,i),iql,iq,i
 1024 format(a10,a5,1x,10i5,i7)
c      write(lfntop,1010) latm(3,i),imol,latm(1,i),latm(2,i),
c     + lseq(1,latm(5,i)),cseq(latm(5,i)),catm(1,i)(1:5),i,iq,ifrac
c 1010 format(5i7,a10,a5,1x,i7,2i5)
c      if(iq.eq.0) then
c      write(lfntop,1011) (qatm(2*j-1,i),j=1,3),(qatm(2*j,i),j=1,3)
c 1011 format(6f12.6)
c      else
c      do 101 j=1,3
c      qq(1,j)=qatm(2*j-1,i)
c      qq(1,3+j)=qatm(2*j,i)
c  101 continue
c      do 102 j=1,6
c      qq(2,j)=qq(1,j)
c      qq(3,j)=qq(1,j)
c  102 continue
c      do 103 j=1,3
c      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'S') then
c      qq(1,j)=0.0d0
c      qq(1,j+3)=0.0d0
c      qq(2,j)=0.0d0
c      qq(2,j+3)=0.0d0
c      endif
c  103 continue
c      write(lfntop,1011) (qq(1,j),j=1,6)
c      write(lfntop,1011) (qq(2,j),j=1,6)
c      write(lfntop,1011) (qq(3,j),j=1,6)
c      endif
      qtot(1)=qtot(1)+qatm(1,iql(1))
      qtot(2)=qtot(2)+qatm(3,iql(2))
      qtot(3)=qtot(3)+qatm(5,iql(3))
      if(isgm.ne.lseq(1,latm(5,i))) then
      if(isgm.ne.0) then
      if(util_print('topology',print_high)) then
      write(lfnout,2001) lseq(1,latm(5,i-1)),cseq(latm(5,i-1)),qsum
 2001 format(i7,1x,a10,t40,3f12.6)
      endif
      endif
      isgm=lseq(1,latm(5,i))
      qsum(1)=0.0d0
      qsum(2)=0.0d0
      qsum(3)=0.0d0
      endif
      qsum(1)=qsum(1)+qatm(1,iql(1))
      qsum(2)=qsum(2)+qatm(3,iql(2))
      qsum(3)=qsum(3)+qatm(5,iql(3))
   10 continue
      if(util_print('topology',print_high)) then
      write(lfnout,2001) lseq(1,latm(5,natm)),cseq(latm(5,natm)),qsum
      write(lfnout,2002)
 2002 format(t25,'------------------------------------',/)
      endif
      if(util_print('topology',print_default)) then
      write(lfnout,2003) qtot
 2003 format(' Total charge',t40,3f12.6,/)
      endif
      endif
c
      do 11 i=1,nbnd
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 111 j=1,3
      if(latt(1,lats(j,latm(3,lbnd(1,i)))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)))).ne.1) lbnd(4,i)=0
  111 continue
      endif
      write(lfntop,1012) lbnd(1,i),lbnd(2,i),lbnd(4,i),i
 1012 format(4i7)
      write(lfntop,1013) (rbnd(j,i),j=1,6)
 1013 format(3(0pf12.6,1pe12.5))
   11 continue
c
      do 12 i=1,nang
      write(lfntop,1014) lang(1,i),lang(2,i),lang(3,i),lang(5,i),i
 1014 format(5i7)
      write(lfntop,1015) (rang(j,i),j=1,6)
 1015 format(3(0pf10.6,1pe12.5))
   12 continue
c
      do 13 i=1,ndih
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 14 k=1,nd
      write(lfntop,1016) (ldih(j,i),j=1,4),ldih(6,i),i
 1016 format(6i7)
      write(lfntop,1017)
     +  (kdih(k,j,i),rdih(k,2*j-1,i),rdih(k,2*j,i),j=1,3)
 1017 format(3(i3,0pf10.6,1pe12.5))
   14 continue
   13 continue
c
      do 15 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1018) limp(2,i),limp(3,i),limp(1,i),limp(4,i),
     + limp(6,i),i
      write(lfntop,1019) (limp(6+j,i),rimp(2*j-1,i),rimp(2*j,i),j=1,3)
      else
      write(lfntop,1018) (limp(j,i),j=1,4),limp(6,i),i
 1018 format(6i7)
      write(lfntop,1019) (0,rimp(2*j-1,i),rimp(2*j,i),j=1,3)
 1019 format(3(i3,0pf10.6,1pe12.5))
      endif
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
      pre_wrttop=.true.
      return
 9999 continue
      pre_wrttop=.false.
      return
      end
