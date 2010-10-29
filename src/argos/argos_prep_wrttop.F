      logical function argos_prep_wrttop(lfnout,title,lfntop,filtop,
     + lfnhdb,lfnhop,filhop,lfnbsg,filbsg,
     + releps,q14fac,nparms,mparms,
     + lseq,cseq,mseq,nseq,
     + mqu,latt,lats,catt,patt,ratt,matt,natt,mats,nats,latm,catm,qatm,
     + matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,
     + l3rd,m3rd,n3rd,lexc,mexc,nexc,lnon,mnon,nnon,
     + natmt,nbndt,nangt,ndiht,nimpt,n3rdt,wcorr,slvnam,itopol)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "argos_prep_common.fh"
c
      integer lfnout,lfntop,nparms,mparms,lfnhdb,lfnhop,lfnbsg,itopol
      character*80 title(2,3)
      character*255 filtop,filhop,filbsg
      real*8 releps,q14fac
      integer mseq,nseq,mqu
      integer lseq(6,mseq)
      character*10 cseq(mseq)
      integer matt,natt,mats,nats
      integer latt(3,matt),lats(nparms,mats)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(11,matm)
      character*6 catm(mparms,matm)
      real*8 qatm(nparms,2,mqu)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(nparms,2,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(nparms,4,mang)
      integer mdih,ndih
      integer ldih(6,mdih),nfdih(nparms,mdih),kdih(6,nparms,mdih)
      real*8 rdih(6,nparms,2,mdih)
      integer mimp,nimp
      integer limp(6,mimp),kimp(nparms,mimp)
      real*8 rimp(nparms,2,mimp)
      integer m3rd,n3rd
      integer l3rd(2,m3rd)
      integer mexc,nexc
      integer lexc(2,mexc)
      integer mnon,nnon
      integer lnon(2,mnon)
      integer natmt,nbndt,nangt,ndiht,nimpt,n3rdt
      real*8 wcorr(10)
      character*3 slvnam
c
      integer i,j,k,kk,l,length,nd,ifrac,imol,iq,isgm,masks,maskd
      character*10 topdat,toptim
      character*18 user
#ifdef USE_POSIXF
      integer ilen,ierror
#endif
      real*8 qsum(10),qtot(10),qatom(10,2)
      logical lprt
c
      integer nqu,iql(3),nbsg(30,2)
c
      integer ip,ipfr,ipto,ipnum
c
      ipfr=1
      ipto=nparms
      if(itopol.eq.0) ipto=1
      if(itopol.eq.1) ipto=1
      ipnum=ipto-ipfr+1
c
      nqu=0
c
      do 2033 i=1,30
      nbsg(i,1)=0
      nbsg(i,2)=0
 2033 continue
      nbsg(28,1)=1
      nbsg(28,2)=1
c
c     generate list with unique charges and polarizabilities
c
      do 331 i=1,natmt
      do 332 j=1,nqu
      kk=0
      do 432 k=ipfr,ipto
      if(qatm(k,1,i).eq.qatm(k,1,j).and.
     + qatm(k,2,i).eq.qatm(k,2,j)) kk=kk+1
  432 continue
      if(kk.eq.ipnum) then
      latm(7,i)=j
      latm(8,i)=j
      latm(9,i)=j
      goto 331
      endif
  332 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call md_abort('Increase mqu',0)
      latm(7,i)=nqu
      do 333 j=ipfr,ipto
      qatm(j,1,nqu)=qatm(j,1,i)
      qatm(j,2,nqu)=qatm(j,2,i)
  333 continue
  331 continue
      nqu=nqu+1
      do 334 j=ipfr,ipto
      qatm(j,1,nqu)=0.0d0
      qatm(j,2,nqu)=0.0d0
  334 continue
      do 335 i=1,natm
      do 336 j=ipfr,ipto
      qatom(j,1)=qatm(j,1,latm(7,i))
      qatom(j,2)=qatm(j,2,latm(7,i))
      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'D') then
      qatom(j,1)=0.0d0
      qatom(j,2)=0.0d0
      endif
  336 continue
      do 337 j=1,nqu
      kk=0
      do 437 k=ipfr,ipto
      if(qatom(k,1).eq.qatm(k,1,j).and.qatom(k,2).eq.qatm(k,2,j))
     +  kk=kk+1
  437 continue
      if(kk.eq.ipnum) then
      latm(7,i)=j
      latm(8,i)=j
      latm(9,i)=j
      goto 335
      endif
  337 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call md_abort('Increase mqu',0)
      latm(7,i)=nqu
      do 338 j=ipfr,ipto
      qatm(j,1,nqu)=qatom(j,1)
      qatm(j,2,nqu)=qatom(j,2)
  338 continue
      latm(7,i)=nqu
      latm(8,i)=nqu
      latm(9,i)=nqu
  335 continue
      do 339 i=1,natm
      do 340 j=ipfr,ipto
      qatom(j,1)=qatm(j,1,latm(7,i))
      qatom(j,2)=qatm(j,2,latm(7,i))
      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'S') then
      qatom(j,1)=0.0d0
      qatom(j,2)=0.0d0
      endif
  340 continue
      do 341 j=1,nqu
      kk=0
      do 441 k=ipfr,ipto
      if(qatom(k,1).eq.qatm(k,1,j).and.qatom(k,2).eq.qatm(k,2,j))
     + kk=kk+1
  441 continue
      if(kk.eq.ipnum) then
      latm(7,i)=j
      latm(8,i)=j
      goto 339
      endif
  341 continue
      nqu=nqu+1
      if(nqu.gt.mqu) call md_abort('Increase mqu',0)
      latm(7,i)=nqu
      latm(8,i)=nqu
      do 342 j=ipfr,ipto
      qatm(j,1,nqu)=qatom(j,1)
      qatm(j,2,nqu)=qatom(j,2)
  342 continue
  339 continue
c
      length=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:length),form='formatted',
     + status='new',err=9999)
c
      call swatch(topdat,toptim)
#ifdef USE_POSIXF
      call pxfgetlogin(user, ilen, ierror)
#else
      call getlog(user)
#endif
      if(user(18:18).ne.' ') user='                  '
c
      write(lfntop,1000) (title(1,i),i=1,3),
     + 4.6,topdat,toptim,ffield(1:10),user
 1000 format(a,/,a,/,a,/,f12.6,3a10,a18)
c
      write(lfntop,1104) ipnum
      write(lfntop,1104) nats
      write(lfntop,1104) nqu
      write(lfntop,1104) nseq
 1104 format(i5)
c
      nbsg(2,1)=0
      nbsg(2,2)=nats
      nbsg(19,1)=0
      nbsg(19,2)=nats
c
      write(lfntop,1004) q14fac,releps,0.0d0,0
 1004 format(3f12.6,i5)
      if(ipnum.eq.1) then
      do 1 i=1,nats
      write(lfntop,1005) i,latt(1,lats(1,i)),catt(1,lats(1,i)),
     + ratt(lats(1,i))
 1005 format(2i5,1x,a6,f12.6)
    1 continue
      do 2 i=1,nats
      do 3 j=i,nats
      write(lfntop,1006) i,j,
     +  (patt(1,k,lats(1,i),lats(1,j)),k=1,2),
     +  (patt(2,k,lats(1,i),lats(1,j)),k=1,2)
 1006 format(2i5,4(1pe12.5))
    3 continue
    2 continue
      do 322 i=1,nqu
      write(lfntop,1304) i,qatm(1,1,i),qatm(1,2,i)
 1304 format(i5,0pf12.6,1pe12.5)
  322 continue
      else
      do 51 i=1,nats
      write(lfntop,1055) i,(latt(1,lats(j,i)),catt(1,lats(j,i)),
     + ratt(lats(j,i)),j=ipfr,ipto)
 1055 format(2i5,1x,a6,f12.6,/,(5x,i5,1x,a6,f12.6))
   51 continue
      do 52 i=1,nats
      do 53 j=i,nats
      write(lfntop,1056) i,j,
     + ((patt(1,k,lats(l,i),lats(l,j)),k=1,2),
     +  (patt(2,k,lats(l,i),lats(l,j)),k=1,2),l=ipfr,ipto)
 1056 format(2i5,4(1pe12.5),/,(10x,4(1pe12.5)))
   53 continue
   52 continue
      do 5322 i=1,nqu
      write(lfntop,1354) i,(qatm(j,1,i),qatm(j,2,i),j=ipfr,ipto)
 1354 format(i5,0pf12.6,1pe12.5,/,(5x,0pf12.6,1pe12.5))
 5322 continue
      endif
      nbsg(29,1)=natmt-natm
      nbsg(29,2)=lseq(3,1)
      do 5342 i=1,nseq
      write(lfntop,1355) i,cseq(i),(lseq(j,i),j=1,6)
 1355 format(i7,a10,i7,i5,2i7,2i5)
      if(i.gt.1) nbsg(29,2)=max(nbsg(29,2),lseq(3,i)-lseq(3,i-1))
 5342 continue
c
      nd=0
      do 29 i=ndih+1,ndiht
      j=0
      do 529 k=ipfr,ipto
      j=max(j,nfdih(k,i))
  529 continue
      nd=nd+j
   29 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1008) natmt-natm,nbndt-nbnd,nangt-nang,
     + nd+nimpt-nimp,0,n3rdt-n3rd,nnon
      else
      write(lfntop,1008) natmt-natm,nbndt-nbnd,nangt-nang,nd,
     + nimpt-nimp,n3rdt-n3rd,nnon
 1008 format(5i7,2i10)
      endif
      nbsg(1,1)=natmt-natm
      nbsg(1,2)=natm
      nbsg(11,1)=(natmt-natm)*(natmt-natm-1)/2-(n3rdt-n3rd)-nnon
      nbsg(12,1)=1
      nbsg(12,2)=nseq
      nd=0
      do 30 i=1,ndih
      j=0
      do 530 k=ipfr,ipto
      j=max(j,nfdih(k,i))
  530 continue
      nd=nd+j
   30 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1009) natm,nbnd,nang,nd+nimp,0,n3rd,nexc,0,0,0
      else
      write(lfntop,1009) natm,nbnd,nang,nd,nimp,n3rd,nexc,0,0,0
 1009 format(5i7,3i10,2i5)
      endif
      nbsg(11,2)=nexc
c
      if(natmt.gt.natm) then
      do 4 i=natm+1,natmt
      write(lfntop,1007) cseq(nseq+1),catm(1,i)(1:5),latm(3,i),
     + latm(7,i),latm(7,i),latm(7,i),i-natm
 1007 format(a10,a5,1x,25x,4i5,i7)
    4 continue
c
      if(nbndt.gt.nbnd) then
      do 20 i=nbnd+1,nbndt
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 120 j=ipfr,ipto
      if(latt(1,lats(j,latm(3,lbnd(1,i)-natm))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)-natm))).ne.1) lbnd(4,i)=0
  120 continue
      endif
      write(lfntop,1112) lbnd(1,i)-natm,lbnd(2,i)-natm,lbnd(4,i),i-nbnd
 1112 format(4i7)
      write(lfntop,1113) (rbnd(j,1,i),rbnd(j,2,i),j=ipfr,ipto)
 1113 format(0pf12.6,1pe12.5)
      if(latt(1,lats(1,latm(3,lbnd(1,i)-natm))).eq.1.or.
     + latt(1,lats(1,latm(3,lbnd(2,i)-natm))).eq.1) then
      nbsg(3,1)=nbsg(3,1)+1
      else
      nbsg(4,1)=nbsg(4,1)+1
      endif
      if(lbnd(4,i).eq.1) nbsg(13,1)=nbsg(13,1)+1
   20 continue
      nbsg(16,1)=nbndt-nbnd
      endif
c
      if(nangt.gt.nang) then
      do 21 i=nang+1,nangt
      write(lfntop,1114) lang(1,i)-natm,lang(2,i)-natm,lang(3,i)-natm,
     + lang(5,i),i-nang
 1114 format(5i7)
      if(ffield(1:6).ne.'charmm') then
      write(lfntop,1115) (rang(j,1,i),rang(j,2,i),j=ipfr,ipto)
 1115 format(0pf10.6,1pe12.5)
      else
      write(lfntop,1125) (rang(j,1,i),rang(j,2,i),rang(j,3,i),
     + rang(j,4,i),j=ipfr,ipto)
 1125 format(0pf10.6,1pe12.5,0pf10.6,1pe12.5)
      endif
      if(latt(1,lats(1,latm(3,lang(1,i)-natm))).eq.1.or.
     + latt(1,lats(1,latm(3,lang(3,i)-natm))).eq.1) then
      nbsg(5,1)=nbsg(5,1)+1
      else
      nbsg(6,1)=nbsg(6,1)+1
      endif
   21 continue
      nbsg(17,1)=nangt-nang
      endif
c
      if(ndiht.gt.ndih) then
      do 22 i=ndih+1,ndiht
      nd=0
      do 522 k=ipfr,ipto
      nd=max(nd,nfdih(k,i))
  522 continue
      do 23 k=1,nd
      write(lfntop,1116) (ldih(j,i)-natm,j=1,4),ldih(6,i),i-ndih
 1116 format(6i7)
      write(lfntop,1117)
     +  (kdih(k,j,i),rdih(k,j,1,i),rdih(k,j,2,i),j=ipfr,ipto)
 1117 format(i3,0pf10.6,1pe12.5)
      if(latt(1,lats(1,latm(3,ldih(1,i)-natm))).eq.1.or.
     + latt(1,lats(1,latm(3,ldih(4,i)-natm))).eq.1) then
      nbsg(7,1)=nbsg(7,1)+1
      else
      nbsg(8,1)=nbsg(8,1)+1
      endif
   23 continue
   22 continue
      nbsg(18,1)=ndiht-ndih
      endif
c
      if(nimpt.gt.nimp) then
      do 24 i=nimp+1,nimpt
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1118) limp(2,i)-natm,limp(3,i)-natm,limp(1,i)-natm,
     + limp(4,i)-natm,limp(6,i),i-nimp
      write(lfntop,1119) (kimp(j,i),rimp(j,1,i),rimp(j,2,i),j=ipfr,ipto)
      else
      write(lfntop,1118) (limp(j,i)-natm,j=1,4),limp(6,i),i-nimp
 1118 format(6i7)
      write(lfntop,1119) (0,rimp(j,1,i),rimp(j,2,i),j=ipfr,ipto)
 1119 format(i3,0pf10.6,1pe12.5)
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
      write(lfntop,1122) cseq(nseq+1),slvnam
 1122 format(a10,a3)
      write(lfntop,2222) (wcorr(i),i=ipfr,ipto)
 2222 format(f12.6)
      endif
c
      ifrac=1
      isgm=0
      lprt=.false.
      do 388 k=ipfr,ipto
      qsum(k)=0.0d0
      qtot(k)=0.0d0
  388 continue
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
c     bits set in iq : 0 set 1 s-w is dummy               (1)
c                      1 set 2 s-w is dummy               (2)
c                      2 set 3 s-w is dummy               (4)
c                      3 set 1 s-s non-self is dummy      (8)
c                      4 set 2 s-s non-self is dummy     (16)
c                      5 set 3 s-s non-self is dummy     (32)
c                      6 set 1 s-s self is dummy         (64)
c                      7 set 2 s-s self is dummy        (128)
c                      8 set 3 s-s self is dummy        (256)
c
      masks=9
      maskd=73
      do 1110 l=ipfr,ipto
      if(catt(1,lats(l,latm(3,i)))(6:6).eq.'S') iq=ior(iq,masks)
      if(catt(1,lats(l,latm(3,i)))(6:6).eq.'D') iq=ior(iq,maskd)
      masks=2*masks
      maskd=2*maskd
 1110 continue
c
      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'Q'.or.
     + catt(1,lats(1,latm(3,i)))(6:6).eq.'H') then
      write(lfntop,1024) cseq(latm(5,i)),catm(1,i)(1:5),ifrac,imol,
     + lseq(1,latm(5,i)),latm(1,i),latm(2,i),latm(3,i),nqu,nqu,nqu,0,
     + latm(10,i),latm(11,i),i
      else
      write(lfntop,1024) cseq(latm(5,i)),catm(1,i)(1:5),ifrac,imol,
     + lseq(1,latm(5,i)),latm(1,i),latm(2,i),latm(3,i),iql,iq,
     + latm(10,i),latm(11,i),i
      endif
 1024 format(a10,a5,1x,i3,4i7,6i5,i3,i7)
c 1024 format(a10,a5,1x,12i5,i7)
c
      do 589 k=ipfr,ipto
      qtot(k)=qtot(k)+qatm(k,1,iql(1))
  589 continue
      if(isgm.ne.lseq(1,latm(5,i))) then
      if(isgm.ne.0) then
      if(util_print('topology',print_high)) then
      write(lfnout,2001) lseq(1,latm(5,i-1)),cseq(latm(5,i-1)),
     + (qsum(k),k=ipfr,ipto)
 2001 format(i7,1x,a10,t40,10f12.6)
      endif
      endif
      isgm=lseq(1,latm(5,i))
      do 389 k=ipfr,ipto
      qsum(k)=0.0d0
  389 continue
      endif
      do 399 k=ipfr,ipto
      qsum(k)=qsum(k)+qatm(k,1,iql(1))
  399 continue
   10 continue
      if(util_print('topology',print_high)) then
      write(lfnout,2001) lseq(1,latm(5,natm)),cseq(latm(5,natm)),
     + (qsum(k),k=ipfr,ipto)
      write(lfnout,2002)
 2002 format(t40,'------------------------------------',/)
      endif
      if(util_print('topology',print_default)) then
      write(lfnout,2003) (qtot(k),k=ipfr,ipto)
 2003 format(' Total charge',t40,10f12.6,/)
      endif
      endif
c
      do 11 i=1,nbnd
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 111 j=ipfr,ipto
      if(latt(1,lats(j,latm(3,lbnd(1,i)))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)))).ne.1) lbnd(4,i)=0
  111 continue
      endif
      write(lfntop,1012) lbnd(1,i),lbnd(2,i),lbnd(4,i),i
 1012 format(4i7)
      write(lfntop,1013) (rbnd(j,1,i),rbnd(j,2,i),j=ipfr,ipto)
 1013 format(0pf12.6,1pe12.5)
      if(latt(1,lats(1,latm(3,lbnd(1,i)))).eq.1.or.
     + latt(1,lats(1,latm(3,lbnd(2,i)))).eq.1) then
      nbsg(3,2)=nbsg(3,2)+1
      else
      nbsg(4,2)=nbsg(4,2)+1
      endif
      if(lbnd(4,i).eq.1) nbsg(13,2)=nbsg(13,2)+1
   11 continue
      nbsg(16,2)=nbnd
c
      do 12 i=1,nang
      write(lfntop,1014) lang(1,i),lang(2,i),lang(3,i),lang(5,i),i
 1014 format(5i7)
      if(ffield(1:6).ne.'charmm') then
      write(lfntop,1015) (rang(j,1,i),rang(j,2,i),j=ipfr,ipto)
 1015 format(0pf10.6,1pe12.5)
      else
      write(lfntop,1035) (rang(j,1,i),rang(j,2,i),rang(j,3,i),
     + rang(j,4,i),j=ipfr,ipto)
 1035 format(0pf10.6,1pe12.5,0pf10.6,1pe12.5)
      endif
      if(latt(1,lats(1,latm(3,lang(1,i)))).eq.1.or.
     + latt(1,lats(1,latm(3,lang(3,i)))).eq.1) then
      nbsg(5,2)=nbsg(5,2)+1
      else
      nbsg(6,2)=nbsg(6,2)+1
      endif
   12 continue
      nbsg(17,2)=nang
c
      do 13 i=1,ndih
      nd=0
      do 513 k=ipfr,ipto
      nd=max(nd,nfdih(k,i))
  513 continue
      do 14 k=1,nd
      write(lfntop,1016) (ldih(j,i),j=1,4),ldih(6,i),i
 1016 format(6i7)
      write(lfntop,1017)
     +  (kdih(k,j,i),rdih(k,j,1,i),rdih(k,j,2,i),j=ipfr,ipto)
 1017 format(i3,0pf10.6,1pe12.5)
   14 continue
      if(latt(1,lats(1,latm(3,ldih(1,i)))).eq.1.or.
     + latt(1,lats(1,latm(3,ldih(4,i)))).eq.1) then
      nbsg(7,2)=nbsg(7,2)+1
      else
      nbsg(8,2)=nbsg(8,2)+1
      endif
   13 continue
      nbsg(18,2)=ndih
c
      do 15 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1018) limp(2,i),limp(3,i),limp(1,i),limp(4,i),
     + limp(6,i),i
      write(lfntop,1019) (kimp(j,i),rimp(j,1,i),rimp(j,2,i),j=ipfr,ipto)
      else
      write(lfntop,1018) (limp(j,i),j=1,4),limp(6,i),i
 1018 format(6i7)
      write(lfntop,1019) (0,rimp(j,1,i),rimp(j,2,i),j=ipfr,ipto)
 1019 format(i3,0pf10.6,1pe12.5)
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
      if(itopol.eq.3) then
      call qhop_setup(lfnhdb,lfntop,
     + mparms,ipnum,matm,natm,mseq,nseq,
     + catm,latm,cseq,lseq,mbnd,nbnd,lbnd,rbnd,mang,nang,lang,rang)
      lqhopmem=.true.
      endif
c
      close(unit=lfntop)
c
      length=index(filbsg,' ')-1
      open(unit=lfnbsg,file=filbsg(1:length),form='formatted',
     + status='unknown',err=9999)
      write(lfnbsg,2022) (nbsg(i,1),i=1,30)
      write(lfnbsg,2022) (nbsg(i,2),i=1,30)
 2022 format(12i7)
      close(unit=lfnbsg)
c
      argos_prep_wrttop=.true.
      return
 9999 continue
      argos_prep_wrttop=.false.
      return
      end
