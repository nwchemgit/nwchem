      logical function pre_wrttop(lfnout,ffield,lfntop,filtop,
     + releps,q14fac,lseq,cseq,mseq,nseq,
     + latt,lats,catt,patt,ratt,matt,natt,nats,latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp,l3rd,m3rd,n3rd,lexc,mexc,nexc,lnon,mnon,nnon,
     + natmt,nbndt,nangt,ndiht,nimpt,n3rdt,wcorr)
c
c $Id: pre_wrttop.F,v 1.15 1998-10-05 20:19:00 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfntop
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
      integer mnon,nnon
      integer lnon(2,mnon)
      integer natmt,nbndt,nangt,ndiht,nimpt,n3rdt
      real*8 wcorr(3)
c
      integer i,j,k,l,length,nd,ifrac,imol,iq,iql(3),isgm,kn
      character*10 topdat,toptim
      real*8 qq(3,6),qsum(3),qtot(3)
      logical lprt
c
      integer nbu,nhu,ndu,nou,nqu,ll
c
      nqu=0
      nbu=0
      nhu=0
      ndu=0
      nou=0
c
c     generate list with unique charges and polarizabilities
c
      do 1 i=1,natm
      do 2 j=1,nqu
      if(qatm(1,i).eq.qatm(1,j).and.qatm(2,i).eq.qatm(2,j).and.
     + qatm(3,i).eq.qatm(3,j).and.qatm(4,i).eq.qatm(4,j).and.
     + qatm(5,i).eq.qatm(5,j).and.qatm(6,i).eq.qatm(6,j)) then
      latm(7,i)=j
      goto 1
      endif
    2 continue
      nqu=nqu+1
      latm(7,i)=nqu
      do 3 j=1,6
      qatm(j,nqu)=qatm(j,i)
    3 continue
    1 continue
      nqu=nqu+1
      do 4 j=1,6
      qatm(j,nqu)=0.0d0
    4 continue
c
c     generate list with unique bonds
c
      do 5 i=1,nbnd
      do 6 j=1,nbu
      if(rbnd(1,i).eq.rbnd(1,j).and.rbnd(2,i).eq.rbnd(2,j).and.
     + rbnd(3,i).eq.rbnd(3,j).and.rbnd(4,i).eq.rbnd(4,j).and.
     + rbnd(5,i).eq.rbnd(5,j).and.rbnd(6,i).eq.rbnd(6,j)) then
      lbnd(3,i)=j
      goto 5
      endif
    6 continue
      nbu=nbu+1
      lbnd(3,i)=nbu
      do 7 j=1,6
      rbnd(j,nbu)=rbnd(j,i)
    7 continue
    5 continue
c
c     generate list with unique angles
c
      do 8 i=1,nang
      do 9 j=1,nhu
      if(rang(1,i).eq.rang(1,j).and.rang(2,i).eq.rang(2,j).and.
     + rang(3,i).eq.rang(3,j).and.rang(4,i).eq.rang(4,j).and.
     + rang(5,i).eq.rang(5,j).and.rang(6,i).eq.rang(6,j)) then
      lang(4,i)=j
      goto 8
      endif
    9 continue
      nhu=nhu+1
      lang(4,i)=nhu
      do 10 j=1,6
      rang(j,nhu)=rang(j,i)
   10 continue
    8 continue
c
c     generate list with unique dihedrals
c
      do 11 i=1,ndih
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 12 j=1,ndu
      l=0
      do 13 k=1,nd
      if(rdih(k,1,i).eq.rdih(k,1,j).and.rdih(k,2,i).eq.rdih(k,2,j).and.
     + rdih(k,3,i).eq.rdih(k,3,j).and.rdih(k,4,i).eq.rdih(k,4,j).and.
     + rdih(k,5,i).eq.rdih(k,5,j).and.rdih(k,6,i).eq.rdih(k,6,j).and.
     + kdih(k,1,i).eq.kdih(k,1,j).and.kdih(k,2,i).eq.kdih(k,2,j).and.
     + kdih(k,3,i).eq.kdih(k,3,j)) l=l+1
   13 continue
      if(l.eq.nd) then
      ldih(5,i)=j
      goto 11
      endif
   12 continue
      ndu=ndu+1
      ldih(5,i)=ndu
      do 14 k=1,6
      do 15 j=1,6
      rdih(k,j,ndu)=rdih(k,j,i)
   15 continue
      kdih(k,1,ndu)=kdih(k,1,i)
      kdih(k,2,ndu)=kdih(k,2,i)
      kdih(k,3,ndu)=kdih(k,3,i)
   14 continue
   11 continue
c
c     generate list with unique angles
c
      do 16 i=1,nimp
      do 17 j=1,nou
      if(rimp(1,i).eq.rimp(1,j).and.rimp(2,i).eq.rimp(2,j).and.
     + rimp(3,i).eq.rimp(3,j).and.rimp(4,i).eq.rimp(4,j).and.
     + rimp(5,i).eq.rimp(5,j).and.rimp(6,i).eq.rimp(6,j).and.
     + limp(7,i).eq.limp(7,j).and.limp(8,i).eq.limp(8,j).and.
     + limp(9,i).eq.limp(9,j)) then
      limp(5,i)=j
      goto 16
      endif
   17 continue
      nou=nou+1
      limp(5,i)=nou
      do 18 j=1,6
      rimp(j,nou)=rimp(j,i)
   18 continue
      limp(7,nou)=limp(7,i)
      limp(8,nou)=limp(8,i)
      limp(9,nou)=limp(9,i)
   16 continue
c
c     open the topology file
c
      length=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:length),form='formatted',
     + status='new',err=9999)
c
      call arg_swatch(topdat,toptim)
c
c     write topology file header
c
      write(lfntop,1000) filtop(1:length),3.3,topdat,toptim,ffield(1:10)
 1000 format('TOPOLOGY FILE ',a,/,
     + 'GENERATED AUTOMICALLY BY NWCHEM PREPARE MODULE',//,
     + f12.6,3a10)
c     
      write(lfntop,1001) nats,q14fac,releps,0.0d0,0,nqu,nbu,nhu,ndu,nou
 1001 format(i5,3f12.6,6i5)
c
c     write atom type parameters
c
      do 19 i=1,nats
      write(lfntop,1002) (catt(1,lats(j,i)),ratt(lats(j,i)),j=1,3),
     + (latt(1,lats(j,i)),j=1,3)
 1002 format(3(a6,f12.6),3i5)
   19 continue
c
c     write atom type van der Waals parameters
c
      do 20 i=1,nats
      do 21 j=i,nats
      write(lfntop,1003) ((patt(1,k,lats(l,i),lats(l,j)),l=1,3),k=1,2)
      write(lfntop,1003) ((patt(2,k,lats(l,i),lats(l,j)),l=1,3),k=1,2)
 1003 format(6(1pe12.5))
   21 continue
   20 continue
c
c     write charge type parameters
c
      do 22 i=1,nqu
      write(lfntop,1004) (qatm(j,i),j=1,6),i
 1004 format(3(0pf12.6,1pe12.5),i5)
   22 continue
c
c     write bond type parameters
c
      do 23 i=1,nbu
      write(lfntop,1005) (rbnd(j,i),j=1,6),i
 1005 format(3(0pf12.6,1pe12.5),i5)
   23 continue
c
c     write angle type parameters
c
      do 24 i=1,nhu
      write(lfntop,1006) (rang(j,i),j=1,6),i
 1006 format(3(0pf12.6,1pe12.5),i5)
   24 continue
c
c     write dihedral type parameters
c
      do 25 i=1,ndu
      kn=1
  125 continue
      if(kdih(kn,1,i).lt.0.or.kdih(kn,2,i).lt.0.or.
     + kdih(kn,3,i).lt.0) then
      kn=kn+1
      if(kn.lt.6) goto 125
      endif
      write(lfntop,1007)
     +  ((kdih(k,j,i),j=1,3),(rdih(k,j,i),j=1,6),i,k=1,kn)
 1007 format((3i3,3(0pf12.6,1pe12.5),i5))
   25 continue
c
c     write improper dihedral angle type parameters
c
      do 26 i=1,nou
      write(lfntop,1008) (limp(j,i),j=7,9),(rimp(j,i),j=1,6),i
 1008 format(3i3,3(0pf12.6,1pe12.5),i5)
   26 continue
c
c     write solvent: number of atoms, bonds, angles, dihedrals, impropers
c
      write(lfntop,1009) natmt-natm,nbndt-nbnd,nangt-nang,ndiht-ndih,
     + nimpt-nimp,n3rdt-n3rd,nnon
 1009 format(5i7,2i10)
c
c     write solute: number of atoms, bonds, angles, dihedrals, impropers
c
      nd=0
      do 27 i=1,ndih
      nd=nd+max(ldih(7,i),ldih(8,i),ldih(9,i))
   27 continue
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1010) natm,nbnd,nang,nd+nimp,0,n3rd,nexc
      else
      write(lfntop,1010) natm,nbnd,nang,nd,nimp,n3rd,nexc
 1010 format(5i7,2i10)
      endif
c
      if(natmt.gt.natm) then
      do 28 i=natm+1,natmt
      write(lfntop,1011) latm(3,i),0,0,0,0,cseq(nseq+1),catm(1,i)(1:5),
     + latm(7,i),i-natm
 1011 format(5i7,a10,a5,1x,2i7)
c      write(lfntop,1012) (qatm(2*j-1,i),j=1,3),(qatm(2*j,i),j=1,3)
c 1012 format(6f12.6)
   28 continue
c
      if(nbndt.gt.nbnd) then
      do 29 i=nbnd+1,nbndt
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 30 j=1,3
      if(latt(1,lats(j,latm(3,lbnd(1,i)-natm))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)-natm))).ne.1) lbnd(4,i)=0
   30 continue
      endif
      write(lfntop,1013) lbnd(1,i)-natm,lbnd(2,i)-natm,lbnd(4,i),
     + lbnd(3,i),i-nbnd
 1013 format(5i7)
c      write(lfntop,1014) (rbnd(j,i),j=1,6)
c 1014 format(3(0pf12.6,1pe12.5))
   29 continue
      endif
c
      if(nangt.gt.nang) then
      do 31 i=nang+1,nangt
      write(lfntop,1015) lang(1,i)-natm,lang(2,i)-natm,lang(3,i)-natm,
     + lang(5,i),lang(4,i),i-nang
 1015 format(6i7)
      write(lfntop,1016) (rang(j,i),j=1,6)
 1016 format(3(0pf10.6,1pe12.5))
   31 continue
      endif
c
      if(ndiht.gt.ndih) then
      do 32 i=ndih+1,ndiht
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 33 k=1,nd
      write(lfntop,1017) (ldih(j,i)-natm,j=1,4),ldih(6,i),ldih(5,i),
     + i-ndih
 1017 format(7i7)
c      write(lfntop,1018)
c     +  (kdih(k,j,i),rdih(k,2*j-1,i),rdih(k,2*j,i),j=1,3)
c 1018 format(3(i3,0pf10.6,1pe12.5))
   33 continue
   32 continue
      endif
c
      if(nimpt.gt.nimp) then
      do 34 i=nimp+1,nimpt
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1019) limp(2,i)-natm,limp(3,i)-natm,limp(1,i)-natm,
     + limp(4,i)-natm,limp(6,i),limp(5,i),i-nimp
c      write(lfntop,1020) (limp(6+j,i),rimp(2*j-1,i),rimp(2*j,i),j=1,3)
      else
      write(lfntop,1019) (limp(j,i)-natm,j=1,4),limp(6,i),limp(5,i),
     + i-nimp
 1019 format(6i7)
c      write(lfntop,1020) (0,rimp(2*j-1,i),rimp(2*j,i),j=1,3)
c 1020 format(3(i3,0pf10.6,1pe12.5))
      endif
   34 continue
      endif
c
      if(n3rdt.gt.n3rd) then
      write(lfntop,1021) (l3rd(1,i)-natm,i=n3rd+1,n3rdt)
      write(lfntop,1021) (l3rd(2,i)-natm,i=n3rd+1,n3rdt)
 1021 format(11i7)
      endif
c
      if(nnon.gt.0) then
      write(lfntop,1022) (lnon(1,i)-natm,i=1,nnon)
      write(lfntop,1022) (lnon(2,i)-natm,i=1,nnon)
 1022 format(11i7)
      endif
c
      write(lfntop,1023) cseq(nseq+1),wcorr
 1023 format(a10,3f12.6)
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
      if(util_print('files',print_default)) then
      write(lfnout,2000)
 2000 format(' Segment charges',t25,7x,'set 1',7x,'set 2',7x,'set 3',/)
      endif
c
      do 35 i=1,natm
c
      imol=latm(6,i)
      if(imol.lt.0) then
      imol=-imol
      ifrac=ifrac+1
      endif
      iq=0
      iql(1)=latm(7,i)
      iql(2)=latm(7,i)
      iql(3)=latm(7,i)
      if(catt(1,lats(1,latm(3,i)))(6:6).eq.'S') then
      iq=3
      iql(1)=nqu
      endif
      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'S') then
      iq=3
      iql(2)=nqu
      endif
      if(catt(1,lats(3,latm(3,i)))(6:6).eq.'S') then
      iq=3
      iql(3)=nqu
      endif
      write(lfntop,1024) latm(3,i),imol,latm(1,i),latm(2,i),
     + lseq(1,latm(5,i)),cseq(latm(5,i)),catm(1,i)(1:5),i,iq,ifrac,
     + iql
 1024 format(5i7,a10,a5,1x,i7,5i5)
c
c      if(iq.eq.0) then
c      write(lfntop,1025) (qatm(2*j-1,i),j=1,3),(qatm(2*j,i),j=1,3)
c 1025 format(6f12.6)
c      else
c      do 136 j=1,3
c      qq(1,j)=qatm(2*j-1,i)
c      qq(1,3+j)=qatm(2*j,i)
c   36 continue
c      do 37 j=1,6
c      qq(2,j)=qq(1,j)
c      qq(3,j)=qq(1,j)
c   37 continue
c      do 38 j=1,3
c      if(catt(1,lats(j,latm(3,i)))(6:6).eq.'S') then
c      qq(1,j)=0.0d0
c      qq(1,j+3)=0.0d0
c      qq(2,j)=0.0d0
c      qq(2,j+3)=0.0d0
c      endif
c   38 continue
c      write(lfntop,1026) (qq(1,j),j=1,6)
c      write(lfntop,1026) (qq(2,j),j=1,6)
c      write(lfntop,1026) (qq(3,j),j=1,6)
c 1026 format(6f12.6)
c      endif
c
      qtot(1)=qtot(1)+qatm(1,iql(1))
      qtot(2)=qtot(2)+qatm(3,iql(2))
      qtot(3)=qtot(3)+qatm(5,iql(3))
      if(isgm.ne.lseq(1,latm(5,i))) then
      if(isgm.ne.0) then
      if(util_print('files',print_default)) then
      write(lfnout,2001) lseq(1,latm(5,i-1)),cseq(latm(5,i-1)),qsum
 2001 format(i7,1x,a10,t25,3f12.6)
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
   35 continue
c
      if(util_print('files',print_default)) then
      write(lfnout,2001) lseq(1,latm(5,natm)),cseq(latm(5,natm)),qsum
      write(lfnout,2002) qtot
 2002 format(t25,'------------------------------------',/,
     + ' Total charge',t25,3f12.6,/)
      endif
c
      do 36 i=1,nbnd
      if(lbnd(4,i).eq.0) then
      lbnd(4,i)=1
      do 37 j=1,3
      if(latt(1,lats(j,latm(3,lbnd(1,i)))).ne.1.and.
     + latt(1,lats(j,latm(3,lbnd(2,i)))).ne.1) lbnd(4,i)=0
   37 continue
      endif
      write(lfntop,1027) lbnd(1,i),lbnd(2,i),lbnd(4,i),lbnd(3,i),i
 1027 format(5i7)
c      write(lfntop,1028) (rbnd(j,i),j=1,6)
c 1028 format(3(0pf12.6,1pe12.5))
   36 continue
c
      do 38 i=1,nang
      write(lfntop,1029) lang(1,i),lang(2,i),lang(3,i),lang(5,i),
     + lang(4,i),i
 1029 format(6i7)
c      write(lfntop,1030) (rang(j,i),j=1,6)
c 1030 format(3(0pf10.6,1pe12.5))
   38 continue
c
c
      do 39 i=1,ndih
      nd=max(ldih(7,i),ldih(8,i),ldih(9,i))
c      do 40 k=1,nd
      write(lfntop,1031) (ldih(j,i),j=1,4),ldih(6,i),ldih(5,i),i
 1031 format(7i7)
c      write(lfntop,1032)
c     +  (kdih(k,j,i),rdih(k,2*j-1,i),rdih(k,2*j,i),j=1,3)
c 1032 format(3(i3,0pf10.6,1pe12.5))
c   40 continue
   39 continue
c
      do 41 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      write(lfntop,1033) limp(2,i),limp(3,i),limp(1,i),limp(4,i),
     + limp(6,i),limp(5,i),i
c      write(lfntop,1034) (limp(6+j,i),rimp(2*j-1,i),rimp(2*j,i),j=1,3)
      else
      write(lfntop,1033) (limp(j,i),j=1,4),limp(6,i),limp(5,i),i
 1033 format(7i7)
c      write(lfntop,1034) (0,rimp(2*j-1,i),rimp(2*j,i),j=1,3)
c 1034 format(3(i3,0pf10.6,1pe12.5))
      endif
   41 continue
c
      if(n3rd.gt.0) then
      write(lfntop,1035) (l3rd(1,i),i=1,n3rd)
      write(lfntop,1035) (l3rd(2,i),i=1,n3rd)
 1035 format(11i7)
      endif
c
      if(nexc.gt.0) then
      write(lfntop,1036) (lexc(1,i),i=1,nexc)
      write(lfntop,1036) (lexc(2,i),i=1,nexc)
 1036 format(11i7)
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
