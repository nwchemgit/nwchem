      logical function argos_prepare_mklist(lfnout,lfnsgm,lfnmat,nparms,
     + mparms,
     + lseq,cseq,mseq,nseq,
     + llnk,clnk,mlnk,nlnk,
     + latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,
     + lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,
     + natmt,nbndt,nangt,ndiht,nimpt,wcorr,
     + nval,ival,rval,ndx,itopol)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "argos_prepare_common.fh"
c
      logical argos_prepare_sgmfil,argos_prepare_rdsgm
      integer argos_prepare_atnum
      external argos_prepare_sgmfil,argos_prepare_rdsgm,
     + argos_prepare_atnum
c
      integer lfnout,lfnsgm,lfnmat,nparms,mparms,itopol
      integer mseq,nseq
      integer lseq(6,mseq)
      character*10 cseq(mseq)
      integer mlnk,nlnk
      integer llnk(3,mlnk)
      character*4 clnk(2,mlnk)
      integer matm,natm
      integer latm(11,matm)
      character*6 catm(mparms,matm)
      real*8 qatm(nparms,2,matm)
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
      integer natmt,nbndt,nangt,ndiht,nimpt
      real*8 wcorr(10)
      integer nval
      integer ival(nval),ndx(nval)
      real*8 rval(nval)
c
      integer i,j,k,l,m,link1,link2,isgm,jsgm,nb,na,nd,ni,nn,iatm,ia(10)
      integer no,ndi,ltemp,igrp,ipgrp,nia(3)
      real*8 rtemp
      character*4 atl1,atl2
      character*255 filsgm,filmat
      integer in,jn,ln,ir,itmp,nsgmpar,ipdef,lto
c
c     lang(1,i) : atom index i in angle i-j-k
c     lang(2,i) : atom index j in angle i-j-k
c     lang(3,i) : atom index k in angle i-j-k
c     lang(4,i) : constraint type
c     lang(5,i) : parameter origin
c
      if(util_print('where',print_debug)) then
      write(lfnout,1100)
 1100 format('ARGOS_PREPARE_MKLIST ENTERED ',a)
      endif
c
      igrp=0
      ipgrp=0
      do 1 i=1,nseq
c
      lseq(3,i)=natm+1
c
      filmat=cseq(i)(1:index(cseq(i),' ')-1)//'.mat '
c
c     generate segment filename
c     -------------------------
c
      if(.not.argos_prepare_sgmfil(cseq(i),lfnsgm,filsgm))
     + call md_abort('argos_prepare_sgmfil failed',9998)
c
c     append the atom and bonded lists
c     --------------------------------
c
      if(.not.argos_prepare_rdsgm(lfnout,nparms,mparms,
     + lseq(4,i),i,igrp,ipgrp,
     + lfnsgm,filsgm,lfnmat,filmat,
     + latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,limp,kimp,rimp,mimp,nimp,
     + wcorr,nsgmpar,ipdef,itopol))
     + call md_abort('argos_prepare_rdsgm failed',9999)
c
      lseq(5,i)=nsgmpar
      lseq(6,i)=ipdef
    1 continue
      lseq(3,nseq+1)=natm+1
c
      if(util_print('where',print_debug)) then
      write(lfnout,1101)
 1101 format('ARGOS_PREPARE_MKLIST LOOP 1 DONE ')
      endif
c
c     save original number of bonds in nb
c
      nb=nbnd
c
c     add chain links
c     ---------------
c
      do 2 i=1,nseq-1
      if((lseq(2,i).eq.-1.or.lseq(2,i).eq.0).and.lseq(2,i+1).ge.0) then
      link2=0
      do 3 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(3,j).eq.2) link2=j
    3 continue
      link1=0
      do 4 j=lseq(3,i+1),lseq(3,i+2)-1
      if(latm(3,j).eq.1) link1=j
    4 continue
      if(link1.gt.0.and.link2.gt.0) then
      nbnd=nbnd+1
      if(nbnd.gt.mbnd) call md_abort('increase mbnd',9999)
      lbnd(1,nbnd)=link1
      lbnd(2,nbnd)=link2
      lbnd(3,nbnd)=0
      lbnd(4,nbnd)=0
      if(util_print('topology',print_debug)) then
      write(lfnout,1004) lseq(1,i),lseq(1,i+1)
 1004 format(' Main chain link ',2i5)
      endif
      else
c
c     possible chain links appears to be missing here
c
      endif
      endif
    2 continue
c
      if(util_print('where',print_debug)) then
      write(lfnout,1102)
 1102 format('ARGOS_PREPARE_MKLIST LOOP 2 DONE ')
      endif
c
c     add explicit links
c     ------------------
c
      do 5 i=1,nlnk
c
c     find segment numbers involved in explicit link
c
      isgm=0
      jsgm=0
      do 6 j=1,nseq
      if(lseq(1,j).eq.llnk(1,i)) isgm=j
      if(lseq(1,j).eq.llnk(2,i)) jsgm=j
    6 continue
      if(isgm.eq.0.or.jsgm.eq.0) then
      call md_abort('mklist: illegal segment in explicit link',9999)
      endif
c
c     find the atoms involved in explicit link
c
      link1=0
      do 7 j=lseq(3,isgm),lseq(3,isgm+1)-1
      if(catm(1,j)(1:4).eq.clnk(1,i)) link1=j
    7 continue
      link2=0
      do 8 j=lseq(3,jsgm),lseq(3,jsgm+1)-1
      if(catm(1,j)(1:4).eq.clnk(2,i)) link2=j
    8 continue
      if(link1.eq.0.or.link2.eq.0) then
      call md_abort('mklist: illegal atom in explicit link',9999)
      endif
c
c     test if this appears to be a chain link already present
c
      if((latm(3,link1).eq.2.and.latm(3,link2).eq.1.and.
     + isgm.eq.jsgm-1).or.
     + (latm(3,link1).eq.1.and.latm(3,link2).eq.2.and.
     + jsgm.eq.isgm-1)) then
      do 9 j=nb+1,nbnd
      if((lbnd(1,j).eq.link1.and.lbnd(2,j).eq.link2).or.
     + (lbnd(1,j).eq.link2.and.lbnd(2,j).eq.link1)) then
c
c     explicit link appears to duplicate automatic chain link
c
      goto 5
      endif
    9 continue
      endif
c
c     test if bonding allowed when non-forced link
c
      if(latm(3,link1).eq.0.or.latm(3,link2).eq.0) then
      atl1=clnk(1,i)
      atl2=clnk(2,i)
      if(atl1(1:1).eq.' ') atl1(1:1)='_'
      if(atl2(1:1).eq.' ') atl2(1:1)='_'
      if(llnk(3,i).eq.0) then
      if(util_print('topology',print_low)) then
      write(lfnout,1000) llnk(1,i),atl1,llnk(2,i),atl2,
     + link1,link2,latm(3,link1),latm(3,link2)
 1000 format(' Skipped link ',i5,':',a4,i5,':',a4,' = ',i5,'-',i5,
     + ' (',i1,',',i1,')')
      endif
c      write(lfnout,1000) llnk(1,i),clnk(1,i),latm(3,link1),
c     + llnk(2,i),clnk(2,i),latm(3,link2)
c 1000 format(' Link ',i5,':',a4,'(',i1,') - ',i5,':',a4,'(',i1,')')
c      call md_abort('mklist: explicit link to nonbonding atom',9999)
      goto 5
      else
      if(util_print('topology',print_low)) then
      write(lfnout,1001) llnk(1,i),atl1,llnk(2,i),atl2,
     + link1,link2,latm(3,link1),latm(3,link2)
 1001 format(' Forced link  ',i5,':',a4,i5,':',a4,' = ',i5,'-',i5,
     + ' (',i1,',',i1,')')
      endif
c      write(lfnout,1001) llnk(1,i),clnk(1,i),latm(3,link1),
c     + llnk(2,i),clnk(2,i),latm(3,link2),link1,link2
c 1001 format(' Forced link ',
c     + i5,':',a4,'(',i1,') - ',i5,':',a4,'(',i1,') : ',i5,'-',i5)
      endif
      endif
c
c     accept explicit link
c
      nbnd=nbnd+1
      if(nbnd.gt.mbnd) call md_abort('increase mbnd',9999)
      lbnd(1,nbnd)=link1
      lbnd(2,nbnd)=link2
      lbnd(3,nbnd)=0
      lbnd(4,nbnd)=0
c
    5 continue
c
      if(util_print('where',print_debug)) then
      write(lfnout,1105)
 1105 format('ARGOS_PREPARE_MKLIST LOOP 5 DONE ')
      endif
c
      if(util_print('topology',print_low)) then
      write(lfnout,1003)
 1003 format(' ')
      endif
c
c     append the angle list
c     ---------------------
c
      na=nang
      do 13 i=nb+1,nbnd
      do 14 j=1,i-1
      if(lbnd(1,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) call md_abort('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(2,j)
      lang(4,nang)=0
      lang(5,nang)=0
      else
      if(lbnd(2,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) call md_abort('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(2,j)
      lang(4,nang)=0
      lang(5,nang)=0
      else
      if(lbnd(1,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) call md_abort('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(1,j)
      lang(4,nang)=0
      lang(5,nang)=0
      else
      if(lbnd(2,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) call md_abort('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(1,j)
      lang(4,nang)=0
      lang(5,nang)=0
      endif
      endif
      endif
      endif
   14 continue
   13 continue
c
c     append the dihedral list
c     ------------------------
c
      nd=ndih
      do 15 i=na+1,nang
      do 16 j=1,i-1
      if(lang(2,i).eq.lang(1,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) call md_abort('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(3,j)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      else
      if(lang(2,i).eq.lang(3,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) call md_abort('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(1,j)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(3,j)) then
      if(ndih.ge.mdih) call md_abort('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(1,j)) then
      if(ndih.ge.mdih) call md_abort('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(3,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      endif
      endif
      endif
      endif
   16 continue
   15 continue
c
c     append the improper dihedral list
c     ---------------------------------
c
      ni=nimp
      do 17 i=nb+1,nbnd
      do 18 j=1,2
      iatm=lbnd(j,i)
      if(latm(4,iatm).gt.0) then
      nn=0
      do 19 k=1,nbnd
      if(lbnd(1,k).eq.iatm) then
      nn=nn+1
      ia(nn)=lbnd(2,k)
      endif
      if(lbnd(2,k).eq.iatm) then
      nn=nn+1
      ia(nn)=lbnd(1,k)
      endif
   19 continue
      nia(1)=0
      nia(2)=0
      nia(3)=0
      do 79 k=1,nbnd
      if(lbnd(1,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(2,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(1,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(2,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(1,k).eq.ia(3)) nia(3)=nia(3)+1
      if(lbnd(2,k).eq.ia(3)) nia(3)=nia(3)+1
   79 continue
      do 80 k=1,2
      do 81 l=2,3
      if(nia(k).lt.nia(l)) then
      ltemp=nia(k)
      nia(k)=nia(l)
      nia(l)=ltemp
      ltemp=ia(k)
      ia(k)=ia(l)
      ia(l)=ltemp
      endif
   81 continue
   80 continue
c
c     planar
c
      if(latm(4,iatm).eq.1.and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call md_abort('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
      limp(6,nimp)=0
      endif
c
c     chiral
c
      if((latm(4,iatm).eq.2.or.latm(4,iatm).eq.3).and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call md_abort('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
      limp(6,nimp)=0
      endif
c
c     aromatic ring
c
      if(ffield(1:5).eq.'amber'.or.ffield(1:5).eq.'charm') then
      if(latm(4,iatm).eq.4.and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call md_abort('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
      limp(6,nimp)=0
      endif
      endif
c
      endif
   18 continue
   17 continue
c
c     dihedral selection
c     ------------------
c
      if(ffield(1:6).eq.'gromos') then
      ndi=ndih
      ndih=nd
      do 20 i=nd+1,ndi
      do 21 j=nd+1,ndih
      if(ldih(2,i).eq.ldih(2,j).and.ldih(3,i).eq.ldih(3,j)) then
      if(argos_prepare_atnum(catm(1,ldih(1,
     + i)))+argos_prepare_atnum(catm(1,ldih(4,i))).lt.
     + argos_prepare_atnum(catm(1,ldih(1,
     + j)))+argos_prepare_atnum(catm(1,ldih(4,j)))) then
      ldih(2,i)=ldih(2,j)
      ldih(4,i)=ldih(4,j)
      endif
      goto 20
      endif
   21 continue
      ndih=ndih+1
      ldih(1,ndih)=ldih(1,i)
      ldih(2,ndih)=ldih(2,i)
      ldih(3,ndih)=ldih(3,i)
      ldih(4,ndih)=ldih(4,i)
      ldih(5,ndih)=ldih(5,i)
      ldih(6,ndih)=ldih(6,i)
   20 continue
      endif
c
c     switch selected propers to impropers
c     ------------------------------------
c
      if(ffield(1:6).eq.'gromos'.or.ffield(1:6).eq.'charmm') then
      no=nimp
      do 22 i=nd+1,ndih
      do 23 j=1,nang
      if((ldih(1,i).eq.lang(1,j).and.ldih(4,i).eq.lang(3,j)).or.
     + (ldih(1,i).eq.lang(3,j).and.ldih(4,i).eq.lang(2,j))) then
      do 24 k=no+1,nimp
      if(limp(1,k).eq.i) goto 23
   24 continue
      nimp=nimp+1
      limp(1,nimp)=i
      endif
   23 continue
   22 continue
c
      do 25 i=nd+1,ndih
      do 26 j=1,i-1
      if((ldih(1,i).eq.ldih(1,j).and.ldih(4,i).eq.ldih(4,j)).or.
     + (ldih(1,i).eq.ldih(4,j).and.ldih(4,i).eq.ldih(4,j))) then
      do 27 k=no+1,nimp
      if(limp(1,k).eq.i) goto 28
   27 continue
      nimp=nimp+1
      limp(1,nimp)=i
   28 continue
      do 29 k=no+1,nimp
      if(limp(1,k).eq.j) goto 26
   29 continue
      nimp=nimp+1
      limp(1,nimp)=j
      endif
   26 continue
   25 continue
c
      do 30 i=no+1,nimp
      k=limp(1,i)
      do 31 j=1,4
      limp(j,i)=ldih(j,k)
   31 continue
      ldih(1,k)=0
   30 continue
c
      no=ndih
      ndih=0
      do 32 i=1,no
      if(ldih(1,i).gt.0) then
      ndih=ndih+1
      do 33 j=1,4
      ldih(j,ndih)=ldih(j,i)
   33 continue
      endif
   32 continue
c
      endif
c
      if(util_print('where',print_debug)) then
      write(lfnout,1107)
 1107 format('ARGOS_PREPARE_MKLIST START SORTING ')
      endif
c
c     sort the bond list
c     ------------------
c
      if(nbnd.gt.1) then
      do 34 i=1,nbnd
      if(lbnd(1,i).gt.lbnd(2,i)) then
      ltemp=lbnd(1,i)
      lbnd(1,i)=lbnd(2,i)
      lbnd(2,i)=ltemp
      endif
   34 continue
c
      do 35 i=1,nbnd
      ndx(i)=i
   35 continue
      ln=nbnd/2+1
      ir=nbnd
   36 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 37
      endif
      endif
      in=ln
      jn=ln+ln
   38 continue
      if(jn.le.ir) then
      if(jn.lt.ir) then
      if(lbnd(1,ndx(jn)).gt.lbnd(1,ndx(jn+1)).or.
     + (lbnd(1,ndx(jn)).eq.lbnd(1,ndx(jn+1)).and.
     + lbnd(2,ndx(jn)).gt.lbnd(2,ndx(jn+1)))) then
      else
      jn=jn+1
      endif
      endif
      if(lbnd(1,itmp).gt.lbnd(1,ndx(jn)).or.
     + (lbnd(1,itmp).eq.lbnd(1,ndx(jn)).and.
     + lbnd(2,itmp).gt.lbnd(2,ndx(jn)))) then
      jn=ir+1
      else
      ndx(in)=ndx(jn)
      in=jn
      jn=jn+jn
      endif
      goto 38
      endif
      ndx(in)=itmp
      goto 36
   37 continue
c
      do 39 k=1,4
      do 40 i=1,nbnd
      ival(i)=lbnd(k,i)
   40 continue
      do 41 i=1,nbnd
      lbnd(k,i)=ival(ndx(i))
   41 continue
   39 continue
      do 42 k=1,nparms
      do 442 l=1,2
      do 43 i=1,nbnd
      rval(i)=rbnd(k,l,i)
   43 continue
      do 44 i=1,nbnd
      rbnd(k,l,i)=rval(ndx(i))
   44 continue
  442 continue
   42 continue
      endif
c
      if(util_print('where',print_debug)) then
      write(lfnout,1121)
 1121 format('ARGOS_PREPARE_MKLIST BONDS SORTED ')
      endif
c
c     sort the angle lists
c     --------------------
c
      if(nang.gt.1) then
      do 45 i=1,nang
      if(lang(1,i).gt.lang(3,i)) then
      ltemp=lang(1,i)
      lang(1,i)=lang(3,i)
      lang(3,i)=ltemp
      endif
   45 continue
c
      do 46 i=1,nang
      ndx(i)=i
   46 continue
      ln=nang/2+1
      ir=nang
   47 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 48
      endif
      endif
      in=ln
      jn=ln+ln
   49 continue
      if(jn.le.ir) then
      if(jn.lt.ir) then
      if(lang(2,ndx(jn)).gt.lang(2,ndx(jn+1)).or.
     + (lang(2,ndx(jn)).eq.lang(2,ndx(jn+1)).and.
     + lang(1,ndx(jn)).gt.lang(1,ndx(jn+1))).or.
     + (lang(2,ndx(jn)).eq.lang(2,ndx(jn+1)).and.
     + lang(1,ndx(jn)).eq.lang(1,ndx(jn+1)).and.
     + lang(3,ndx(jn)).gt.lang(3,ndx(jn+1)))) then
      else
      jn=jn+1
      endif
      endif
      if(lang(2,itmp).gt.lang(2,ndx(jn)).or.
     + (lang(2,itmp).eq.lang(2,ndx(jn)).and.
     + lang(1,itmp).gt.lang(1,ndx(jn))).or.
     + (lang(2,itmp).eq.lang(2,ndx(jn)).and.
     + lang(1,itmp).eq.lang(1,ndx(jn)).and.
     + lang(3,itmp).gt.lang(3,ndx(jn)))) then
      jn=ir+1
      else
      ndx(in)=ndx(jn)
      in=jn
      jn=jn+jn
      endif
      goto 49
      endif
      ndx(in)=itmp
      goto 47
   48 continue
c
      do 50 k=1,5
      do 51 i=1,nang
      ival(i)=lang(k,i)
   51 continue
      do 52 i=1,nang
      lang(k,i)=ival(ndx(i))
   52 continue
   50 continue
      lto=2
      if(ffield(1:6).eq.'charmm') lto=4
      do 54 k=1,nparms
      do 554 l=1,lto
      do 55 i=1,nang
      rval(i)=rang(k,l,i)
   55 continue
      do 56 i=1,nang
      rang(k,l,i)=rval(ndx(i))
   56 continue
  554 continue
   54 continue
      endif
c
      if(util_print('where',print_debug)) then
      write(lfnout,1122)
 1122 format('ARGOS_PREPARE_MKLIST ANGLES SORTED ')
      endif
c
c     sort dihedral list
c     ------------------
c
      if(ndih.gt.1) then
      do 57 i=1,ndih
      if(ldih(2,i).gt.ldih(3,i)) then
      ltemp=ldih(2,i)
      ldih(2,i)=ldih(3,i)
      ldih(3,i)=ltemp
      ltemp=ldih(1,i)
      ldih(1,i)=ldih(4,i)
      ldih(4,i)=ltemp
      endif
   57 continue
c
      do 58 i=1,ndih
      ndx(i)=i
   58 continue
      ln=ndih/2+1
      ir=ndih
   59 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 60
      endif
      endif
      in=ln
      jn=ln+ln
   61 continue
      if(jn.le.ir) then
      if(jn.lt.ir) then
      if(ldih(2,ndx(jn)).gt.ldih(2,ndx(jn+1)).or.
     + (ldih(2,ndx(jn)).eq.ldih(2,ndx(jn+1)).and.
     + ldih(3,ndx(jn)).gt.ldih(3,ndx(jn+1))).or.
     + (ldih(2,ndx(jn)).eq.ldih(2,ndx(jn+1)).and.
     + ldih(3,ndx(jn)).eq.ldih(3,ndx(jn+1)).and.
     +  ldih(1,ndx(jn)).gt.ldih(1,ndx(jn+1))).or.
     + (ldih(2,ndx(jn)).eq.ldih(2,ndx(jn+1)).and.
     + ldih(3,ndx(jn)).eq.ldih(3,ndx(jn+1)).and.
     +  ldih(1,ndx(jn)).eq.ldih(1,ndx(jn+1)).and.
     + ldih(4,ndx(jn)).gt.ldih(4,ndx(jn+1)))) then
      else
      jn=jn+1
      endif
      endif
      if(ldih(2,itmp).gt.ldih(2,ndx(jn)).or.
     + (ldih(2,itmp).eq.ldih(2,ndx(jn)).and.
     + ldih(3,itmp).gt.ldih(3,ndx(jn))).or.
     + (ldih(2,itmp).eq.ldih(2,ndx(jn)).and.
     + ldih(3,itmp).eq.ldih(3,ndx(jn)).and.
     +  ldih(1,itmp).gt.ldih(1,ndx(jn))).or.
     + (ldih(2,itmp).eq.ldih(2,ndx(jn)).and.
     + ldih(3,itmp).eq.ldih(3,ndx(jn)).and.
     +  ldih(1,itmp).eq.ldih(1,ndx(jn)).and.
     + ldih(4,itmp).gt.ldih(4,ndx(jn)))) then
      jn=ir+1
      else
      ndx(in)=ndx(jn)
      in=jn
      jn=jn+jn
      endif
      goto 61
      endif
      ndx(in)=itmp
      goto 59
   60 continue
      do 62 k=1,6
      do 63 i=1,ndih
      ival(i)=ldih(k,i)
   63 continue
      do 64 i=1,ndih
      ldih(k,i)=ival(ndx(i))
   64 continue
   62 continue
      do 65 l=1,6
      do 66 k=1,nparms
      do 566 m=1,2
      do 67 i=1,ndih
      rval(i)=rdih(l,k,m,i)
   67 continue
      do 68 i=1,ndih
      rdih(l,k,m,i)=rval(ndx(i))
   68 continue
  566 continue
   66 continue
      do 69 k=1,nparms
      do 70 i=1,ndih
      ival(i)=kdih(l,k,i)
   70 continue
      do 71 i=1,ndih
      kdih(l,k,i)=ival(ndx(i))
   71 continue
   69 continue
   65 continue
      do 569 k=1,nparms
      do 570 i=1,ndih
      ival(i)=nfdih(k,i)
  570 continue
      do 571 i=1,ndih
      nfdih(k,i)=ival(ndx(i))
  571 continue
  569 continue
      endif
c
      if(util_print('where',print_debug)) then
      write(lfnout,1123)
 1123 format('ARGOS_PREPARE_MKLIST TORSIONS SORTED ')
      endif
c
c
c     sort improper list
c     ------------------
c
      if(nimp.gt.1) then
      do 72 i=1,nimp-1
      do 73 j=i+1,nimp
      if(limp(1,i).gt.limp(1,j).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).gt.limp(2,j)).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).eq.limp(2,j).and.
     + limp(3,i).gt.limp(3,j)).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).eq.limp(2,j).and.
     + limp(3,i).eq.limp(3,j).and.limp(4,i).gt.limp(4,j))) then
      do 74 k=1,6
      ltemp=limp(k,i)
      limp(k,i)=limp(k,j)
      limp(k,j)=ltemp
   74 continue
      do 574 k=1,nparms
      ltemp=kimp(k,i)
      kimp(k,i)=kimp(k,j)
      kimp(k,j)=ltemp
  574 continue
      do 75 k=1,nparms
      do 575 l=1,2
      rtemp=rimp(k,l,i)
      rimp(k,l,i)=rimp(k,l,j)
      rimp(k,l,j)=rtemp
  575 continue
   75 continue
      endif
   73 continue
   72 continue
      endif
c
      if(util_print('where',print_debug)) then
      write(lfnout,1108)
 1108 format('ARGOS_PREPARE_MKLIST DONE SORTING ')
      endif
cc
cc     check number of segments involved in angles
cc     -------------------------------------------
cc
c      do 76 i=1,nang
c      if(lseq(1,latm(5,lang(1,i))).ne.lseq(1,latm(5,lang(2,i))).and.
c     + lseq(1,latm(5,lang(1,i))).ne.lseq(1,latm(5,lang(3,i))).and.
c     + lseq(1,latm(5,lang(2,i))).ne.lseq(1,latm(5,lang(3,i)))) then
c      lang(5,i)=-1
c      endif
c   76 continue
cc
cc     check number of segments involved in dihedrals
cc     ----------------------------------------------
cc
c      do 77 i=1,ndih
c      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(2,i))).and.
c     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
c     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(3,i)))) then
c      ldih(6,i)=-1
c      else
c      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(2,i))).and.
c     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
c     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
c      ldih(6,i)=-1
c      else
c      if(lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
c     + lseq(1,latm(5,ldih(1,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
c     + lseq(1,latm(5,ldih(3,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
c      ldih(6,i)=-1
c      else
c      if(lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(3,i))).and.
c     + lseq(1,latm(5,ldih(2,i))).ne.lseq(1,latm(5,ldih(4,i))).and.
c     + lseq(1,latm(5,ldih(3,i))).ne.lseq(1,latm(5,ldih(4,i)))) then
c      ldih(6,i)=-1
c      endif
c      endif
c      endif
c      endif
c   77 continue
cc
cc     check number of segments involved in improper dihedrals
cc     -------------------------------------------------------
cc
c      do 78 i=1,nimp
c      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(2,i))).and.
c     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(3,i))).and.
c     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(3,i)))) then
c      limp(6,i)=-1
c      else
c      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(2,i))).and.
c     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(4,i))).and.
c     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(4,i)))) then
c      limp(6,i)=-1
c      else
c      if(lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(3,i))).and.
c     + lseq(1,latm(5,limp(1,i))).ne.lseq(1,latm(5,limp(4,i))).and.
c     + lseq(1,latm(5,limp(3,i))).ne.lseq(1,latm(5,limp(4,i)))) then
c      limp(6,i)=-1
c      else
c      if(lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(3,i))).and.
c     + lseq(1,latm(5,limp(2,i))).ne.lseq(1,latm(5,limp(4,i))).and.
c     + lseq(1,latm(5,limp(3,i))).ne.lseq(1,latm(5,limp(4,i)))) then
c      limp(6,i)=-1
c      endif
c      endif
c      endif
c      endif
c   78 continue
c
c     read solvent segment
c     --------------------
c
      natmt=natm
      nbndt=nbnd
      nangt=nang
      ndiht=ndih
      nimpt=nimp
      if(cseq(nseq+1)(1:1).ne.' ') then
      filmat=cseq(i)(1:index(cseq(nseq+1),' ')-1)//'.mat '
      lseq(3,nseq+1)=natm+1
      if(.not.argos_prepare_sgmfil(cseq(nseq+1),lfnsgm,filsgm))
     + call md_abort('argos_prepare_sgmfil failed',9999)
      if(.not.argos_prepare_rdsgm(lfnout,nparms,mparms,
     + 0,0,igrp,ipgrp,lfnsgm,
     + filsgm,lfnmat,filmat,latm,catm,qatm,matm,natmt,
     + lbnd,rbnd,mbnd,nbndt,lang,rang,mang,nangt,
     + ldih,nfdih,kdih,rdih,mdih,ndiht,limp,kimp,rimp,mimp,nimpt,
     + wcorr,nsgmpar,ipdef,itopol))
     + call md_abort('argos_prepare_rdsgm failed',9999)
      endif
      lseq(3,nseq+2)=natmt+1
c
      argos_prepare_mklist=.true.
c
      if(util_print('where',print_debug)) then
      write(lfnout,1110)
 1110 format('ARGOS_PREPARE_MKLIST DONE')
      endif
c
      return
      end
