      logical function pre_mklist(lfnout,lprint,lfnsgm,
     + ffield,dir_s,dir_x,dir_u,dir_t,
     + lseq,cseq,mseq,nseq,
     + llnk,clnk,mlnk,nlnk,
     + latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,
     + lang,rang,mang,nang,
     + ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp,lmat,
     + natmt,nbndt,nangt,ndiht,nimpt,wcorr)
c
      implicit none
c
      logical pre_sgmfil,pre_rdsgm
      integer pre_atnum
      external pre_sgmfil,pre_rdsgm,pre_atnum
c
      integer lfnout,lfnsgm
      integer lprint
      character*80 ffield
      character*255 dir_s,dir_x,dir_u,dir_t
      integer mseq,nseq
      integer lseq(4,mseq),lmat(10,mseq)
      character*10 cseq(mseq)
      integer mlnk,nlnk
      integer llnk(2,mlnk)
      character*4 clnk(2,mlnk)
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
      integer natmt,nbndt,nangt,ndiht,nimpt
      real*8 wcorr(3)
c
      integer i,j,k,l,link1,link2,isgm,jsgm,nb,na,nd,ni,nn,iatm,ia(10)
      integer no,ndi,ltemp,igrp,ipgrp,nia(3)
      real*8 rtemp
      character*255 filsgm
c
      igrp=0
      ipgrp=0
      do 1 i=1,nseq
c
      lseq(3,i)=natm+1
c
c     generate segment filename
c     -------------------------
c
      if(.not.pre_sgmfil(cseq(i),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_sgmfil failed',9999)
c
c     append the atom and bonded lists
c     --------------------------------
c
      if(.not.pre_rdsgm(lfnout,lprint,lseq(4,i),i,igrp,ipgrp,
     + lfnsgm,filsgm,
     + latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp,wcorr))
     + call errquit('pre_rdsgm failed',9999)
c
    1 continue
      lseq(3,nseq+1)=natm+1
c
c     add links
c     ---------
c
      do 2 i=1,nseq
      do 3 j=1,10
      lmat(j,i)=0
    3 continue
    2 continue
c
      do 5 i=1,nseq-1
      if((lseq(2,i).eq.-1.or.lseq(2,i).eq.0).and.lseq(2,i+1).ge.0) then
      link2=0
      do 6 j=lseq(3,i),lseq(3,i+1)-1
      if(latm(3,j).eq.2) link2=j
    6 continue
      link1=0
      do 7 j=lseq(3,i+1),lseq(3,i+2)-1
      if(latm(3,j).eq.1) link1=j
    7 continue
      if(link1.gt.0.and.link2.gt.0) then
      lmat(2,i)=link1
      lmat(1,i+1)=link2
      endif
      endif
    5 continue
c
      do 8 i=1,nlnk
      isgm=0
      jsgm=0
      do 9 j=1,nseq
      if(lseq(1,j).eq.llnk(1,i)) isgm=j
      if(lseq(1,j).eq.llnk(2,i)) jsgm=j
    9 continue
      if(isgm.eq.0.or.jsgm.eq.0)
     + call errquit('error 1 in pre_mklist',9999)
      link1=0
      do 10 j=lseq(3,isgm),lseq(3,isgm+1)-1
      if(catm(1,j)(1:4).eq.clnk(1,i)) link1=j
   10 continue
      link2=0
      do 11 j=lseq(3,jsgm),lseq(3,jsgm+1)-1
      if(catm(1,j)(1:4).eq.clnk(2,i)) link2=j
   11 continue
      if(link1.eq.0.or.link2.eq.0)
     + call errquit('error 2 in pre_mklist',9999)
      if(latm(3,link1).eq.0)
     + call errquit('error 3 in pre_mklist',9999)
      if(latm(3,link2).eq.0)
     + call errquit('error 4 in pre_mklist',9999)
      lmat(latm(3,link1),isgm)=link2
      lmat(latm(3,link2),jsgm)=link1
    8 continue
c
      nb=nbnd
      do 12 i=1,natm
      if(latm(3,i).gt.0) then
      if(lmat(latm(3,i),latm(5,i)).eq.0)
     + call errquit('error 5 in pre_mklist',9999)
      if(i.lt.lmat(latm(3,i),latm(5,i))) then
      nbnd=nbnd+1
      if(nbnd.gt.mbnd) call errquit('increase mbnd',9999)
      lbnd(1,nbnd)=i
      lbnd(2,nbnd)=lmat(latm(3,i),latm(5,i))
      lbnd(3,nbnd)=0
      lbnd(4,nbnd)=0
      endif
      endif
   12 continue
c
c     append the angle list
c     ---------------------
c
      na=nang
      do 13 i=1,nb
      do 14 j=nb+1,nbnd
      if(lbnd(1,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) call errquit('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(2,j)
      lang(4,nang)=0
      else
      if(lbnd(2,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) call errquit('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(2,j)
      lang(4,nang)=0
      else
      if(lbnd(1,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) call errquit('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(1,j)
      lang(4,nang)=0
      else
      if(lbnd(2,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) call errquit('increase mang',9999)
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(1,j)
      lang(4,nang)=0
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
      do 15 i=1,nang
      do 16 j=na+1,nang
      if(i.gt.na.and.j.le.i) goto 16
      if(lang(2,i).eq.lang(1,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) call errquit('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(3,j)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      ldih(7,ndih)=0
      ldih(8,ndih)=0
      ldih(9,ndih)=0
      else
      if(lang(2,i).eq.lang(3,j).and.lang(3,i).eq.lang(2,j)) then
      if(ndih.ge.mdih) call errquit('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,i)
      ldih(2,ndih)=lang(2,i)
      ldih(3,ndih)=lang(3,i)
      ldih(4,ndih)=lang(1,j)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      ldih(7,ndih)=0
      ldih(8,ndih)=0
      ldih(9,ndih)=0
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(3,j)) then
      if(ndih.ge.mdih) call errquit('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(1,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      ldih(7,ndih)=0
      ldih(8,ndih)=0
      ldih(9,ndih)=0
      else
      if(lang(1,i).eq.lang(2,j).and.lang(2,i).eq.lang(1,j)) then
      if(ndih.ge.mdih) call errquit('increase mdih',9999)
      ndih=ndih+1
      ldih(1,ndih)=lang(3,j)
      ldih(2,ndih)=lang(1,i)
      ldih(3,ndih)=lang(2,i)
      ldih(4,ndih)=lang(3,i)
      ldih(5,ndih)=0
      ldih(6,ndih)=0
      ldih(7,ndih)=0
      ldih(8,ndih)=0
      ldih(9,ndih)=0
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
      do 54 k=1,nbnd
      if(lbnd(1,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(2,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(1,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(2,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(1,k).eq.ia(3)) nia(3)=nia(3)+1
      if(lbnd(2,k).eq.ia(3)) nia(3)=nia(3)+1
   54 continue
      do 55 k=1,2
      do 56 l=2,3
      if(nia(k).lt.nia(l)) then
      ltemp=nia(k)
      nia(k)=nia(l)
      nia(l)=ltemp
      ltemp=ia(k)
      ia(k)=ia(l)
      ia(l)=ltemp
      endif
   56 continue
   55 continue
c
c     planar
c
      if(latm(4,iatm).eq.1.and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call errquit('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
      endif
c
c     chiral
c
      if((latm(4,iatm).eq.2.or.latm(4,iatm).eq.3).and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call errquit('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
      endif
c
c     aromatic ring
c
      if(ffield(1:5).eq.'amber'.or.ffield(1:5).eq.'charm') then
      if(latm(4,iatm).eq.4.and.nn.eq.3) then
      nimp=nimp+1
      if(nimp.ge.mimp) call errquit('increase mimp',9999)
      limp(1,nimp)=iatm
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      limp(5,nimp)=0
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
      if(pre_atnum(catm(1,ldih(1,i)))+pre_atnum(catm(1,ldih(4,i))).lt.
     + pre_atnum(catm(1,ldih(1,j)))+pre_atnum(catm(1,ldih(4,j)))) then
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
      ldih(7,ndih)=ldih(7,i)
      ldih(8,ndih)=ldih(8,i)
      ldih(9,ndih)=ldih(9,i)
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
c     sort the bond list
c     ------------------
c
      do 34 i=1,nbnd
      if(lbnd(1,i).gt.lbnd(2,i)) then
      ltemp=lbnd(1,i)
      lbnd(1,i)=lbnd(2,i)
      lbnd(2,i)=ltemp
      endif
   34 continue
      do 35 i=1,nbnd-1
      do 36 j=i+1,nbnd
      if(lbnd(1,i).gt.lbnd(1,j).or.
     + (lbnd(1,i).eq.lbnd(1,j).and.lbnd(2,i).gt.lbnd(2,j))) then
      do 37 k=1,4
      ltemp=lbnd(k,i)
      lbnd(k,i)=lbnd(k,j)
      lbnd(k,j)=ltemp
   37 continue
      do 38 k=1,6
      rtemp=rbnd(k,i)
      rbnd(k,i)=rbnd(k,j)
      rbnd(k,j)=rtemp
   38 continue
      endif
   36 continue
   35 continue
c
c     sort the angle lists
c     --------------------
c
      do 39 i=1,nang
      if(lang(1,i).gt.lang(3,i)) then
      ltemp=lang(1,i)
      lang(1,i)=lang(3,i)
      lang(3,i)=ltemp
      endif
   39 continue
c
      do 40 i=1,nang-1
      do 41 j=i+1,nang
      if(lang(2,i).gt.lang(2,j).or.
     + (lang(2,i).eq.lang(2,j).and.lang(1,i).gt.lang(1,j)).or.
     + (lang(2,i).eq.lang(2,j).and.lang(1,i).eq.lang(1,j).and.
     + lang(3,i).gt.lang(3,j))) then
      do 42 k=1,5
      ltemp=lang(k,i)
      lang(k,i)=lang(k,j)
      lang(k,j)=ltemp
   42 continue
      do 43 k=1,6
      rtemp=rang(k,i)
      rang(k,i)=rang(k,j)
      rang(k,j)=rtemp
   43 continue
      endif
   41 continue
   40 continue
c
c     sort dihedral list
c     ------------------
c
      do 44 i=1,ndih
      if(ldih(2,i).gt.ldih(3,i)) then
      ltemp=ldih(2,i)
      ldih(2,i)=ldih(3,i)
      ldih(3,i)=ltemp
      ltemp=ldih(1,i)
      ldih(1,i)=ldih(4,i)
      ldih(4,i)=ltemp
      endif
   44 continue
c
      do 45 i=1,ndih-1
      do 46 j=i+1,ndih
      if(ldih(2,i).gt.ldih(2,j).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).gt.ldih(3,j)).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).eq.ldih(3,j).and.
     +  ldih(1,i).gt.ldih(1,j)).or.
     + (ldih(2,i).eq.ldih(2,j).and.ldih(3,i).eq.ldih(3,j).and.
     +  ldih(1,i).eq.ldih(1,j).and.ldih(4,i).gt.ldih(4,j))) then
      do 47 k=1,9
      ltemp=ldih(k,i)
      ldih(k,i)=ldih(k,j)
      ldih(k,j)=ltemp
   47 continue
      do 147 l=1,6
      do 48 k=1,6
      rtemp=rdih(l,k,i)
      rdih(l,k,i)=rdih(l,k,j)
      rdih(l,k,j)=rtemp
   48 continue
      do 148 k=1,3
      ltemp=kdih(l,k,i)
      kdih(l,k,i)=kdih(l,k,j)
      kdih(l,k,j)=ltemp
  148 continue
  147 continue
      endif
   46 continue
   45 continue
c
c     sort improper list
c     ------------------
c
      do 49 i=1,nimp-1
      do 50 j=i+1,nimp
      if(limp(1,i).gt.limp(1,j).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).gt.limp(2,j)).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).eq.limp(2,j).and.
     + limp(3,i).gt.limp(3,j)).or.
     + (limp(1,i).eq.limp(1,j).and.limp(2,i).eq.limp(2,j).and.
     + limp(3,i).eq.limp(3,j).and.limp(4,i).gt.limp(4,j))) then
      do 51 k=1,9
      ltemp=limp(k,i)
      limp(k,i)=limp(k,j)
      limp(k,j)=ltemp
   51 continue
      do 52 k=1,6
      rtemp=rimp(k,i)
      rimp(k,i)=rimp(k,j)
      rimp(k,j)=rtemp
   52 continue
      endif
   50 continue
   49 continue
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
      lseq(3,nseq+1)=natm+1
      if(.not.pre_sgmfil(cseq(nseq+1),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_sgmfil failed',9999)
      if(.not.pre_rdsgm(lfnout,lprint,0,0,igrp,ipgrp,lfnsgm,filsgm,
     + latm,catm,qatm,matm,natmt,
     + lbnd,rbnd,mbnd,nbndt,lang,rang,mang,nangt,
     + ldih,kdih,rdih,mdih,ndiht,limp,rimp,mimp,nimpt,wcorr))
     + call errquit('pre_rdsgm failed',9999)
      endif
      lseq(3,nseq+2)=natmt+1
c
      pre_mklist=.true.
      return
      end
