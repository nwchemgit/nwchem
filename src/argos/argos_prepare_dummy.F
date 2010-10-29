      logical function argos_prepare_dummy(lfnout,latt,catt,patt,ratt,
     + matt,natt,
     + latm,catm,matm,natm,lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,limp,kimp,rimp,mimp,nimp,
     + nparms,mparms)
c
c $Id$
c
      implicit none
c
#include "argos_prepare_common.fh"
c
      integer lfnout
      integer matt,natt,nparms,mparms
      integer latt(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(11,matm)
      character*6 catm(mparms,matm)
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
c
      integer i,j,k,l,kpar
c
      do 1 i=1,nbnd
c      write(*,'(3i5)') (lbnd(j,i),j=1,3)
      k=-lbnd(3,i)
c      if(k.eq.0) call md_abort('Error in argos_prepare_dummy',1)
      if(k.ne.0) then
      kpar=0
      do 2 j=nparms,1,-1
      if(iand(k,2**(j-1)).ne.0) kpar=j
    2 continue
      if(kpar.eq.0) call md_abort('Error in argos_prepare_dummy',2)
      do 3 j=1,nparms
      if(iand(k,2**(j-1)).eq.0) then
      if(catm(j+1,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(j+1,lbnd(2,i))(6:6).eq.'D') then
      rbnd(j,1,i)=rbnd(kpar,1,i)
      rbnd(j,2,i)=rbnd(kpar,2,i)
      k=ior(k,2**(j-1))
      endif
      endif
    3 continue
      lbnd(3,i)=-k
      endif
    1 continue
c
      do 4 i=1,nang
      k=-lang(4,i)
c      print*,i,k
c      write(*,'(4i5)') (lang(j,i),j=1,4)
c      if(k.eq.0) call md_abort('Error in argos_prepare_dummy',3)
      if(k.gt.0) then
      kpar=0
      do 5 j=nparms,1,-1
      if(iand(k,2**(j-1)).ne.0) kpar=j
    5 continue
      if(kpar.eq.0) call md_abort('Error in argos_prepare_dummy',4)
      do 6 j=1,nparms
c      if(iand(k,2**(j-1)).eq.0) then
      if(catm(j+1,lang(1,i))(6:6).eq.'D'.or.
     +   catm(j+1,lang(2,i))(6:6).eq.'D'.or.
     +   catm(j+1,lang(3,i))(6:6).eq.'D') then
      rang(j,1,i)=rang(kpar,1,i)
      rang(j,2,i)=0.0d0
      if(ffield(1:6).eq.'charmm') then
      rang(j,3,i)=rang(kpar,3,i)
      rang(j,4,i)=0.0d0
      endif
c      print*,'Dummy',i,j
cx      rang(j,2,i)=rang(kpar,2,i)
      k=ior(k,2**(j-1))
      endif
c      endif
    6 continue
      lang(4,i)=-k
      endif
    4 continue
c
      do 7 i=1,ndih
      k=-ldih(5,i)
c      if(k.eq.0) call md_abort('Error in argos_prepare_dummy',5)
      if(k.gt.0) then
      kpar=0
      do 8 j=nparms,1,-1
      if(iand(k,2**(j-1)).ne.0) kpar=j
    8 continue
      if(kpar.eq.0) call md_abort('Error in argos_prepare_dummy',6)
      do 9 j=1,nparms
      if(iand(k,2**(j-1)).eq.0) then
      if(catm(j+1,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(j+1,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(j+1,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(j+1,ldih(4,i))(6:6).eq.'D') then
      nfdih(j,i)=nfdih(kpar,i)
      do 10 l=1,nfdih(j,i)
      kdih(l,j,i)=kdih(l,kpar,i)
      rdih(l,j,1,i)=rdih(l,kpar,1,i)
      rdih(l,j,2,i)=rdih(l,kpar,2,i)
   10 continue
      k=ior(k,2**(j-1))
      endif
      endif
    9 continue
      ldih(5,i)=-k
      endif
    7 continue
c
      do 11 i=1,nimp
      k=-limp(5,i)
c      if(k.eq.0) call md_abort('Error in argos_prepare_dummy',7)
      if(k.gt.0) then
      kpar=0
      do 12 j=nparms,1,-1
      if(iand(k,2**(j-1)).ne.0) kpar=j
   12 continue
      if(kpar.eq.0) call md_abort('Error in argos_prepare_dummy',8)
      do 13 j=1,nparms
      if(iand(k,2**(j-1)).eq.0) then
      if(catm(j+1,limp(1,i))(6:6).eq.'D'.or.
     +   catm(j+1,limp(2,i))(6:6).eq.'D'.or.
     +   catm(j+1,limp(3,i))(6:6).eq.'D'.or.
     +   catm(j+1,limp(4,i))(6:6).eq.'D') then
      kimp(j,i)=kimp(kpar,i)
      rimp(j,1,i)=rimp(kpar,1,i)
      rimp(j,2,i)=rimp(kpar,2,i)
      k=ior(k,2**(j-1))
      endif
      endif
   13 continue
      limp(5,i)=-k
      endif
   11 continue
c
      argos_prepare_dummy=.true.
      return
 9999 continue
      argos_prepare_dummy=.false.
      return
      end
