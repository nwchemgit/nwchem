      logical function pre_dummy(lfnout,latt,catt,patt,ratt,matt,natt,
     + latm,catm,matm,natm,lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp)
c
c $Id: pre_dummy.F,v 1.1 2001-10-19 17:55:01 d3j191 Exp $
c
      implicit none
c
      integer lfnout
      integer matt,natt
      integer latt(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(9,matm)
      character*6 catm(4,matm)
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
c
      integer i,j,k
c
      do 1 i=1,nbnd
      k=-lbnd(3,i)
      if(k.eq.1) then
      if(catm(3,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(3,lbnd(2,i))(6:6).eq.'D') then
      rbnd(3,i)=rbnd(1,i)
      rbnd(4,i)=rbnd(2,i)
      k=k+2
      endif
      if(catm(4,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(4,lbnd(2,i))(6:6).eq.'D') then
      rbnd(5,i)=rbnd(1,i)
      rbnd(6,i)=rbnd(2,i)
      k=k+4
      endif
      elseif(k.eq.2) then
      if(catm(2,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(2,lbnd(2,i))(6:6).eq.'D') then
      rbnd(1,i)=rbnd(3,i)
      rbnd(2,i)=rbnd(4,i)
      k=k+1
      endif
      if(catm(4,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(4,lbnd(2,i))(6:6).eq.'D') then
      rbnd(5,i)=rbnd(3,i)
      rbnd(6,i)=rbnd(4,i)
      k=k+4
      endif
      elseif(k.eq.3) then
      if(catm(4,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(4,lbnd(2,i))(6:6).eq.'D') then
      rbnd(5,i)=rbnd(1,i)
      rbnd(6,i)=rbnd(2,i)
      k=k+4
      endif
      elseif(k.eq.4) then
      if(catm(2,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(2,lbnd(2,i))(6:6).eq.'D') then
      rbnd(1,i)=rbnd(5,i)
      rbnd(2,i)=rbnd(6,i)
      k=k+1
      endif
      if(catm(3,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(3,lbnd(2,i))(6:6).eq.'D') then
      rbnd(3,i)=rbnd(5,i)
      rbnd(4,i)=rbnd(6,i)
      k=k+2
      endif
      elseif(k.eq.5) then
      if(catm(3,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(3,lbnd(2,i))(6:6).eq.'D') then
      rbnd(3,i)=rbnd(1,i)
      rbnd(4,i)=rbnd(2,i)
      k=k+2
      endif
      elseif(k.eq.6) then
      if(catm(2,lbnd(1,i))(6:6).eq.'D'.or.
     +   catm(2,lbnd(2,i))(6:6).eq.'D') then
      rbnd(1,i)=rbnd(3,i)
      rbnd(2,i)=rbnd(4,i)
      k=k+1
      endif
      endif
      lbnd(3,i)=-k
    1 continue
c
      do 2 i=1,nang
      k=-lang(4,i)
      if(k.eq.1) then
      if(catm(3,lang(1,i))(6:6).eq.'D'.or.
     +   catm(3,lang(2,i))(6:6).eq.'D'.or.
     +   catm(3,lang(3,i))(6:6).eq.'D') then
      rang(3,i)=rang(1,i)
      rang(4,i)=rang(2,i)
      k=k+2
      endif
      if(catm(4,lang(1,i))(6:6).eq.'D'.or.
     +   catm(4,lang(2,i))(6:6).eq.'D'.or.
     +   catm(4,lang(3,i))(6:6).eq.'D') then
      rang(5,i)=rang(1,i)
      rang(6,i)=rang(2,i)
      k=k+4
      endif
      elseif(k.eq.2) then
      if(catm(2,lang(1,i))(6:6).eq.'D'.or.
     +   catm(2,lang(2,i))(6:6).eq.'D'.or.
     +   catm(2,lang(3,i))(6:6).eq.'D') then
      rang(1,i)=rang(3,i)
      rang(2,i)=rang(4,i)
      k=k+1
      endif
      if(catm(4,lang(1,i))(6:6).eq.'D'.or.
     +   catm(4,lang(2,i))(6:6).eq.'D'.or.
     +   catm(4,lang(3,i))(6:6).eq.'D') then
      rang(5,i)=rang(3,i)
      rang(6,i)=rang(4,i)
      k=k+4
      endif
      elseif(k.eq.3) then
      if(catm(4,lang(1,i))(6:6).eq.'D'.or.
     +   catm(4,lang(2,i))(6:6).eq.'D'.or.
     +   catm(4,lang(3,i))(6:6).eq.'D') then
      rang(5,i)=rang(1,i)
      rang(6,i)=rang(2,i)
      k=k+4
      endif
      elseif(k.eq.4) then
      if(catm(2,lang(1,i))(6:6).eq.'D'.or.
     +   catm(2,lang(2,i))(6:6).eq.'D'.or.
     +   catm(2,lang(3,i))(6:6).eq.'D') then
      rang(1,i)=rang(5,i)
      rang(2,i)=rang(6,i)
      k=k+1
      endif
      if(catm(3,lang(1,i))(6:6).eq.'D'.or.
     +   catm(3,lang(2,i))(6:6).eq.'D'.or.
     +   catm(3,lang(3,i))(6:6).eq.'D') then
      rang(3,i)=rang(5,i)
      rang(4,i)=rang(6,i)
      k=k+2
      endif
      elseif(k.eq.5) then
      if(catm(3,lang(1,i))(6:6).eq.'D'.or.
     +   catm(3,lang(2,i))(6:6).eq.'D'.or.
     +   catm(3,lang(3,i))(6:6).eq.'D') then
      rang(3,i)=rang(1,i)
      rang(4,i)=rang(2,i)
      k=k+2
      endif
      elseif(k.eq.6) then
      if(catm(2,lang(1,i))(6:6).eq.'D'.or.
     +   catm(2,lang(2,i))(6:6).eq.'D'.or.
     +   catm(2,lang(3,i))(6:6).eq.'D') then
      rang(1,i)=rang(3,i)
      rang(2,i)=rang(4,i)
      k=k+1
      endif
      endif
      lang(4,i)=-k
    2 continue
c
      do 3 i=1,ndih
      k=-ldih(5,i)
      if(k.eq.1) then
      if(catm(3,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(4,i))(6:6).eq.'D') then
      ldih(8,i)=ldih(7,i)
      kdih(ldih(8,i),2,i)=kdih(ldih(7,i),1,i)
      rdih(ldih(8,i),3,i)=rdih(ldih(7,i),1,i)
      rdih(ldih(8,i),4,i)=rdih(ldih(7,i),2,i)
      k=k+2
      endif
      if(catm(4,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(4,i))(6:6).eq.'D') then
      ldih(9,i)=ldih(7,i)
      kdih(ldih(9,i),3,i)=kdih(ldih(7,i),1,i)
      rdih(ldih(9,i),5,i)=rdih(ldih(7,i),1,i)
      rdih(ldih(9,i),6,i)=rdih(ldih(7,i),2,i)
      k=k+4
      endif
      elseif(k.eq.2) then
      if(catm(2,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(4,i))(6:6).eq.'D') then
      ldih(7,i)=ldih(8,i)
      kdih(ldih(7,i),1,i)=kdih(ldih(8,i),2,i)
      rdih(ldih(7,i),1,i)=rdih(ldih(8,i),3,i)
      rdih(ldih(7,i),2,i)=rdih(ldih(8,i),4,i)
      k=k+1
      endif
      if(catm(4,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(4,i))(6:6).eq.'D') then
      ldih(9,i)=ldih(8,i)
      kdih(ldih(9,i),3,i)=kdih(ldih(8,i),2,i)
      rdih(ldih(9,i),5,i)=rdih(ldih(8,i),3,i)
      rdih(ldih(9,i),6,i)=rdih(ldih(8,i),4,i)
      k=k+4
      endif
      elseif(k.eq.3) then
      if(catm(4,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(4,ldih(4,i))(6:6).eq.'D') then
      ldih(9,i)=ldih(7,i)
      kdih(ldih(9,i),3,i)=kdih(ldih(7,i),1,i)
      rdih(ldih(9,i),5,i)=rdih(ldih(7,i),1,i)
      rdih(ldih(9,i),6,i)=rdih(ldih(7,i),2,i)
      k=k+4
      endif
      elseif(k.eq.4) then
      if(catm(2,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(4,i))(6:6).eq.'D') then
      ldih(7,i)=ldih(9,i)
      kdih(ldih(7,i),1,i)=kdih(ldih(9,i),3,i)
      rdih(ldih(7,i),1,i)=rdih(ldih(9,i),5,i)
      rdih(ldih(7,i),2,i)=rdih(ldih(9,i),6,i)
      k=k+1
      endif
      if(catm(3,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(4,i))(6:6).eq.'D') then
      ldih(8,i)=ldih(9,i)
      kdih(ldih(8,i),2,i)=kdih(ldih(9,i),3,i)
      rdih(ldih(8,i),3,i)=rdih(ldih(9,i),5,i)
      rdih(ldih(8,i),4,i)=rdih(ldih(9,i),6,i)
      k=k+2
      endif
      elseif(k.eq.5) then
      if(catm(3,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(3,ldih(4,i))(6:6).eq.'D') then
      ldih(8,i)=ldih(7,i)
      kdih(ldih(8,i),2,i)=kdih(ldih(7,i),1,i)
      rdih(ldih(8,i),3,i)=rdih(ldih(7,i),1,i)
      rdih(ldih(8,i),4,i)=rdih(ldih(7,i),2,i)
      k=k+2
      endif
      elseif(k.eq.6) then
      if(catm(2,ldih(1,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(2,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(3,i))(6:6).eq.'D'.or.
     +   catm(2,ldih(4,i))(6:6).eq.'D') then
      ldih(7,i)=ldih(8,i)
      kdih(ldih(7,i),1,i)=kdih(ldih(8,i),2,i)
      rdih(ldih(7,i),1,i)=rdih(ldih(8,i),3,i)
      rdih(ldih(7,i),2,i)=rdih(ldih(8,i),4,i)
      k=k+1
      endif
      endif
      ldih(5,i)=-k
    3 continue
c
      do 4 i=1,nimp
      k=-limp(5,i)
      if(k.eq.1) then
      if(catm(3,limp(1,i))(6:6).eq.'D'.or.
     +   catm(3,limp(2,i))(6:6).eq.'D'.or.
     +   catm(3,limp(3,i))(6:6).eq.'D'.or.
     +   catm(3,limp(4,i))(6:6).eq.'D') then
      limp(8,i)=limp(7,i)
      rimp(3,i)=rimp(1,i)
      rimp(4,i)=rimp(2,i)
      k=k+2
      endif
      if(catm(4,limp(1,i))(6:6).eq.'D'.or.
     +   catm(4,limp(2,i))(6:6).eq.'D'.or.
     +   catm(4,limp(3,i))(6:6).eq.'D'.or.
     +   catm(4,limp(4,i))(6:6).eq.'D') then
      limp(9,i)=limp(7,i)
      rimp(5,i)=rimp(1,i)
      rimp(6,i)=rimp(2,i)
      k=k+4
      endif
      elseif(k.eq.2) then
      if(catm(2,limp(1,i))(6:6).eq.'D'.or.
     +   catm(2,limp(2,i))(6:6).eq.'D'.or.
     +   catm(2,limp(3,i))(6:6).eq.'D'.or.
     +   catm(2,limp(4,i))(6:6).eq.'D') then
      limp(7,i)=limp(8,i)
      rimp(1,i)=rimp(3,i)
      rimp(2,i)=rimp(4,i)
      k=k+1
      endif
      if(catm(4,limp(1,i))(6:6).eq.'D'.or.
     +   catm(4,limp(2,i))(6:6).eq.'D'.or.
     +   catm(4,limp(3,i))(6:6).eq.'D'.or.
     +   catm(4,limp(4,i))(6:6).eq.'D') then
      limp(9,i)=limp(8,i)
      rimp(5,i)=rimp(3,i)
      rimp(6,i)=rimp(4,i)
      k=k+4
      endif
      elseif(k.eq.3) then
      if(catm(4,limp(1,i))(6:6).eq.'D'.or.
     +   catm(4,limp(2,i))(6:6).eq.'D'.or.
     +   catm(4,limp(3,i))(6:6).eq.'D'.or.
     +   catm(4,limp(4,i))(6:6).eq.'D') then
      limp(9,i)=limp(7,i)
      rimp(5,i)=rimp(1,i)
      rimp(6,i)=rimp(2,i)
      k=k+4
      endif
      elseif(k.eq.4) then
      if(catm(2,limp(1,i))(6:6).eq.'D'.or.
     +   catm(2,limp(2,i))(6:6).eq.'D'.or.
     +   catm(2,limp(3,i))(6:6).eq.'D'.or.
     +   catm(2,limp(4,i))(6:6).eq.'D') then
      limp(7,i)=limp(9,i)
      rimp(1,i)=rimp(5,i)
      rimp(2,i)=rimp(6,i)
      k=k+1
      endif
      if(catm(3,limp(1,i))(6:6).eq.'D'.or.
     +   catm(3,limp(2,i))(6:6).eq.'D'.or.
     +   catm(3,limp(3,i))(6:6).eq.'D'.or.
     +   catm(3,limp(4,i))(6:6).eq.'D') then
      limp(8,i)=limp(9,i)
      rimp(3,i)=rimp(5,i)
      rimp(4,i)=rimp(6,i)
      k=k+2
      endif
      elseif(k.eq.5) then
      if(catm(3,limp(1,i))(6:6).eq.'D'.or.
     +   catm(3,limp(2,i))(6:6).eq.'D'.or.
     +   catm(3,limp(3,i))(6:6).eq.'D'.or.
     +   catm(3,limp(4,i))(6:6).eq.'D') then
      limp(8,i)=limp(7,i)
      rimp(3,i)=rimp(1,i)
      rimp(4,i)=rimp(2,i)
      k=k+2
      endif
      elseif(k.eq.6) then
      if(catm(2,limp(1,i))(6:6).eq.'D'.or.
     +   catm(2,limp(2,i))(6:6).eq.'D'.or.
     +   catm(2,limp(3,i))(6:6).eq.'D'.or.
     +   catm(2,limp(4,i))(6:6).eq.'D') then
      limp(7,i)=limp(8,i)
      rimp(1,i)=rimp(3,i)
      rimp(2,i)=rimp(4,i)
      k=k+1
      endif
      endif
      limp(5,i)=-k
    4 continue
c
      pre_dummy=.true.
      return
 9999 continue
      pre_dummy=.false.
      return
      end
