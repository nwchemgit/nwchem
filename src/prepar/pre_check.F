      logical function pre_check(lfnout,lprint,
     + latt,catt,matt,natt,catm,matm,natm,lbnd,mbnd,nbnd,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp)
c
      implicit none
c
      integer lfnout,lprint,matm,natm
      integer matt,natt,mbnd,nbnd,mang,nang,mdih,ndih,mimp,nimp
      integer latt(3,matt),lbnd(4,mbnd),lang(5,mang),ldih(9,mdih)
      integer limp(9,mimp)
      character*6 catt(2,matt),catm(4,matm)
c
      integer i,j,k
c
      pre_check=.true.
c
c     check the atom types
c     --------------------
c
      do 1 i=1,natt
      if(latt(1,i).le.0) then
      if(lprint.ge.1) write(lfnout,1000) catt(1,i)
 1000 format('PARAMETERS COULD NOT BE FOUND FOR ATOM TYPE ',a6)
      pre_check=.false.
      endif
    1 continue
c
      do 2 i=1,nbnd
      if(lbnd(3,i).le.0.and.lbnd(3,i).gt.-7) then
      k=-lbnd(3,i)
      do 3 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 3
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 3
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 3
      if(lprint.ge.1) write(lfnout,1001)
     + catm(1+j,lbnd(1,i)),catm(1+j,lbnd(2,i))
 1001 format('PARAMETERS COULD NOT BE FOUND FOR BOND TYPE ',a6,'-',a6)
    3 continue
      pre_check=.false. 
      endif
    2 continue
c
      do 4 i=1,nang
      if(lang(4,i).le.0.and.lang(4,i).gt.-7) then
      k=-lang(4,i)
      do 5 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 3
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 3
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 3
      if(lprint.ge.1) write(lfnout,1002)
     + catm(1+j,lang(1,i)),catm(1+j,lang(2,i)),catm(1+j,lang(3,i))
 1002 format('PARAMETERS COULD NOT BE FOUND FOR ANGLE TYPE ',
     + a6,'-',a6,'-',a6)
    5 continue
      pre_check=.false. 
      endif
    4 continue
c
      do 6 i=1,ndih
      if(ldih(5,i).le.0.and.ldih(5,i).gt.-7) then
      k=-ldih(5,i)
      do 7 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 3
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 3
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 3
      if(lprint.ge.1) write(lfnout,1003) catm(1+j,ldih(1,i)),
     + catm(1+j,ldih(2,i)),catm(1+j,ldih(3,i)),catm(1+j,ldih(4,i))
 1003 format('PARAMETERS COULD NOT BE FOUND FOR DIHEDRAL TYPE ',
     + a6,'-',a6,'-',a6,'-',a6)
    7 continue
      pre_check=.false. 
      endif
    6 continue
c
      do 8 i=1,nimp
      if(limp(5,i).le.0.and.limp(5,i).gt.-7) then
      k=-limp(5,i)
      do 9 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 3
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 3
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 3
      if(lprint.ge.1) write(lfnout,1004) catm(1+j,limp(1,i)),
     + catm(1+j,limp(2,i)),catm(1+j,limp(3,i)),catm(1+j,limp(4,i))
 1004 format('PARAMETERS COULD NOT BE FOUND FOR IMPROPER TYPE ',
     + a6,'-',a6,'-',a6,'-',a6)
    9 continue
      pre_check=.false. 
      endif
    8 continue
c
      return
      end
