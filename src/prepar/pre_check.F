      logical function pre_check(lfnout,ffield,
     + latt,catt,matt,natt,catm,matm,natm,lbnd,mbnd,nbnd,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp)
c
c $Id: pre_check.F,v 1.10 1999-05-19 20:16:15 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,matm,natm
      integer matt,natt,mbnd,nbnd,mang,nang,mdih,ndih,mimp,nimp
      integer latt(3,matt),lbnd(4,mbnd),lang(5,mang),ldih(9,mdih)
      integer limp(9,mimp)
      character*6 catt(2,matt),catm(4,matm)
      character*80 ffield
c
      integer i,j,k,ii,jj,kk
c
      pre_check=.true.
c
c     check the atom types
c     --------------------
c
 1000 format(/,'Undetermined force field parameters',/)
c
      do 1 i=1,natt
      if(latt(1,i).le.0) then
      do 11 j=1,i-1
      if(latt(1,j).le.0.and.catt(1,i).eq.catt(1,j)) goto 1
   11 continue
      if(util_print('topology',print_none)) then
      if(pre_check) write(lfnout,1000)
      write(lfnout,1001) catt(1,i)
 1001 format('Parameters could not be found for atom type     ',
     + a6)
      endif
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
      do 12 ii=1,i
      if(lbnd(3,ii).le.0.and.lbnd(3,ii).gt.-7) then
      kk=-lbnd(3,i)
      do 13 jj=1,3
      if(jj.eq.1.and.(kk.eq.1.or.kk.eq.3.or.kk.eq.5)) goto 13
      if(jj.eq.2.and.(kk.eq.2.or.kk.eq.3.or.kk.eq.6)) goto 13
      if(jj.eq.3.and.(kk.eq.4.or.kk.eq.5.or.kk.eq.6)) goto 13
      if(i.eq.ii.and.j.eq.jj) goto 12
      if((catm(1+j,lbnd(1,i)).eq.catm(1+jj,lbnd(1,ii)).and.
     + catm(1+j,lbnd(2,i)).eq.catm(1+jj,lbnd(2,ii))).or.
     + (catm(1+j,lbnd(1,i)).eq.catm(1+jj,lbnd(2,ii)).and.
     + catm(1+j,lbnd(2,i)).eq.catm(1+jj,lbnd(1,ii)))) goto 3
   13 continue
      endif
   12 continue
      if(util_print('topology',print_none)) then
      if(pre_check) write(lfnout,1000)
      write(lfnout,1002)
     + catm(1+j,lbnd(1,i)),catm(1+j,lbnd(2,i)),
     + catm(1,lbnd(1,i)),catm(1,lbnd(2,i))
 1002 format('Parameters could not be found for bond type     ',
     + a6,'-',a6,' (',a6,'-',a6,')')
      endif
      pre_check=.false. 
    3 continue
      endif
    2 continue
c
      do 4 i=1,nang
      if(lang(4,i).le.0.and.lang(4,i).gt.-7) then
      k=-lang(4,i)
      do 5 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 5
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 5
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 5
      do 14 ii=1,i
      if(lang(4,ii).le.0.and.lang(4,ii).gt.-7) then
      kk=-lang(4,i)
      do 15 jj=1,3
      if(jj.eq.1.and.(kk.eq.1.or.kk.eq.3.or.kk.eq.5)) goto 15
      if(jj.eq.2.and.(kk.eq.2.or.kk.eq.3.or.kk.eq.6)) goto 15
      if(jj.eq.3.and.(kk.eq.4.or.kk.eq.5.or.kk.eq.6)) goto 15
      if(i.eq.ii.and.j.eq.jj) goto 14
      if(catm(1+j,lang(2,i)).eq.catm(1+jj,lang(2,ii)).and.
     + ((catm(1+j,lang(1,i)).eq.catm(1+jj,lang(1,ii)).and.
     + catm(1+j,lang(3,i)).eq.catm(1+jj,lang(3,ii))).or.
     + (catm(1+j,lang(1,i)).eq.catm(1+jj,lang(3,ii)).and.
     + catm(1+j,lang(3,i)).eq.catm(1+jj,lang(1,ii))))) goto 5
   15 continue
      endif
   14 continue
      if(util_print('topology',print_none)) then
      if(pre_check) write(lfnout,1000)
      write(lfnout,1003)
     + catm(1+j,lang(1,i)),catm(1+j,lang(2,i)),catm(1+j,lang(3,i)),
     + catm(1,lang(1,i)),catm(1,lang(2,i)),catm(1,lang(3,i))
 1003 format('Parameters could not be found for angle type    ',
     + a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,')')
      endif
      pre_check=.false. 
    5 continue
      endif
    4 continue
c
      do 6 i=1,ndih
      if(ldih(5,i).le.0.and.ldih(5,i).gt.-7) then
      k=-ldih(5,i)
      do 7 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 7
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 7
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 7
      do 16 ii=1,i
      if(ldih(5,ii).le.0.and.ldih(5,ii).gt.-7) then
      kk=-ldih(5,i)
      do 17 jj=1,3
      if(jj.eq.1.and.(kk.eq.1.or.kk.eq.3.or.kk.eq.5)) goto 17
      if(jj.eq.2.and.(kk.eq.2.or.kk.eq.3.or.kk.eq.6)) goto 17
      if(jj.eq.3.and.(kk.eq.4.or.kk.eq.5.or.kk.eq.6)) goto 17
      if(i.eq.ii.and.j.eq.jj) goto 16
      if((catm(1+j,ldih(1,i)).eq.catm(1+jj,ldih(1,ii)).and.
     + catm(1+j,ldih(2,i)).eq.catm(1+jj,ldih(2,ii)).and.
     + catm(1+j,ldih(3,i)).eq.catm(1+jj,ldih(3,ii)).and.
     + catm(1+j,ldih(4,i)).eq.catm(1+jj,ldih(4,ii))).or.
     + (catm(1+j,ldih(1,i)).eq.catm(1+jj,ldih(4,ii)).and.
     + catm(1+j,ldih(2,i)).eq.catm(1+jj,ldih(3,ii)).and.
     + catm(1+j,ldih(3,i)).eq.catm(1+jj,ldih(2,ii)).and.
     + catm(1+j,ldih(4,i)).eq.catm(1+jj,ldih(1,ii)))) goto 7
   17 continue
      endif
   16 continue
      if(util_print('topology',print_none)) then
      if(pre_check) write(lfnout,1000)
      write(lfnout,1004) catm(1+j,ldih(1,i)),catm(1+j,ldih(2,i)),
     + catm(1+j,ldih(3,i)),catm(1+j,ldih(4,i)),catm(1,ldih(1,i)),
     + catm(1,ldih(2,i)),catm(1,ldih(3,i)),catm(1,ldih(4,i))
 1004 format('Parameters could not be found for dihedral type ',
     + a6,'-',a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,'-',a6,')')
      endif
      pre_check=.false. 
    7 continue
      endif
    6 continue
c
      do 8 i=1,nimp
      if(limp(5,i).le.0.and.limp(5,i).gt.-7) then
      k=-limp(5,i)
      do 9 j=1,3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 9
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 9
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 9
      do 18 ii=1,i
      if(limp(5,ii).le.0.and.limp(5,ii).gt.-7) then
      kk=-limp(5,i)
      do 19 jj=1,3
      if(jj.eq.1.and.(kk.eq.1.or.kk.eq.3.or.kk.eq.5)) goto 19
      if(jj.eq.2.and.(kk.eq.2.or.kk.eq.3.or.kk.eq.6)) goto 19
      if(jj.eq.3.and.(kk.eq.4.or.kk.eq.5.or.kk.eq.6)) goto 19
      if(i.eq.ii.and.j.eq.jj) goto 18
      if(catm(1+j,limp(1,i)).eq.catm(1+jj,limp(1,ii)).and.
     + ((catm(1+j,limp(2,i)).eq.catm(1+jj,limp(2,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(3,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(4,ii))).or.
     + (catm(1+j,limp(2,i)).eq.catm(1+jj,limp(2,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(4,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(3,ii))).or.
     + (catm(1+j,limp(2,i)).eq.catm(1+jj,limp(3,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(2,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(4,ii))).or.
     + (catm(1+j,limp(2,i)).eq.catm(1+jj,limp(3,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(4,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(2,ii))).or.
     + (catm(1+j,limp(2,i)).eq.catm(1+jj,limp(4,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(2,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(3,ii))).or.
     + (catm(1+j,limp(2,i)).eq.catm(1+jj,limp(4,ii)).and.
     + catm(1+j,limp(3,i)).eq.catm(1+jj,limp(3,ii)).and.
     + catm(1+j,limp(4,i)).eq.catm(1+jj,limp(2,ii))))) goto 9
   19 continue
      endif
   18 continue
      if(util_print('topology',print_none)) then
      if(pre_check) write(lfnout,1000)
      if(ffield(1:5).eq.'amber') then
      write(lfnout,1005) catm(1+j,limp(2,i)),catm(1+j,limp(3,i)),
     + catm(1+j,limp(1,i)),catm(1+j,limp(4,i)),
     + catm(1,limp(2,i)),catm(1,limp(3,i)),
     + catm(1,limp(1,i)),catm(1,limp(4,i))
      else
      write(lfnout,1005) catm(1+j,limp(1,i)),catm(1+j,limp(2,i)),
     + catm(1+j,limp(3,i)),catm(1+j,limp(4,i)),
     + catm(1,limp(1,i)),catm(1,limp(2,i)),
     + catm(1,limp(3,i)),catm(1,limp(4,i))
      endif
 1005 format('Parameters could not be found for improper type ',
     + a6,'-',a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,'-',a6,')')
      endif
      pre_check=.false.
    9 continue 
      endif
    8 continue
c
      return
      end
