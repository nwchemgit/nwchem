      logical function argos_prepare_check(lfnout,nparms,mparms,ignore,
     + latt,catt,matt,natt,catm,latm,matm,natm,lbnd,mbnd,nbnd,
     + lang,mang,nang,ldih,mdih,ndih,limp,mimp,nimp,
     + lseq,cseq,mseq,nseq)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "argos_prepare_common.fh"
c
      integer lfnout,matm,natm,nparms,mparms,ignore
      integer matt,natt,mbnd,nbnd,mang,nang,mdih,ndih,mimp,nimp
      integer latt(3,matt),lbnd(4,mbnd),lang(5,mang),ldih(6,mdih)
      integer limp(6,mimp),latm(11,matm)
      character*6 catt(2,matt),catm(mparms,matm)
      integer mseq,nseq
      integer lseq(6,mseq)
      character*10 cseq(mseq)
c
      integer i,j,k,ii,jj,kk,mask
c
      argos_prepare_check=.true.
c
      mask=0
      do 123 j=1,nparms
      mask=ior(mask,2**(j-1))
  123 continue
c
c     check the atom types
c     --------------------
c
 1000 format(/,' Undetermined force field parameters',/)
c
      do 1 i=1,natt
      if(latt(1,i).le.0) then
      do 11 j=1,i-1
      if(latt(1,j).le.0.and.catt(1,i).eq.catt(1,j)) goto 1
   11 continue
      if(util_print('topology',print_none)) then
      if(argos_prepare_check) write(lfnout,1000)
      write(lfnout,1001) catt(1,i)
 1001 format(' Parameters could not be found for atom type     ',
     + a6)
      endif
      argos_prepare_check=.false.
      endif
    1 continue
c
      do 2 i=1,nbnd
      if(lbnd(3,i).le.0.and.lbnd(3,i).gt.-mask.and.lbnd(4,i).eq.0) then
      k=-lbnd(3,i)
      do 3 j=1,nparms
      if((catm(1+j,lbnd(1,i))(6:6).eq.'Q'.or.
     + catm(1+j,lbnd(1,i))(6:6).eq.'H').and.
     + (catm(1+j,lbnd(2,i))(6:6).eq.'Q'.or.
     + catm(1+j,lbnd(2,i))(6:6).eq.'H').and.ignore.ne.0) goto 3
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 3
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 3
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 3
      do 12 ii=1,i
      if(lbnd(3,ii).le.0.and.lbnd(3,ii).gt.-mask) then
      kk=-lbnd(3,i)
      do 13 jj=1,nparms
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
      if(argos_prepare_check) write(lfnout,1000)
      write(lfnout,1002)
     + catm(1+j,lbnd(1,i)),catm(1+j,lbnd(2,i)),
     + catm(1,lbnd(1,i)),catm(1,lbnd(2,i)),
     + cseq(latm(5,lbnd(1,i))),lseq(1,latm(5,lbnd(1,i))),
     + lbnd(1,i),lbnd(2,i),i,j,k
 1002 format(' Parameters could not be found for bond type     ',
     + a6,'-',a6,' (',a6,'-',a6,')',a6,i6,2i7,i7,2i3)
      endif
      argos_prepare_check=.false.
    3 continue
      endif
    2 continue
c
      do 4 i=1,nang
      if(lang(4,i).le.0.and.lang(4,i).gt.-mask.and.lang(5,i).eq.0) then
      k=-lang(4,i)
      do 5 j=1,nparms
      if((catm(1+j,lang(1,i))(6:6).eq.'Q'.or.
     + catm(1+j,lang(1,i))(6:6).eq.'H').and.
     + (catm(1+j,lang(2,i))(6:6).eq.'Q'.or.
     + catm(1+j,lang(2,i))(6:6).eq.'H').and.
     + (catm(1+j,lang(3,i))(6:6).eq.'Q'.or.
     + catm(1+j,lang(3,i))(6:6).eq.'H').and.ignore.ne.0) goto 5
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 5
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 5
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 5
      do 14 ii=1,i
      if(lang(4,ii).le.0.and.lang(4,ii).gt.-mask) then
      kk=-lang(4,i)
      do 15 jj=1,nparms
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
      if(argos_prepare_check) write(lfnout,1000)
      write(lfnout,1003)
     + catm(1+j,lang(1,i)),catm(1+j,lang(2,i)),catm(1+j,lang(3,i)),
     + catm(1,lang(1,i)),catm(1,lang(2,i)),catm(1,lang(3,i)),
     + cseq(latm(5,lang(1,i))),lseq(1,latm(5,lang(1,i))),
     + lang(1,i),lang(2,i),lang(3,i),i,j,k
 1003 format(' Parameters could not be found for angle type    ',
     + a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,')',a6,i6,3i7,i7,2i3)
      endif
      argos_prepare_check=.false.
    5 continue
      endif
    4 continue
c
      do 6 i=1,ndih
      if(ldih(5,i).le.0.and.ldih(5,i).gt.-mask.and.ldih(6,i).eq.0) then
      k=-ldih(5,i)
      do 7 j=1,nparms
      if((catm(1+j,ldih(1,i))(6:6).eq.'Q'.or.
     + catm(1+j,ldih(1,i))(6:6).eq.'H').and.
     + (catm(1+j,ldih(2,i))(6:6).eq.'Q'.or.
     + catm(1+j,ldih(2,i))(6:6).eq.'H').and.
     + (catm(1+j,ldih(3,i))(6:6).eq.'Q'.or.
     + catm(1+j,ldih(3,i))(6:6).eq.'H').and.
     + (catm(1+j,ldih(4,i))(6:6).eq.'Q'.or.
     + catm(1+j,ldih(4,i))(6:6).eq.'H').and.ignore.ne.0) goto 7
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 7
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 7
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 7
      do 16 ii=1,i
      if(ldih(5,ii).le.0.and.ldih(5,ii).gt.-mask) then
      kk=-ldih(5,i)
      do 17 jj=1,nparms
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
      if(argos_prepare_check) write(lfnout,1000)
      write(lfnout,1004) catm(1+j,ldih(1,i)),catm(1+j,ldih(2,i)),
     + catm(1+j,ldih(3,i)),catm(1+j,ldih(4,i)),catm(1,ldih(1,i)),
     + catm(1,ldih(2,i)),catm(1,ldih(3,i)),catm(1,ldih(4,i)),
     + cseq(latm(5,ldih(1,i))),lseq(1,latm(5,ldih(1,i))),
     + ldih(1,i),ldih(2,i),ldih(3,i),ldih(4,i),i,j,k
 1004 format(' Parameters could not be found for dihedral type ',
     + a6,'-',a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,'-',a6,')',
     + a6,i6,4i7,i7,2i3)
      endif
      argos_prepare_check=.false.
    7 continue
      endif
    6 continue
c
      do 8 i=1,nimp
      if(limp(5,i).le.0.and.limp(5,i).gt.-mask.and.limp(6,i).eq.0) then
      k=-limp(5,i)
      do 9 j=1,nparms
      if((catm(1+j,limp(1,i))(6:6).eq.'Q'.or.
     + catm(1+j,limp(1,i))(6:6).eq.'H').and.
     + (catm(1+j,limp(2,i))(6:6).eq.'Q'.or.
     + catm(1+j,limp(2,i))(6:6).eq.'H').and.
     + (catm(1+j,limp(3,i))(6:6).eq.'Q'.or.
     + catm(1+j,limp(3,i))(6:6).eq.'H').and.
     + (catm(1+j,limp(4,i))(6:6).eq.'Q'.or.
     + catm(1+j,limp(4,i))(6:6).eq.'H').and.ignore.ne.0) goto 9
      if(j.eq.1.and.(k.eq.1.or.k.eq.3.or.k.eq.5)) goto 9
      if(j.eq.2.and.(k.eq.2.or.k.eq.3.or.k.eq.6)) goto 9
      if(j.eq.3.and.(k.eq.4.or.k.eq.5.or.k.eq.6)) goto 9
      do 18 ii=1,i
      if(limp(5,ii).le.0.and.limp(5,ii).gt.-mask) then
      kk=-limp(5,i)
      do 19 jj=1,nparms
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
      if(argos_prepare_check) write(lfnout,1000)
      if(ffield(1:5).eq.'amber') then
      write(lfnout,1005) catm(1+j,limp(2,i)),catm(1+j,limp(3,i)),
     + catm(1+j,limp(1,i)),catm(1+j,limp(4,i)),
     + catm(1,limp(2,i)),catm(1,limp(3,i)),
     + catm(1,limp(1,i)),catm(1,limp(4,i)),
     + cseq(latm(5,limp(1,i))),lseq(1,latm(5,limp(1,i))),
     + limp(2,i),limp(3,i),limp(1,i),limp(4,i),i,j,k
      else
      write(lfnout,1005) catm(1+j,limp(1,i)),catm(1+j,limp(2,i)),
     + catm(1+j,limp(3,i)),catm(1+j,limp(4,i)),
     + catm(1,limp(1,i)),catm(1,limp(2,i)),
     + catm(1,limp(3,i)),catm(1,limp(4,i)),
     + cseq(latm(5,limp(1,i))),lseq(1,latm(5,limp(1,i))),
     + limp(1,i),limp(2,i),limp(3,i),limp(4,i),i,j,k
      endif
 1005 format(' Parameters could not be found for improper type ',
     + a6,'-',a6,'-',a6,'-',a6,' (',a6,'-',a6,'-',a6,'-',a6,')',
     + a6,i6,4i7,i7,2i3)
      endif
      argos_prepare_check=.false.
    9 continue
      endif
    8 continue
c
      return
      end
 
