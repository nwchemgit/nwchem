      logical function argos_prep_admin(lfncmd,filcmd,lfnout,lfnfrg)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
#include "argos_prep_common.fh"
c
      logical argos_prep_frgsiz,argos_prep_frgrd,argos_prep_master
      external argos_prep_frgsiz,argos_prep_frgrd,argos_prep_master
c
      integer lfncmd,lfnout,lfnfrg
      character*255 filcmd,card,string
      character*255 dir,file1,file2
      integer i,length
      integer natm,nparms,ipardef,nbnd,nzmat,matm,mbnd,mzmat
      integer l_latm,i_latm,l_catm,i_catm,l_qatm,i_qatm
      integer l_bnd,i_bnd,l_fil,i_fil
      integer i_izmat,l_izmat,i_zmat,l_zmat,i_ndx,l_ndx,i_ndx2,l_ndx2
c
      logical ladmin
c
      ladmin=.false.
c
c     check if the command file exists
c     --------------------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=99)
c
    1 continue
c
      read(lfncmd,3000,end=9,err=999) card
 3000 format(a)
c
      if(card(1:6).eq.'master') then
      ladmin=.true.
      string=card(8:255)
      dir='./ '
      file1=string
      length=index(string,' ')-1
      do 2 i=length,1,-1
      if(string(i:i).eq.'/') then
      dir=string(1:i)//' '
      file1=string(i+1:length)//' '
      goto 3
      endif
    2 continue
    3 continue
c
      write(filmst,'(a,a)') dir(1:index(dir,' ')-1),
     + file1(1:index(file1,' '))
      length=index(filmst,' ')-1
c
      if(.not.argos_prep_frgsiz(lfnmst,filmst,lfnout,
     + natm,nbnd,nzmat,nparms,ipardef))
     + call md_abort('Error in argos_prep_admin',0)
c
      matm=natm
      mbnd=2*matm
      if(mbnd.lt.nbnd) mbnd=nbnd
      mzmat=nzmat+1
c
c     -----------------------------
c
c     integer latm(1,matm) :
c                  2       :
c                  3       :
c                  4       :
c                  5       :
c
      if(.not.ma_push_get(mt_int,6*matm,'latm',l_latm,i_latm))
     + call md_abort('mksgm: Memory allocation failed for latm',9999)
      if(.not.ma_push_get(mt_int,matm,'ndx',l_ndx,i_ndx))
     + call md_abort('mksgm: Memory allocation failed for ndx',9999)
      if(.not.ma_push_get(mt_int,matm,'ndx',l_ndx2,i_ndx2))
     + call md_abort('mksgm: Memory allocation failed for ndx2',9999)
c
c    char*255 fil
c
      if(.not.ma_push_get(mt_byte,255*nparms,'fil',l_fil,i_fil))
     + call md_abort('mksgm: Memory allocation failed for fil',9999)
c
c     char*6  catm(1,matm) : atom name
c                  2       : atom type
c
      if(.not.ma_push_get(mt_byte,6*(nparms+1)*matm,'catm',
     + l_catm,i_catm))
     + call md_abort('mksgm: Memory allocation failed for catm',9999)
c
c     real*8 qatm(1,matm) : partial atomic charges
c                 2       : polarizabilities
c
      if(.not.ma_push_get(mt_dbl,2*nparms*matm,'qatm',l_qatm,i_qatm))
     + call md_abort('mksgm: Memory allocation failed for qatm',9999)
c
c     integer lbnd(1:2,mbnd) : bond indices
c
      if(.not.ma_push_get(mt_int,2*mbnd,'bnd',l_bnd,i_bnd))
     + call md_abort('mksgm: Memory allocation failed for bnd',9999)
c
      if(mzmat.gt.0) then
      if(.not.ma_push_get(mt_int,4*mzmat,'izmat',l_izmat,i_izmat))
     + call md_abort('mksgm: Memory allocation failed for izmat',mzmat)
      if(.not.ma_push_get(mt_dbl,3*mzmat,'zmat',l_zmat,i_zmat))
     + call md_abort('mksgm: Memory allocation failed for zmat',mzmat)
      endif
c
      if(.not.argos_prep_frgrd(lfnmst,filmst,lfnout,nparms,nparms+1,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_izmat),
     + dbl_mb(i_zmat),mzmat,nzmat,byte_mb(i_fil)))
     + call md_abort('mksgm: argos_prep_frgrd failed',9999)
c
      if(.not.argos_prep_master(lfnfrg,lfnout,nparms,nparms+1,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_izmat),
     + dbl_mb(i_zmat),mzmat,nzmat,dir,byte_mb(i_fil),
     + int_mb(i_ndx),int_mb(i_ndx2)))
     + call md_abort('mksgm: argos_prep_master failed',9999)
c
      if(mzmat.gt.0) then
      if(.not.ma_pop_stack(l_zmat))
     + call md_abort('mksgm: Memory deallocation failed for zmat',9999)
      if(.not.ma_pop_stack(l_izmat))
     + call md_abort('mksgm: Memory deallocation failed for izmat',9999)
      endif
      if(.not.ma_pop_stack(l_bnd))
     + call md_abort('mksgm: Memory deallocation failed for bnd',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call md_abort('mksgm: Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call md_abort('mksgm: Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_fil))
     + call md_abort('mksgm: Memory deallocation failed for fil',9999)
      if(.not.ma_pop_stack(l_ndx2))
     + call md_abort('mksgm: Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('mksgm: Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_latm))
     + call md_abort('mksgm: Memory deallocation failed for latm',9999)
c
      open(unit=lfnmst,file=filmst(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2001) filmst(1:length)
 2001 format(' Reading master fragment',t40,a,/)
      endif
c
    4 continue
c
      read(lfnmst,1000,end=9999,err=9999) card
 1000 format(a)
c
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 4
c
      read(card,1001) natm,nparms,ipardef
 1001 format(3i5)
c
      close(lfnmst,status='keep')
c
      endif
c
      goto 1
c
    9 continue
      close(unit=lfncmd)
c
   99 continue
  999 continue
      argos_prep_admin=ladmin
      return
 9999 continue
      call md_abort('Error in argos_prep_master',0)
      return
      end
