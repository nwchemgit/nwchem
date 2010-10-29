      logical function argos_prepare_mksgm(lfnout,lfnfrg,filfrg,lfnsgm,
     + filsgm)
c
c $Id$
c
c     function to make segment file
c
c     in  : lfnfrg    = fragment file logical file number
c           lfnsgm    = segment file logical file number
c           lfnout    = output file logical file number
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
#include "argos_prepare_common.fh"
c
      logical argos_prepare_frgsiz,argos_prepare_frgrd,
     + argos_prepare_ang,argos_prepare_dih,argos_prepare_imp
      logical argos_prepare_dihsel,argos_prepare_dihimp,
     + argos_prepare_wrtsgm
      external argos_prepare_frgsiz,argos_prepare_frgrd,
     + argos_prepare_ang,argos_prepare_dih,argos_prepare_imp
      external argos_prepare_dihsel,argos_prepare_dihimp,
     + argos_prepare_wrtsgm
c
      integer lfnout,lfnfrg,lfnsgm
      character*255 filfrg,filsgm
c
      integer length
      integer matm,mbnd,mang,mdih,mimp,mzmat
      integer natm,nbnd,nang,ndih,nimp,nzmat
c
      integer l_latm,i_latm,l_catm,i_catm,l_qatm,i_qatm
      integer l_bnd,i_bnd,l_ang,i_ang,l_dih,i_dih,l_imp,i_imp
      integer i_izmat,l_izmat,i_zmat,l_zmat,l_fil,i_fil
c
      integer nparms,ipardef
c
      length=index(filsgm,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filsgm(1:length)
 2000 format(//,'SEGMENT FILE ',a)
      endif
c
      if(.not.argos_prepare_frgsiz(lfnfrg,filfrg,lfnout,natm,nbnd,nzmat,
     + nparms,ipardef))
     + call md_abort('argos_prepare_frgsiz failed',9999)
c
c     allocate memory for bonds, angles, torsions and impropers
c     ---------------------------------------------------------
c
      matm=natm
      mbnd=4*matm
      if(mbnd.lt.nbnd) mbnd=nbnd
      mang=6*matm
      mdih=12*matm
      mimp=matm
      nang=0
      ndih=0
      nimp=0
      mzmat=nzmat+1
c
c     allocate memory for atom list
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
c     integer lang(1:3,mang) : angle indices
c
      if(.not.ma_push_get(mt_int,3*mang,'ang',l_ang,i_ang))
     + call md_abort('mksgm: Memory allocation failed for ang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c
      if(.not.ma_push_get(mt_int,4*mdih,'dih',l_dih,i_dih))
     + call md_abort('mksgm: Memory allocation failed for dih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c
      if(.not.ma_push_get(mt_int,4*mimp,'imp',l_imp,i_imp))
     + call md_abort('mksgm: Memory allocation failed for imp',9999)
c
      if(mzmat.gt.0) then
      if(.not.ma_push_get(mt_int,4*mzmat,'izmat',l_izmat,i_izmat))
     + call md_abort('mksgm: Memory allocation failed for izmat',mzmat)
      if(.not.ma_push_get(mt_dbl,3*mzmat,'zmat',l_zmat,i_zmat))
     + call md_abort('mksgm: Memory allocation failed for zmat',mzmat)
      endif
c
c     read atom list and bond list from fragment file
c     -----------------------------------------------
c
      if(.not.argos_prepare_frgrd(lfnfrg,filfrg,lfnout,nparms,nparms+1,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_izmat),
     + dbl_mb(i_zmat),mzmat,nzmat,byte_mb(i_fil)))
     + call md_abort('mksgm: argos_prepare_frgrd failed',9999)
c
c     make the angle list
c     -------------------
c
      if(.not.argos_prepare_ang(int_mb(i_bnd),mbnd,nbnd,int_mb(i_ang),
     + mang,nang))
     + call md_abort('mksgm: argos_prepare_ang failed',9999)
c
c     make the dihedral list
c     ----------------------
c
      if(.not.argos_prepare_dih(int_mb(i_ang),mang,nang,int_mb(i_dih),
     + mdih,ndih))
     + call md_abort('mksgm: argos_prepare_dih failed',9999)
c
c     make the improper list
c     ----------------------
c
      if(.not.argos_prepare_imp(int_mb(i_latm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_imp),mimp,nimp))
     + call md_abort('mksgm: argos_prepare_imp failed',9999)
c
c     select dihedrals
c     ----------------
c
      if(.not.argos_prepare_dihsel(byte_mb(i_catm),nparms,nparms+1,matm,
     + natm,
     + int_mb(i_dih),mdih,ndih))
     + call md_abort('mksgm: argos_prepare_dihsel failed',9999)
c
c     switch propers to impropers
c     ---------------------------
c
      if(.not.argos_prepare_dihimp(int_mb(i_latm),matm,natm,
     + int_mb(i_ang),mang,nang,
     + int_mb(i_dih),mdih,ndih,int_mb(i_imp),mimp,nimp))
     + call md_abort('mksgm: argos_prepare_dihsel failed',9999)
c
c     write segment file
c     ------------------
c
      if(.not.argos_prepare_wrtsgm(lfnout,lfnsgm,filsgm,nparms,nparms+1,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_ang),mang,nang,
     + int_mb(i_dih),mdih,ndih,int_mb(i_imp),mimp,nimp,
     + int_mb(i_izmat),dbl_mb(i_zmat),mzmat,nzmat,ipardef))
     + call md_abort('mksgm: argos_prepare_wrtsgm failed',9999)
c
c     deallocate memory
c     -----------------
c
      if(mzmat.gt.0) then
      if(.not.ma_pop_stack(l_zmat))
     + call md_abort('mksgm: Memory deallocation failed for zmat',9999)
      if(.not.ma_pop_stack(l_izmat))
     + call md_abort('mksgm: Memory deallocation failed for izmat',9999)
      endif
      if(.not.ma_pop_stack(l_imp))
     + call md_abort('mksgm: Memory deallocation failed for imp',9999)
      if(.not.ma_pop_stack(l_dih))
     + call md_abort('mksgm: Memory deallocation failed for dih',9999)
      if(.not.ma_pop_stack(l_ang))
     + call md_abort('mksgm: Memory deallocation failed for ang',9999)
      if(.not.ma_pop_stack(l_bnd))
     + call md_abort('mksgm: Memory deallocation failed for bnd',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call md_abort('mksgm: Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call md_abort('mksgm: Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_fil))
     + call md_abort('mksgm: Memory deallocation failed for fil',9999)
      if(.not.ma_pop_stack(l_latm))
     + call md_abort('mksgm: Memory deallocation failed for latm',9999)
c
      argos_prepare_mksgm=.true.
      return
c
 9999 continue
      argos_prepare_mksgm=.false.
      return
      end
 
