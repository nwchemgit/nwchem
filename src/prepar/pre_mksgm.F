      logical function pre_mksgm(lfnout,ffield,
     + lfnfrg,filfrg,lfnsgm,filsgm)
c
c $Id: pre_mksgm.F,v 1.10 2000-07-13 23:13:02 edo Exp $
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
c
      logical pre_frgsiz,pre_frgrd,pre_ang,pre_dih,pre_imp
      logical pre_dihsel,pre_dihimp,pre_wrtsgm
      external pre_frgsiz,pre_frgrd,pre_ang,pre_dih,pre_imp
      external pre_dihsel,pre_dihimp,pre_wrtsgm
c
      integer lfnout,lfnfrg,lfnsgm
      character*80 ffield
      character*255 filfrg,filsgm
c
      integer length
      integer matm,mbnd,mang,mdih,mimp
      integer natm,nbnd,nang,ndih,nimp
c
      integer l_latm,i_latm,l_catm,i_catm,l_qatm,i_qatm
      integer l_bnd,i_bnd,l_ang,i_ang,l_dih,i_dih,l_imp,i_imp
c
      length=index(filsgm,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filsgm(1:length)
 2000 format(//,'SEGMENT FILE ',a)
      endif
c
      if(.not.pre_frgsiz(lfnfrg,filfrg,lfnout,natm,nbnd))
     + call errquit('pre_frgsiz failed',9999)
c
c     allocate memory for bonds, angles, torsions and impropers
c     ---------------------------------------------------------
c
      matm=natm
      mbnd=2*matm
      if(mbnd.lt.nbnd) mbnd=nbnd
      mang=2*matm
      mdih=3*matm
      mimp=matm
      nang=0
      ndih=0
      nimp=0
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
      if(.not.ma_push_get(mt_int,5*matm,'latm',l_latm,i_latm))
     + call errquit('mksgm: Memory allocation failed for latm',9999)
c
c     char*6  catm(1,matm) : atom name
c                  2       : atom type
c
      if(.not.ma_push_get(mt_byte,12*matm,'catm',l_catm,i_catm))
     + call errquit('mksgm: Memory allocation failed for catm',9999)
c
c     real*8 qatm(1,matm) : partial atomic charges
c                 2       : polarizabilities
c
      if(.not.ma_push_get(mt_dbl,2*matm,'qatm',l_qatm,i_qatm))
     + call errquit('mksgm: Memory allocation failed for qatm',9999)
c
c     integer lbnd(1:2,mbnd) : bond indices
c
      if(.not.ma_push_get(mt_int,2*mbnd,'bnd',l_bnd,i_bnd))
     + call errquit('mksgm: Memory allocation failed for bnd',9999)
c
c     integer lang(1:3,mang) : angle indices
c
      if(.not.ma_push_get(mt_int,3*mang,'ang',l_ang,i_ang))
     + call errquit('mksgm: Memory allocation failed for ang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c
      if(.not.ma_push_get(mt_int,4*mdih,'dih',l_dih,i_dih))
     + call errquit('mksgm: Memory allocation failed for dih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c
      if(.not.ma_push_get(mt_int,4*mimp,'imp',l_imp,i_imp))
     + call errquit('mksgm: Memory allocation failed for imp',9999)
c
c     read atom list and bond list from fragment file
c     -----------------------------------------------
c
      if(.not.pre_frgrd(lfnfrg,filfrg,lfnout,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd))
     + call errquit('mksgm: pre_frgrd failed',9999)
c
c     make the angle list
c     -------------------
c
      if(.not.pre_ang(int_mb(i_bnd),mbnd,nbnd,int_mb(i_ang),mang,nang))
     + call errquit('mksgm: pre_ang failed',9999)
c
c     make the dihedral list
c     ----------------------
c
      if(.not.pre_dih(int_mb(i_ang),mang,nang,int_mb(i_dih),mdih,ndih))
     + call errquit('mksgm: pre_dih failed',9999)
c
c     make the improper list
c     ----------------------
c
      if(.not.pre_imp(ffield,int_mb(i_latm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_imp),mimp,nimp))
     + call errquit('mksgm: pre_imp failed',9999)
c
c     select dihedrals
c     ----------------
c
      if(.not.pre_dihsel(ffield,byte_mb(i_catm),matm,natm,
     + int_mb(i_dih),mdih,ndih))
     + call errquit('mksgm: pre_dihsel failed',9999)
c
c     switch propers to impropers
c     ---------------------------
c
      if(.not.pre_dihimp(ffield,int_mb(i_latm),matm,natm,
     + int_mb(i_ang),mang,nang,
     + int_mb(i_dih),mdih,ndih,int_mb(i_imp),mimp,nimp))
     + call errquit('mksgm: pre_dihsel failed',9999)
c
c     write segment file
c     ------------------
c
      if(.not.pre_wrtsgm(lfnout,ffield,lfnsgm,filsgm,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_bnd),mbnd,nbnd,int_mb(i_ang),mang,nang,
     + int_mb(i_dih),mdih,ndih,int_mb(i_imp),mimp,nimp))
     + call errquit('mksgm: pre_wrtsgm failed',9999)
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_imp))
     + call errquit('mksgm: Memory deallocation failed for imp',9999)
      if(.not.ma_pop_stack(l_dih))
     + call errquit('mksgm: Memory deallocation failed for dih',9999)
      if(.not.ma_pop_stack(l_ang))
     + call errquit('mksgm: Memory deallocation failed for ang',9999)
      if(.not.ma_pop_stack(l_bnd))
     + call errquit('mksgm: Memory deallocation failed for bnd',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call errquit('mksgm: Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call errquit('mksgm: Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_latm))
     + call errquit('mksgm: Memory deallocation failed for latm',9999)
c
      pre_mksgm=.true.
      return
c
 9999 continue
      pre_mksgm=.false.
      return
      end
