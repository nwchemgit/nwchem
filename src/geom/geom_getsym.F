      subroutine geom_getsym0(rtdb,geom,oautosym,usesymmol,
     ,     ncenter,isys,
     ,     name,tags,
     ,     threquiv,scale,coords,charge,velocities)
C     $Id: geom_getsym.F,v 1.1 2004-04-29 03:13:31 edo Exp $
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "geom.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      character*(*) name,tags(*)
      logical oautosym,usesymmol
      integer rtdb,geom
      integer isys
      integer ncenter
      double precision scale,coords(3,*),charge(*),
     ,     velocities(*),threquiv
c
      logical oprint,oprint_sym
      integer ncenter_before_auto
      integer nata              ! num. of unique atoms 
      integer k_coordsym,l_coordsym,nops
      logical geom_symmol
      external geom_symmol
      character*16 groupname
c
      oprint=ga_nodeid().eq.0
      oprint_sym=util_print('autosym symmetry information',
     ,     print_high).and.oprint
      if(oautosym) then
c     
c     ----- detect symmetry -----
c     
            if (.not. geom_cart_get(geom,ncenter,tags,coords,charge))
     $           call errquit('geom_input: geom_cart_get failed', 0,
     &       GEOM_ERR)
            if (.not. geom_vel_get(geom, velocities))
     $        call errquit('geom_input: geom_vel_get failed', 0,
     &       GEOM_ERR)
            ncenter_before_auto = ncenter
            if(usesymmol) then
               if (.not. ma_push_get(mt_dbl,3*ncenter,'coords',
     $              l_coordsym,k_coordsym))
     $              call errquit('geominput: ma',0, MA_ERR)
               usesymmol=geom_symmol(geom,tags,dbl_mb(k_coordsym),
     $              groupname,nops,threquiv)
               call dcopy(3*ncenter,dbl_mb(k_coordsym),1,coords,1)
               if(.not.ma_pop_stack(l_coordsym)) call 
     $           errquit(' o',0,GEOM_ERR)
            endif
            if(.not.usesymmol) then
               call geom_auto_sym(rtdb,geom,coords,charge,tags,ncenter,
     $              threquiv,groupname,velocities)
            endif
            if (geom_group_set(geom,groupname)) then
               if (.not.geom_cart_set(geom,ncenter,tags,coords,charge))
     $              call errquit('geom_input: geom_cart_set failed', 0,
     &       GEOM_ERR)
               if (.not. geom_vel_set(geom, velocities))
     $              call errquit('geom_input: geom_vel_set failed', 0,
     &       GEOM_ERR)
               if (ncenter_before_auto .ne. ncenter) call errquit
     $              ('geom_input: autosym bug : number of atoms wrong',
     $              ncenter, GEOM_ERR)
            else
               write(LuOut,*) ' autosym detected unknown group ', 
     &                        groupname
               call errquit('geom_input: autosym: invalid group',0,
     &       INPUT_ERR)
            endif
         endif

c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c     
         nata=ncenter
         call sym_nwc(geom,rtdb,nata,oprint_sym,scale,threquiv,nops,
     &        usesymmol)
c     
c     Check that if we used autosym that we ended up with the 
c     same no. of atoms ... if we don't then autosym and nwchemsym
c     don't agree on the orientation of point group elements
c     
         if (oautosym) then
            if (.not. geom_ncent(geom,ncenter)) call errquit
     $           ('geom_input: geom_cent?',0, GEOM_ERR)
            if (ncenter_before_auto .ne. ncenter) call errquit
     $           ('geom_input: autosym bug : too many atoms',ncenter,
     &       GEOM_ERR)
         endif
c     
         if (oprint_sym.and.(isys.eq.0)) call sym_print_all(geom,
     $        .true.,.true.,.true.,.true.,.true.)
c
         return
         end
      subroutine geom_getsym(rtdb,geom,name)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "inp.fh"
      integer rtdb,geom
c
      double precision threquiv,scale
      logical usesymmol
      character*(*) name        ! for name of geometry
      logical geom_get_user_scale,oautosym,s
      external geom_get_user_scale
      character*256 tmp
      integer k
c
      threquiv=1d-2 ! hardwired for now
      if(.not.rtdb_get(rtdb,'geom:symmol',mt_log,1,
     $     usesymmol)) usesymmol=.false.
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit('geom_input: failed getting user scale',0,
     &       GEOM_ERR)
      oautosym=.true.
c
      call geom_getsym0(rtdb,geom,oautosym,usesymmol,
     ,     ncenter(geom),isystype(geom),
     ,     name,tags(1,geom),
     ,     threquiv,scale,coords(1,1,geom),
     ,        charge(1,geom),velocities(1,1,geom))
c
c     store stuff to rtdb
c
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      s = .true.
      tmp(k:) = ' ' 
      tmp(k:) = ':no. unique centers'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group number'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group name'
      s = s .and. rtdb_cput(rtdb, tmp, 1, group_name(geom))
      return
      end
