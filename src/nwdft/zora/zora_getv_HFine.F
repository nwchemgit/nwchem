      subroutine zora_getv_HFine_slow(
     &                           rtdb, 
     &                           g_dens,        ! in  : atomic density 
     &                           ofinite,       ! in  : = .true. if Gaussian Nucl. Model requested
     &                           atmass,        ! in  : atomic mass
     &                           xyz_NMRcoords, ! in  : nuclear coordinates
     &                           g_zpso,        ! out : ZPSO term
     &                           g_fcsd,        ! out : FC+SD (v,u) term
     &                           nexc)
c
C$Id: zora_getv_EPR.F,v 1.3 2010-11-01 FA Exp $
c Adapted from zora_getv_so     

      implicit none
#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "cgridfile.fh"
#include "grid_cube.fh"
c
c     == arguments ==
      integer rtdb
      integer g_dens(2)
      integer g_zpso(3),g_fcsd(3,3)
      integer nexc
c
c     == local variables ==
      integer i,j,k,t,n,ind,nij
      double precision rho_n
      double precision tmat
      double precision dummy(2)
      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwght,iqwght,nqpts,ncube,
     &     ictr_buf,iqpts
      double precision rad,ke
      integer lbas_cset_info, ibas_cset_info,
     &        lbas_cent_info, ibas_cent_info,
     &        ldocset, idocset,
     &        l_rchi_atom,i_rchi_atom,
     &        l_rq,i_rq,lniz, iniz,
     &        lchi_ao, ichi_ao,
     &        ldelchi_ao, idelchi_ao
      integer lzpso(3),izpso(3),
     &        lfcsd(3,3),ifcsd(3,3)
      integer inntsize,ddblsize,ok
      double precision xyz_NMRcoords(3),atmass
      double precision chi_cntr(3,nbf_ao)

      logical grid_file_rewind,ofinite
      external grid_file_rewind,calc_zora_HFine_slow,
     &         ga_antisymmetrize
c
c     == allocate memory ==
      do i=1,3
       if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lzpso(i),izpso(i)))
     &    call errquit('zora_getv_HFine: zpso',911,MA_ERR)
       do j=1,3
        if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lfcsd',lfcsd(i,j),ifcsd(i,j)))
     &    call errquit('zora_getv_HFine: fcsd',911,MA_ERR)
       enddo
      enddo ! end-loop-i
c     == preliminaries ==
       do j=1,3
         do i= 1, nbf_ao*nbf_ao
          dbl_mb(izpso(j)+i-1)=0.d0
         enddo  
       enddo ! end-loop-j
       do k=1,3
        do j=1,3
         do i= 1, nbf_ao*nbf_ao
         dbl_mb(ifcsd(k,j)+i-1)=0.d0
         enddo
        enddo
       enddo
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      call grid_quadv0(rtdb,g_dens,g_zpso,nexc,rho_n,dummy,tmat)
c     == ao basis set info used by xc_eval_basis ==
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('zora_getv_sf:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('zora_getv_sf: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("zora_getv_sf:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('zora_getv_sf: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('zora_getv_sf: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('zora_getv_sf: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('zora_getv_sf: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('zora_getv_sf: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('zora_getv_sf: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('zora_getv_sf: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('zora_getv_sf: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
c     == loop over records in the grid file ==
      do iqsh = 1, n_rec_in_file
c
c       == define the current range of radial shells and integration center ==
        call grid_file_read(n_per_rec, nqpts, ictr_buf,
     &        rad,dbl_mb(irqbuf),nsubb)

        if(nqpts.gt.buffer_size) 
     &    call errquit(' buffersize exceed by qpts ',nqpts, UNKNOWN_ERR)
c
c        == loop over a subset of the grid ==
         istep=0
         do  ncube=1,nsubb
c
c           put buf into currently used arrays qxyz and qwght
            call grid_repack(dbl_mb(irqbuf), dbl_mb(iqxyz),
     &           dbl_mb(iqwght), nqpts, rad,istep)

            if(nqpts.ne.0) then
c
c              == compute the basis functions over the grid ==
               if(.not.MA_Push_get(MT_dbl, ncenters, 'rchi_atom',
     &             l_rchi_atom,i_rchi_atom))
     &             call errquit("zora_getv:rchi_atom",0, MA_ERR)
c
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &             l_rq,i_rq))
     &             call errquit("zora_getv_sf:rq",0, MA_ERR)
c
c              == delchi ==
               if (.not.MA_Push_Get(mt_dbl, 3*nqpts*nbf_ao,
     &             'delchi_ao', ldelchi_ao, idelchi_ao))
     &             call errquit('zora_getv: delchi_ao',0, MA_ERR)
c
c              == chi ==
               if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &             'chi_ao', lchi_ao, ichi_ao))
     &             call errquit('zora_getv: chi_ao',0, MA_ERR)
               call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
               call xc_eval_basis(ao_bas_han, 1, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
                call calc_zora_HFine_slow(
     &                                  ao_bas_han,geom,ipol,g_dens,     
     &                                  dbl_mb(ichi_ao),  
     &                                  dbl_mb(idelchi_ao),
     &                                  dbl_mb(iqxyz),dbl_mb(iqwght),        
     &                                  nbf_ao,nqpts,ncenters,   
     &                                  ofinite,
     &                                  atmass,
     &                                  xyz_NMRcoords,
     &                                  dbl_mb(izpso(1)),  ! out
     &                                  dbl_mb(izpso(2)),  ! out
     &                                  dbl_mb(izpso(3)),  ! out
     &                                  dbl_mb(ifcsd(1,1)),! out
     &                                  dbl_mb(ifcsd(1,2)),! out
     &                                  dbl_mb(ifcsd(1,3)),! out
     &                                  dbl_mb(ifcsd(2,1)),! out
     &                                  dbl_mb(ifcsd(2,2)),! out
     &                                  dbl_mb(ifcsd(2,3)),! out
     &                                  dbl_mb(ifcsd(3,1)),! out
     &                                  dbl_mb(ifcsd(3,2)),! out
     &                                  dbl_mb(ifcsd(3,3)))! out

c              == delete memory ==
               if(.not.MA_pop_stack(lchi_ao))
     &            call errquit("zora_getv: pop chi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(ldelchi_ao))
     &            call errquit("zora_getv: pop delchi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rq))
     &            call errquit("zora_getv: pop rq", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rchi_atom))
     &            call errquit("zora_getv: pop rchi_atom",100,MA_ERR)
            endif ! nqpts
         enddo ! ncube
      end do ! iqsh
c
c     == delete memory ==
      if(.not.MA_pop_stack(lrqbuf))
     &     call errquit("zora_getv_sf: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("zora_getv_sf: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("zora_getv_sf: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("zora_getv_sf: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("zora_getv_sf: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("zora_getv_sf: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("zora_getv_sf: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("zora_getv_sf: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("zora_getv_sf: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("zora_getv_sf: pop bas_cent_info", 100, MA_ERR)
c
c     == tally up over all the nodes ==
      nij = 0
      do i = 1, nbf_ao
        do j = 1, nbf_ao
         nij = nij + 1
         do k=1,3
          call ga_dgop(msg_excrho, dbl_mb(izpso(k)+nij-1), 1, '+')
          do t=1,3
           call ga_dgop(msg_excrho,dbl_mb(ifcsd(k,t)+nij-1), 1, '+')
          enddo ! end-loop-t
         enddo ! end-loop-k
        enddo ! end-loop-j
      enddo ! end-loop-i
c
c     == pack into a ga all HFine AOs ==
      do i=1,3 
       call ga_zero(g_zpso(i)) 
       call ga_put(g_zpso(i),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(izpso(i)),nbf_ao)
       do j=1,3
        call ga_zero(g_fcsd(i,j)) 
        call ga_put(g_fcsd(i,j),
     &              1,nbf_ao,1,nbf_ao,dbl_mb(ifcsd(i,j)),nbf_ao)
       enddo ! end-loop-j
      enddo ! end-loop-i
      call ga_sync()
c ---- Free memory
      do i=1,3
        if (.not.ma_free_heap(lzpso(i))) call errquit
     &   ('zora_getv_HFine: ma_free_heap of lzpso failed',
     &    911,MA_ERR)
       do j=1,3
        if (.not.ma_free_heap(lfcsd(i,j))) call errquit
     &   ('zora_getv_HFine: ma_free_heap of lfcsd failed',
     &    911,MA_ERR)
       enddo
      enddo
      return
      end
c ++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++
      subroutine zora_getv_HFine_fast( 
     &                           rtdb, 
     &                           g_dens,        ! in  : atomic density 
     &                           ofinite,       ! in  : = .true. if Gaussian Nucl. Model requested
     &                           zetanuc_arr,   ! in  : sqrt(zetanuc(i)) i=1,natoms
     &                           zetanuc_slc,   ! in  : zetanuc(i)
     &                           Knucl,
     &                           xyz_NMRcoords, ! in  : nuclear coordinates
     &                           g_zpso,        ! out : ZPSO term
     &                           g_fcsd,        ! out : FC+SD (v,u) term
     &                           nexc)
c
C$Id: zora_getv_EPR.F,v 1.3 2010-11-01 FA Exp $
c Adapted from zora_getv_so     

      implicit none
#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "cgridfile.fh"
#include "grid_cube.fh"
c
c     == arguments ==
      integer rtdb
      integer g_dens(2)
      integer g_zpso(6),g_fcsd(3,3)
      integer nexc
c
c     == local variables ==
      integer i,j,k,t,m,n,ind,nij,ac
      double precision rho_n
      double precision tmat
      double precision dummy(2)
      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwght,iqwght,nqpts,ncube,
     &     ictr_buf,iqpts
      double precision rad,ke
      double precision zetanuc_arr(*) ! For Gaussian Nuclear Model
      double precision zetanuc_slc    ! For Gaussian Nuclear Model
      integer lbas_cset_info, ibas_cset_info,
     &        lbas_cent_info, ibas_cent_info,
     &        ldocset, idocset,
     &        l_rchi_atom,i_rchi_atom,
     &        l_rq,i_rq,lniz, iniz,
     &        lchi_ao, ichi_ao,
     &        ldelchi_ao, idelchi_ao
      integer lzpso(3,3),izpso(3,3),
     &        lfcsd(3,3),ifcsd(3,3)
      integer inntsize,ddblsize,ok,atn
      double precision xyz_NMRcoords(3),atmass   
      double precision chi_cntr(3,nbf_ao)
      logical grid_file_rewind,ofinite,Knucl
      integer ind_tmn(2,3) 
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312 
      external grid_file_rewind,calc_zora_HFine_fast,
     &         ga_antisymmetrize
c
c     == allocate memory ==
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)
        if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                'lnmr',lzpso(m,n),izpso(m,n)))
     &    call errquit('zora_getv_NMR: zpso',911,MA_ERR)
        if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                'lnmr',lzpso(n,m),izpso(n,m)))
     &    call errquit('zora_getv_NMR: zpso',911,MA_ERR)
       enddo ! end-loop-t
       do i=1,3
        do j=1,3
          if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lfcsd',lfcsd(i,j),ifcsd(i,j)))
     &    call errquit('zora_getv_HFine: fcsd',911,MA_ERR)
        enddo
       enddo ! end-loop-i
c     == preliminaries ==
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)     
         do i= 1, nbf_ao*nbf_ao
          dbl_mb(izpso(m,n)+i-1)=0.d0
          dbl_mb(izpso(n,m)+i-1)=0.d0
         enddo ! end-loop-i   
       enddo ! end-loop-t
       do k=1,3
        do j=1,3
         do i= 1, nbf_ao*nbf_ao
         dbl_mb(ifcsd(k,j)+i-1)=0.d0
         enddo
        enddo
       enddo
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      call grid_quadv0(rtdb,g_dens,g_zpso,nexc,rho_n,dummy,tmat)
c     == ao basis set info used by xc_eval_basis ==
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('zora_getv_sf:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('zora_getv_sf: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("zora_getv_sf:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('zora_getv_sf: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('zora_getv_sf: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('zora_getv_sf: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('zora_getv_sf: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('zora_getv_sf: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('zora_getv_sf: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('zora_getv_sf: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('zora_getv_sf: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
c     == loop over records in the grid file ==
      do iqsh = 1, n_rec_in_file
c
c       == define the current range of radial shells and integration center ==
        call grid_file_read(n_per_rec, nqpts, ictr_buf,
     &        rad,dbl_mb(irqbuf),nsubb)

        if(nqpts.gt.buffer_size) 
     &    call errquit(' buffersize exceed by qpts ',nqpts, UNKNOWN_ERR)
c
c        == loop over a subset of the grid ==
         istep=0
         do  ncube=1,nsubb
c
c           put buf into currently used arrays qxyz and qwght
            call grid_repack(dbl_mb(irqbuf), dbl_mb(iqxyz),
     &           dbl_mb(iqwght), nqpts, rad,istep)

            if(nqpts.ne.0) then
c
c              == compute the basis functions over the grid ==
               if(.not.MA_Push_get(MT_dbl, ncenters, 'rchi_atom',
     &             l_rchi_atom,i_rchi_atom))
     &             call errquit("zora_getv:rchi_atom",0, MA_ERR)
c
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &             l_rq,i_rq))
     &             call errquit("zora_getv_sf:rq",0, MA_ERR)
c
c              == delchi ==
               if (.not.MA_Push_Get(mt_dbl, 3*nqpts*nbf_ao,
     &             'delchi_ao', ldelchi_ao, idelchi_ao))
     &             call errquit('zora_getv: delchi_ao',0, MA_ERR)
c
c              == chi ==
               if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &             'chi_ao', lchi_ao, ichi_ao))
     &             call errquit('zora_getv: chi_ao',0, MA_ERR)
               call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
               call xc_eval_basis(ao_bas_han, 1, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
                call calc_zora_HFine_fast(
     &                                ao_bas_han,geom,ipol,g_dens,     
     &                                dbl_mb(ichi_ao),  
     &                                dbl_mb(idelchi_ao),
     &                                dbl_mb(iqxyz),dbl_mb(iqwght),        
     &                                nbf_ao,nqpts,ncenters,   
     &                                ofinite,
     &                                zetanuc_arr,
     &                                zetanuc_slc,
     &                                Knucl,
     &                                xyz_NMRcoords,
     &                                dbl_mb(izpso(1,2)),! out
     &                                dbl_mb(izpso(1,3)),! out
     &                                dbl_mb(izpso(2,1)),! out
     &                                dbl_mb(izpso(2,3)),! out
     &                                dbl_mb(izpso(3,1)),! out
     &                                dbl_mb(izpso(3,2)),! out
     &                                dbl_mb(ifcsd(1,1)),! out
     &                                dbl_mb(ifcsd(1,2)),! out
     &                                dbl_mb(ifcsd(1,3)),! out
     &                                dbl_mb(ifcsd(2,1)),! out
     &                                dbl_mb(ifcsd(2,2)),! out
     &                                dbl_mb(ifcsd(2,3)),! out
     &                                dbl_mb(ifcsd(3,1)),! out
     &                                dbl_mb(ifcsd(3,2)),! out
     &                                dbl_mb(ifcsd(3,3)))! out
c              == delete memory ==
               if(.not.MA_pop_stack(lchi_ao))
     &            call errquit("zora_getv: pop chi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(ldelchi_ao))
     &            call errquit("zora_getv: pop delchi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rq))
     &            call errquit("zora_getv: pop rq", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rchi_atom))
     &            call errquit("zora_getv: pop rchi_atom",100,MA_ERR)
            endif ! nqpts
         enddo ! ncube
      end do ! iqsh
c
c     == delete memory ==
      if(.not.MA_pop_stack(lrqbuf))
     &     call errquit("zora_getv_sf: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("zora_getv_sf: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("zora_getv_sf: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("zora_getv_sf: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("zora_getv_sf: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("zora_getv_sf: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("zora_getv_sf: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("zora_getv_sf: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("zora_getv_sf: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("zora_getv_sf: pop bas_cent_info", 100, MA_ERR)
c
c     == tally up over all the nodes ==
      nij = 0
      do i = 1, nbf_ao
        do j = 1, nbf_ao
         nij = nij + 1
         do k=1,3
          m=ind_tmn(1,k)
          n=ind_tmn(2,k) 
          call ga_dgop(msg_excrho,dbl_mb(izpso(m,n)+nij-1), 1, '+')
          call ga_dgop(msg_excrho,dbl_mb(izpso(n,m)+nij-1), 1, '+')
          do t=1,3
           call ga_dgop(msg_excrho,dbl_mb(ifcsd(k,t)+nij-1), 1, '+')
          enddo ! end-loop-t
         enddo ! end-loop-k
        enddo ! end-loop-j
      enddo ! end-loop-i
c
c     == pack into a ga all HFine AOs ==
      ac=1
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
       call ga_zero(g_zpso(ac)) 
       call ga_put(g_zpso(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(izpso(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      do t=1,3
       m=ind_tmn(2,t)
       n=ind_tmn(1,t) 
       call ga_zero(g_zpso(ac)) 
       call ga_put(g_zpso(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(izpso(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      call ga_sync()
      do i=1,3 
       do j=1,3
        call ga_zero(g_fcsd(i,j)) 
        call ga_put(g_fcsd(i,j),
     &              1,nbf_ao,1,nbf_ao,dbl_mb(ifcsd(i,j)),nbf_ao)
       enddo ! end-loop-j
      enddo ! end-loop-i
      call ga_sync()
c ----- free memory ----------- START
      do i=1,3
       do j=1,3
        if (.not.ma_free_heap(lfcsd(i,j))) call errquit
     &   ('zora_getv_HFine: ma_free_heap of lfcsd failed',
     &    911,MA_ERR)
       enddo
      enddo
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
      if (.not.ma_free_heap(lzpso(m,n))) call errquit
     &   ('zora_getv_HFine: ma_free_heap of lzpso failed',
     &    911,MA_ERR)
      if (.not.ma_free_heap(lzpso(n,m))) call errquit
     &   ('zora_getv_HFine: ma_free_heap of lzpso failed',
     &    911,MA_ERR)
      enddo ! end-loop-t
c ----- free memory ----------- END
      return
      end

      subroutine zora_getv_NMRHFine_F1ji(
     &          rtdb, 
     &          g_dens,       ! in : atomic density 
     &          g_hfine,      ! out: 
     &          natoms,       ! in: nr. atoms 
     &          ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &          zetanuc_arr,  ! in: sqrt(zetanuc(i)) i=1,natoms for Gaussian Nuclear Model
     &          Knucl,        ! in: = .true. if K_ZORA(V=Nuclear pot. only)
     &          nexc)
c
C$Id: zora_getv_NMRCS_SR.F,v 1.3 2010-11-04 FA Exp $
c Adapted from zora_getv_sf

      implicit none
#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "cgridfile.fh"
#include "grid_cube.fh"
c
c     == arguments ==
      integer rtdb
      integer g_dens(2)
      integer g_hfine(3)
      integer nexc
c     == local variables ==
      integer i,j,k,t,a,m,n,ind,nij,ac
      double precision rho_n
      double precision tmat
      double precision dummy(2)
      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwght,iqwght,nqpts,ncube,
     &     ictr_buf,iqpts
      double precision rad,ke
      integer lbas_cset_info, ibas_cset_info,
     &        lbas_cent_info, ibas_cent_info,
     &        ldocset, idocset,
     &        l_rchi_atom,i_rchi_atom,
     &        l_rq,i_rq,lniz, iniz,
     &        lchi_ao, ichi_ao,
     &        ldelchi_ao, idelchi_ao
      double precision xyz_NMRcoords(3)   
      integer lhfine(3)   ,ihfine(3)
      integer inntsize,ddblsize,ok
      logical grid_file_rewind
      integer natoms
      logical ofinite,Knucl
      double precision zetanuc_arr(natoms)
      external grid_file_rewind,calc_NMRHFine_F1ij,
     &         ga_antisymmetrize,grid_quadv0
c
c     == allocate memory ==
      do i=1,3
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lhfine',lhfine(i),ihfine(i)))
     &    call errquit('zora_getv_NMR: hfine',911,MA_ERR)
      enddo ! end-loop-i
c     == preliminaries ==
       do j=1,3
         do i= 1, nbf_ao*nbf_ao
          dbl_mb(ihfine(j)+i-1) =0.d0
         enddo  
       enddo
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      call grid_quadv0(rtdb,g_dens,g_hfine(1),nexc,rho_n,dummy,tmat)
c     == ao basis set info used by xc_eval_basis ==
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('zora_getv_sf:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('zora_getv_sf: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('zora_getv_sf: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("zora_getv_sf:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('zora_getv_sf: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('zora_getv_sf: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('zora_getv_sf: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('zora_getv_sf: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('zora_getv_sf: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('zora_getv_sf: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('zora_getv_sf: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('zora_getv_sf: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
c     == loop over records in the grid file ==
      do iqsh = 1, n_rec_in_file
c       == define the current range of radial shells and integration center ==
        call grid_file_read(n_per_rec, nqpts, ictr_buf,
     &        rad,dbl_mb(irqbuf),nsubb)
        if(nqpts.gt.buffer_size) 
     &    call errquit(' buffersize exceed by qpts ',nqpts, UNKNOWN_ERR)
c        == loop over a subset of the grid ==
         istep=0
         do  ncube=1,nsubb
c        put buf into currently used arrays qxyz and qwght
            call grid_repack(dbl_mb(irqbuf), dbl_mb(iqxyz),
     &           dbl_mb(iqwght), nqpts, rad,istep)         
            if(nqpts.ne.0) then
c        == compute the basis functions over the grid ==
               if(.not.MA_Push_get(MT_dbl, ncenters, 'rchi_atom',
     &             l_rchi_atom,i_rchi_atom))
     &             call errquit("zora_getv:rchi_atom",0, MA_ERR)
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &             l_rq,i_rq))
     &             call errquit("zora_getv_sf:rq",0, MA_ERR)
c              == delchi ==
               if (.not.MA_Push_Get(mt_dbl, 3*nqpts*nbf_ao,
     &             'delchi_ao', ldelchi_ao, idelchi_ao))
     &             call errquit('zora_getv: delchi_ao',0, MA_ERR)
c              == chi ==
               if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &             'chi_ao', lchi_ao, ichi_ao))
     &             call errquit('zora_getv: chi_ao',0, MA_ERR)
               call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
               call xc_eval_basis(ao_bas_han, 1, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
               call calc_NMRHFine_F1ij(ao_bas_han,geom,ipol,g_dens,     
     &                                 dbl_mb(idelchi_ao),
     &                                 dbl_mb(iqxyz),dbl_mb(iqwght),        
     &                                 nbf_ao,nqpts,ncenters,  
     &                                 ofinite,      ! in: = .true. if Gaussian Nucl. Model of charges requested
     &                                 zetanuc_arr,  ! in: sqrt(zetanuc(i)) i=1,natoms for Gaussian Nuclear Model
     &                                 Knucl,
     &                                 dbl_mb(ihfine(1)),  ! out
     &                                 dbl_mb(ihfine(2)),  ! out
     &                                 dbl_mb(ihfine(3)))  ! out
c              == delete memory ==
               if(.not.MA_pop_stack(lchi_ao))
     &            call errquit("zora_getv: pop chi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(ldelchi_ao))
     &            call errquit("zora_getv: pop delchi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rq))
     &            call errquit("zora_getv: pop rq", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rchi_atom))
     &            call errquit("zora_getv: pop rchi_atom",100,MA_ERR)
            endif ! nqpts
         enddo ! ncube
      end do ! iqsh
c     == delete memory ==
      if(.not.MA_pop_stack(lrqbuf))
     &     call errquit("zora_getv_sf: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("zora_getv_sf: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("zora_getv_sf: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("zora_getv_sf: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("zora_getv_sf: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("zora_getv_sf: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("zora_getv_sf: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("zora_getv_sf: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("zora_getv_sf: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("zora_getv_sf: pop bas_cent_info", 100, MA_ERR)
c
c     == tally up over all the nodes ==
      nij = 0
      do i = 1, nbf_ao
        do j = 1, nbf_ao
         nij = nij + 1
         do k=1,3
          call ga_dgop(msg_excrho,dbl_mb(ihfine(k) +nij-1), 1, '+')
         enddo ! end-loop-k
        enddo ! end-loop-j
      enddo ! end-loop-i
c
c     == pack into a ga all NMR AOs ==
      do i=1,3 
       call ga_zero(g_hfine(i)) 
       call ga_put(g_hfine(i),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(ihfine(i)),nbf_ao)
       call ga_antisymmetrize(g_hfine(i)) 
      enddo ! end-loop-i
      call ga_sync()
c ----- free memory ----------- START
      do i=1,3
      if (.not.ma_free_heap(lhfine(i))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lhfine failed',911,MA_ERR)
      enddo
c ----- free memory ----------- END
      return
      end
