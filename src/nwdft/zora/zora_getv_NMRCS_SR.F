      subroutine zora_getv_NMRCS_SR12(rtdb, 
     &                                g_dens,        ! in : atomic density 
     &                                xyz_NMRcoords, ! in
     &                                g_dia1,        ! out
     &                                g_dia2,g_dia3, ! out
     &                                g_nmr1,g_nmr2, ! out: munu matrix
     &                                nexc)
c
C$Id$
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
#include "modelpotential.fh"
c
c     == arguments ==
      integer rtdb
      integer g_dens(2)
      integer g_dia1,g_dia2(3),g_dia3(3,3)
      integer g_nmr1(6),g_nmr2(18)
      integer nexc
c
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
      integer ldia1,idia1,
     &        ldia2(3),idia2(3),
     &        ldia3(3,3),idia3(3,3)
      integer lnmr1(3,3),  inmr1(3,3),
     &        lnmr2(3,3,3),inmr2(3,3,3)
      integer ind_tmn(2,3)  
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312 
      integer inntsize,ddblsize,ok
      logical grid_file_rewind
      logical dft_mirrdens_start,dorepl
      integer g_dens0(2),g_densrep(2),ii
      integer nn
      external dft_mirrdens_start
      external grid_file_rewind
c
c     model potential parameters
      character*2 gelem(ncenters)
      double precision gexpo(ncenters,50)
      double precision gcoef(ncenters,50)
c
c     == allocate memory ==
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr1(m,n),inmr1(m,n)))
     &    call errquit('zora_getv_NMR: nmr1',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr1(n,m),inmr1(n,m)))
     &    call errquit('zora_getv_NMR: nmr1',911,MA_ERR)
        do a=1,3
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr2(a,m,n),inmr2(a,m,n)))
     &    call errquit('zora_getv_NMR: nmr2',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr2(a,n,m),inmr2(a,n,m)))
     &    call errquit('zora_getv_NMR: nmr2',911,MA_ERR)
        enddo
       enddo ! end-loop-t
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'ldia1',ldia1,idia1))
     &    call errquit('zora_getv_NMR: dia1',911,MA_ERR)
      do i=1,3
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'ldia1',ldia2(i),idia2(i)))
     &    call errquit('zora_getv_NMR: dia2',911,MA_ERR)
       do j=1,3
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'ldia1',ldia3(i,j),idia3(i,j)))
     &    call errquit('zora_getv_NMR: dia3',911,MA_ERR)
       enddo
      enddo ! end-loop-i
c     == preliminaries ==
      call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(idia1),1)
       do k=1,3
        do j=1,3
           call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(idia3(k,j)),1)
        enddo
       enddo
       do j=1,3
          call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(idia2(j)),1)
       enddo
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)     
          call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(inmr1(m,n)),1)
          call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(inmr1(n,m)),1)
         do a=1,3
            call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(inmr2(a,m,n)),1)
            call dfill(nbf_ao*nbf_ao,0d0,dbl_mb(inmr2(a,n,m)),1)
         enddo ! end-loop-a
       enddo ! end-loop-t
c
c     get zora model potential parameters of geometry
      if (use_modelpotential)
     &  call get_modelpotential_params(rtdb,ncenters,gelem,gexpo,gcoef)
c
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      dorepl = dft_mirrdens_start(g_dens,g_densrep,g_dens0,
     i     ipol)
      call grid_quadv0(rtdb,g_dens,g_dia2,nexc,rho_n,dummy,tmat)
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
                call calc_zora_NMRCS_SR(ao_bas_han,geom,ipol,g_dens,     
     &                                  dbl_mb(ichi_ao),  
     &                                  dbl_mb(idelchi_ao),
     &                                  dbl_mb(iqxyz),dbl_mb(iqwght),        
     &                                  nbf_ao,nqpts,ncenters,     
     &                                  xyz_NMRcoords,
     &                                  use_modelpotential,
     &                                  gexpo, gcoef,
     &                                  dbl_mb(idia1),       ! out 
     &                                  dbl_mb(idia2(1)),    ! out
     &                                  dbl_mb(idia2(2)),    ! out
     &                                  dbl_mb(idia2(3)),    ! out
     &                                  dbl_mb(idia3(1,1)),  ! out
     &                                  dbl_mb(idia3(1,2)),  ! out
     &                                  dbl_mb(idia3(1,3)),  ! out
     &                                  dbl_mb(idia3(2,1)),  ! out
     &                                  dbl_mb(idia3(2,2)),  ! out
     &                                  dbl_mb(idia3(2,3)),  ! out
     &                                  dbl_mb(idia3(3,1)),  ! out
     &                                  dbl_mb(idia3(3,2)),  ! out
     &                                  dbl_mb(idia3(3,3)),  ! out
     &                                  dbl_mb(inmr1(1,2)),  ! out
     &                                  dbl_mb(inmr1(1,3)),  ! out
     &                                  dbl_mb(inmr1(2,1)),  ! out
     &                                  dbl_mb(inmr1(2,3)),  ! out
     &                                  dbl_mb(inmr1(3,1)),  ! out
     &                                  dbl_mb(inmr1(3,2)),  ! out
     &                                  dbl_mb(inmr2(1,1,2)),! out
     &                                  dbl_mb(inmr2(1,1,3)),! out
     &                                  dbl_mb(inmr2(1,2,1)),! out
     &                                  dbl_mb(inmr2(1,2,3)),! out
     &                                  dbl_mb(inmr2(1,3,1)),! out
     &                                  dbl_mb(inmr2(1,3,2)),! out
     &                                  dbl_mb(inmr2(2,1,2)),! out
     &                                  dbl_mb(inmr2(2,1,3)),! out
     &                                  dbl_mb(inmr2(2,2,1)),! out
     &                                  dbl_mb(inmr2(2,2,3)),! out
     &                                  dbl_mb(inmr2(2,3,1)),! out
     &                                  dbl_mb(inmr2(2,3,2)),! out
     &                                  dbl_mb(inmr2(3,1,2)),! out
     &                                  dbl_mb(inmr2(3,1,3)),! out
     &                                  dbl_mb(inmr2(3,2,1)),! out
     &                                  dbl_mb(inmr2(3,2,3)),! out
     &                                  dbl_mb(inmr2(3,3,1)),! out
     &                                  dbl_mb(inmr2(3,3,2)))! out
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
      call ga_sync()
      if(dorepl) then
         do ii=1,ipol
            call util_mirrstop(g_densrep(ii))
            g_dens(ii)=g_dens0(ii)
         enddo
      endif
      nn=nbf_ao*nbf_ao
      call ga_dgop(msg_excrho,dbl_mb(idia1),nn, '+')
      do k=1,3
         call ga_dgop(msg_excrho,dbl_mb(idia2(k)),nn, '+')
         do t=1,3
           call ga_dgop(msg_excrho,dbl_mb(idia3(k,t)),nn,'+')
         enddo ! end-loop-t
      enddo                     ! end-loop-k
      do t=1,3
         m=ind_tmn(1,t)
         n=ind_tmn(2,t) 
         call ga_dgop(msg_excrho,dbl_mb(inmr1(m,n)),nn, '+')
         call ga_dgop(msg_excrho,dbl_mb(inmr1(n,m)),nn, '+')
         do a=1,3
            call ga_dgop(msg_excrho,dbl_mb(inmr2(a,m,n)),nn, '+')
            call ga_dgop(msg_excrho,dbl_mb(inmr2(a,n,m)),nn, '+')
          enddo ! end-loop-a
         enddo ! end-loop-t
c
c     == pack into a ga all NMR AOs ==
       call ga_zero(g_dia1) 
       call ga_put(g_dia1,
     &             1,nbf_ao,1,nbf_ao,dbl_mb(idia1),nbf_ao)
       call ga_symmetrize(g_dia1) 
      do i=1,3 
       call ga_zero(g_dia2(i)) 
       call ga_put(g_dia2(i),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(idia2(i)),nbf_ao)
       call ga_symmetrize(g_dia2(i)) 
       do j=1,3
       call ga_zero(g_dia3(i,j)) 
       call ga_put(g_dia3(i,j),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(idia3(i,j)),nbf_ao)
       call ga_symmetrize(g_dia3(i,j)) 
       enddo ! end-loop-j
      enddo ! end-loop-i
      ac=1
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
       call ga_zero(g_nmr1(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr1(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr1(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      do t=1,3
       m=ind_tmn(2,t)
       n=ind_tmn(1,t) 
       call ga_zero(g_nmr1(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr1(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr1(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      ac=1
      do a=1,3
       do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
       call ga_zero(g_nmr2(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr2(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr2(a,m,n)),nbf_ao)
       ac=ac+1
       enddo ! end-loop-t
       do t=1,3
       m=ind_tmn(2,t)
       n=ind_tmn(1,t) 
       call ga_zero(g_nmr2(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr2(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr2(a,m,n)),nbf_ao)
       ac=ac+1
       enddo ! end-loop-t
      enddo ! end-loop-a
      call ga_sync()
c ----- free memory ----------- START
      if (.not.ma_free_heap(ldia1)) call errquit
     &   ('zora_getv_NMR: ma_free_heap of ldia1 failed',911,MA_ERR)
      do i=1,3
      if (.not.ma_free_heap(ldia2(i))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of ldia2 failed',911,MA_ERR)
       do j=1,3
      if (.not.ma_free_heap(ldia3(i,j))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of ldia3 failed',911,MA_ERR)
       enddo ! end-loop-j
      enddo
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
      if (.not.ma_free_heap(lnmr1(m,n))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr1 failed',911,MA_ERR)
      if (.not.ma_free_heap(lnmr1(n,m))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr1 failed',911,MA_ERR)
       do a=1,3
      if (.not.ma_free_heap(lnmr2(a,m,n))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr2 failed',911,MA_ERR)
      if (.not.ma_free_heap(lnmr2(a,n,m))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr2 failed',911,MA_ERR)
       enddo
      enddo
c ----- free memory ----------- END
      return
      end

      subroutine zora_getv_NMRCS_SR34(rtdb, 
     &                                g_dens,        ! in : atomic density 
     &                                g_nmr,         ! out
     &                                g_nmr3,g_nmr4, ! out: munu matrix
     &                                nexc)
c
C$Id$
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
#include "modelpotential.fh"
c
c     == arguments ==
      integer rtdb
      integer g_dens(2)
      integer g_nmr(3),
     &        g_nmr3(3),g_nmr4(6)
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
      integer lnmr(3)   ,inmr(3),
     &        lnmr3(3)  ,inmr3(3),
     &        lnmr4(3,3),inmr4(3,3)
      integer ind_tmn(2,3)  
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312 
      integer inntsize,ddblsize,ok
      logical grid_file_rewind
      logical dft_mirrdens_start,dorepl
      integer g_dens0(2),g_densrep(2),ii
      integer nn
      external grid_file_rewind
c
c     model potential parameters
      character*2 gelem(ncenters)
      double precision gexpo(ncenters,50)
      double precision gcoef(ncenters,50)
c
c     == allocate memory ==
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr4(m,n),inmr4(m,n)))
     &    call errquit('zora_getv_NMR: nmr4',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr4(n,m),inmr4(n,m)))
     &    call errquit('zora_getv_NMR: nmr4',911,MA_ERR)
       enddo ! end-loop-t
      do i=1,3
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr(i),inmr(i)))
     &    call errquit('zora_getv_NMR: nmr',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nbf_ao*nbf_ao,
     &                      'lnmr',lnmr3(i),inmr3(i)))
     &    call errquit('zora_getv_NMR: nmr3',911,MA_ERR)
      enddo ! end-loop-i
c     == preliminaries ==
       do j=1,3
         do i= 1, nbf_ao*nbf_ao
          dbl_mb(inmr(j)+i-1) =0.d0
          dbl_mb(inmr3(j)+i-1)=0.d0
         enddo  
       enddo
       do t=1,3
        m=ind_tmn(1,t)
        n=ind_tmn(2,t)     
         do i= 1, nbf_ao*nbf_ao
          dbl_mb(inmr4(m,n)+i-1)=0.d0
          dbl_mb(inmr4(n,m)+i-1)=0.d0
         enddo ! end-loop-i   
       enddo ! end-loop-t
c
c     get zora model potential parameters of geometry
      if (use_modelpotential)
     &  call get_modelpotential_params(rtdb,ncenters,gelem,gexpo,gcoef)
c
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      dorepl = dft_mirrdens_start(g_dens,g_densrep,g_dens0,
     i     ipol)
      call grid_quadv0(rtdb,g_dens,g_nmr(1),nexc,rho_n,dummy,tmat)
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
                call calc_NMRCS_SR_F1ij(ao_bas_han,geom,ipol,g_dens,     
     &                                  dbl_mb(ichi_ao),  
     &                                  dbl_mb(idelchi_ao),
     &                                  dbl_mb(iqxyz),dbl_mb(iqwght),        
     &                                  nbf_ao,nqpts,ncenters,     
     &                                  use_modelpotential,
     &                                  gexpo, gcoef,
     &                                  dbl_mb(inmr(1)),     ! out
     &                                  dbl_mb(inmr(2)),     ! out
     &                                  dbl_mb(inmr(3)),     ! out
     &                                  dbl_mb(inmr3(1)),    ! out
     &                                  dbl_mb(inmr3(2)),    ! out
     &                                  dbl_mb(inmr3(3)),    ! out
     &                                  dbl_mb(inmr4(1,2)),  ! out
     &                                  dbl_mb(inmr4(1,3)),  ! out
     &                                  dbl_mb(inmr4(2,1)),  ! out
     &                                  dbl_mb(inmr4(2,3)),  ! out
     &                                  dbl_mb(inmr4(3,1)),  ! out
     &                                  dbl_mb(inmr4(3,2)))  ! out
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
      call ga_sync()
      if(dorepl) then
         do ii=1,ipol
            call util_mirrstop(g_densrep(ii))
            g_dens(ii)=g_dens0(ii)
         enddo
      endif
      nn=nbf_ao*nbf_ao
      do k=1,3
         call ga_dgop(msg_excrho,dbl_mb(inmr(k)), nn, '+')
         call ga_dgop(msg_excrho,dbl_mb(inmr3(k)), nn, '+')
      enddo                     ! end-loop-k
      do t=1,3
         m=ind_tmn(1,t)
         n=ind_tmn(2,t) 
         call ga_dgop(msg_excrho,dbl_mb(inmr4(m,n)), nn, '+')
         call ga_dgop(msg_excrho,dbl_mb(inmr4(n,m)), nn, '+')
      enddo                     ! end-loop-t
c
c     == pack into a ga all NMR AOs ==
      do i=1,3 
       call ga_zero(g_nmr(i)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr(i),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr(i)),nbf_ao)
       call ga_symmetrize(g_nmr(i)) 
       call ga_zero(g_nmr3(i)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr3(i),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr3(i)),nbf_ao)
      enddo ! end-loop-i
      ac=1
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
       call ga_zero(g_nmr4(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr4(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr4(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      do t=1,3
       m=ind_tmn(2,t)
       n=ind_tmn(1,t) 
       call ga_zero(g_nmr4(ac)) 
       if(ga_nodeid().eq.0)
     G      call ga_put(g_nmr4(ac),
     &             1,nbf_ao,1,nbf_ao,dbl_mb(inmr4(m,n)),nbf_ao)
       ac=ac+1
      enddo ! end-loop-t
      call ga_sync()

c ----- free memory ----------- START
      do i=1,3
      if (.not.ma_free_heap(lnmr(i))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr failed',911,MA_ERR)
      if (.not.ma_free_heap(lnmr3(i))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr3 failed',911,MA_ERR)
      enddo
      do t=1,3
       m=ind_tmn(1,t)
       n=ind_tmn(2,t) 
      if (.not.ma_free_heap(lnmr4(m,n))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr4 failed',911,MA_ERR)
      if (.not.ma_free_heap(lnmr4(n,m))) call errquit
     &   ('zora_getv_NMR: ma_free_heap of lnmr4 failed',911,MA_ERR)
      enddo ! end-loop-t
c ----- free memory ----------- END
      return
      end

