c
c     == compute the matrix elements for the zora corrections ==
c     == scalar and spin-orbit components (x,y,z)
      subroutine zora_getv_so(rtdb, g_dens,  
     &                          g_zora_sf, g_zora_scale_sf, ! out
     &                          g_zora_so, g_zora_scale_so, ! out
     &                          nexc)
c
C$Id$
c Modified from zora_getv_so() by FA 10-31-10

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
      integer g_zora_sf(2)
      integer g_zora_scale_sf(2)
      integer g_zora_so(3)
      integer g_zora_scale_so(3)
      integer nexc
c
c     == local variables ==
      integer i,j,k,ind,nij
      double precision rho_n
      double precision tmat
      double precision dummy(2)
      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwght,iqwght,nqpts,ncube,
     &     ictr_buf,iqpts
      double precision rad,ke
      integer lbas_cset_info, ibas_cset_info,
     &     lbas_cent_info, ibas_cent_info,
     &     ldocset, idocset,
     &     l_rchi_atom,i_rchi_atom,
     &     l_rq,i_rq,lniz, iniz,
     &     lchi_ao, ichi_ao,
     &     ldelchi_ao, idelchi_ao
      integer lzora(4),izora(4),lscal(4),iscal(4)
      integer inntsize,ddblsize,ok
      logical grid_file_rewind
      external grid_file_rewind,ga_antisymmetrize
c
c     == allocate memory ==
      do i=1,4
       if (.not.MA_Push_Get(mt_dbl,nbf_ao*nbf_ao,
     &                      'zora0',lzora(i),izora(i))) 
     &    call errquit('zora_getv_so: zorai',0, MA_ERR)
       if (.not.MA_Push_Get(mt_dbl,nbf_ao*nbf_ao,
     &                     'scal0',lscal(i),iscal(i))) 
     &    call errquit('zora_getv_so: scali',0, MA_ERR)
      enddo
c     == preliminaries ==
      do i= 1, nbf_ao*nbf_ao
       do j=1,4
         dbl_mb(izora(j)+i-1)=0.d0
         dbl_mb(iscal(j)+i-1)=0.d0
       enddo
      enddo
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
      call grid_quadv0(rtdb,g_dens,g_zora_sf,nexc,rho_n,dummy,tmat)
c
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
c
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
c
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
c
c              == calculate scalar and spin-orbit zora contribution ==
               call calc_zora_so(ao_bas_han, geom, ipol, g_dens,
     &                              dbl_mb(idelchi_ao),
     &                              dbl_mb(iqxyz),dbl_mb(iqwght),
     &                              nbf_ao, nqpts, ncenters,
     &                              dbl_mb(izora(1)),dbl_mb(iscal(1)), ! out
     &                              dbl_mb(izora(2)),dbl_mb(iscal(2)), ! out
     &                              dbl_mb(izora(3)),dbl_mb(iscal(3)), ! out
     &                              dbl_mb(izora(4)),dbl_mb(iscal(4))) ! out
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
         do k=1,4
          call ga_dgop(msg_excrho, dbl_mb(izora(k)+nij-1), 1, '+')
          call ga_dgop(msg_excrho, dbl_mb(iscal(k)+nij-1), 1, '+')
         enddo
        enddo
      enddo
c     == pack into a ga ==
c     == scalar contribution ==
      do i=1,ipol
       call ga_zero(g_zora_sf(i))
       call ga_put(g_zora_sf(i),1,nbf_ao,1,nbf_ao,
     &             dbl_mb(izora(1)),nbf_ao)
       call ga_symmetrize(g_zora_sf(i))
       call ga_zero(g_zora_scale_sf(i))
       call ga_put(g_zora_scale_sf(i),1,nbf_ao,1,nbf_ao,
     &             dbl_mb(iscal(1)),nbf_ao)
       call ga_symmetrize(g_zora_scale_sf(i))
      enddo
c     == spin-orbit contributions ==
      ind=4 ! = 4,3,2 = z,y,x in g_zora_so(i)
      do i=1,3 ! = z,y,x
       call ga_zero(g_zora_so(i)) 
       call ga_put(g_zora_so(i),1,nbf_ao,1,nbf_ao,
     &             dbl_mb(izora(ind)),nbf_ao)
       call ga_antisymmetrize(g_zora_so(i))
       call ga_zero(g_zora_scale_so(i)) 
       call ga_put(g_zora_scale_so(i),1,nbf_ao,1,nbf_ao,
     &             dbl_mb(iscal(ind)),nbf_ao)
       call ga_antisymmetrize(g_zora_scale_so(i))
       ind=ind-1
      enddo
      call ga_sync()
c
      if(.not.MA_chop_stack(lscal))
     &     call errquit("zora_getv_so: pop lscal", 100, MA_ERR)
      if(.not.MA_chop_stack(lzora))
     &     call errquit("zora_getv_so: pop lzora", 100, MA_ERR)
c
      return
      end
