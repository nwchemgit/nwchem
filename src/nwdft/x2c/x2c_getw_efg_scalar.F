c
c     == compute the matrix elements for EFG operator in 
c     == small component basis, for x2c

      subroutine x2c_getw_efg_scalar(rtdb, g_dens, origin, 
     &                          g_w,
     &                          nexc)
c
c     modified from zora_getv_so

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
      double precision origin(3) ! [in] EFG center
      integer g_w(6) ! [out]
      integer nexc
c
c     == local variables ==
      integer i,j,k,iefg,ind,nij
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
      integer lx2c ,ix2c
      integer dbg
      integer inntsize,ddblsize,ok
      logical grid_file_rewind
      external grid_file_rewind,ga_antisymmetrize

      logical me
c
      me = (ga_nodeid().eq.0)
      dbg=0

      if (dbg>0 .and. me) write (luout,*) 'entering getw_efg'

c     == allocate memory ==
c$$$      do i=1,4
c$$$        do iefg = 1,6
c$$$          if (.not.MA_Push_Get(mt_dbl,nbf_ao*nbf_ao,
c$$$     &       'x2ci',lx2c(i,iefg),ix2c(i,iefg))) 
c$$$     &       call errquit('x2c_getv_so: x2ci',i + iefg*10, MA_ERR)
c$$$        end do
c$$$      enddo

      if (.not.MA_Push_Get(mt_dbl, 6*(nbf_ao**2),
     &   'x2ci',lx2c,ix2c)) 
     &   call errquit('x2c_getv_so: x2ci', 6*(nbf_ao**2), MA_ERR)

c     == set output MA arrays to zero

c$$$      do i = 1,4
c$$$        do iefg = 1,6
c$$$          call dfill(nbf_ao**2, 0.0d0, dbl_mb(ix2c(i,iefg)),  1)
c$$$        end do
c$$$      end do

      call dfill(6*(nbf_ao**2), 0.0d0, dbl_mb(ix2c),  1)

c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0
c     g_w is a dummy array in the following call
      if (dbg>0 .and.me)  write (luout,*) 'calling grid_quadv0'
      call grid_quadv0(rtdb,g_dens,g_w,nexc,rho_n,dummy,tmat)
      if (me) write (luout,*) 'rho_n =',rho_n
c
c     == ao basis set info used by xc_eval_basis ==
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('x2c_getv_so:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('x2c_getv_so: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('x2c_getv_so: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
c
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('x2c_getv_so: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("x2c_getv_so:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
c
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('x2c_getv_so: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('x2c_getv_so: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('x2c_getv_so: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('x2c_getv_so: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('x2c_getv_so: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('x2c_getv_so: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('x2c_getv_so: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('x2c_getv_so: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
c     == loop over records in the grid file ==
      do iqsh = 1, n_rec_in_file
c
c       == define current range of rad. shells and integration center ==
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
     &             call errquit("x2c_getv:rchi_atom",0, MA_ERR)
c
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &             l_rq,i_rq))
     &             call errquit("x2c_getv_so:rq",0, MA_ERR)
c
c              == delchi ==
               if (.not.MA_Push_Get(mt_dbl, 3*nqpts*nbf_ao,
     &             'delchi_ao', ldelchi_ao, idelchi_ao))
     &             call errquit('x2c_getv: delchi_ao',0, MA_ERR)
c
c              == chi ==
               if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &             'chi_ao', lchi_ao, ichi_ao))
     &             call errquit('x2c_getv: chi_ao',0, MA_ERR)
               call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
               call xc_eval_basis(ao_bas_han, 1, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
c
c              calculate EFG integrals
               if(dbg>0 .and. me) write (luout,*) 'calling wefg_sc'
               call calc_x2c_wefg_scalar(ao_bas_han, geom, ipol, g_dens,
     &            dbl_mb(idelchi_ao),
     &            dbl_mb(iqxyz),dbl_mb(iqwght),
     &            nbf_ao, nqpts, ncenters,
     &            origin,
     &            dbl_mb(ix2c)) ! out
               if(dbg>0 .and. me) write (luout,*) 'return from wefg_sc'


c              == delete memory ==
               if(.not.MA_pop_stack(lchi_ao))
     &            call errquit("x2c_getv: pop chi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(ldelchi_ao))
     &            call errquit("x2c_getv: pop delchi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rq))
     &            call errquit("x2c_getv: pop rq", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rchi_atom))
     &            call errquit("x2c_getv: pop rchi_atom",100,MA_ERR)
            endif ! nqpts
         enddo ! ncube
      end do ! iqsh
c
c     == delete memory ==
      if(.not.MA_pop_stack(lrqbuf))
     &     call errquit("x2c_getv_so: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("x2c_getv_so: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("x2c_getv_so: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("x2c_getv_so: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("x2c_getv_so: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("x2c_getv_so: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("x2c_getv_so: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("x2c_getv_so: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("x2c_getv_so: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("x2c_getv_so: pop bas_cent_info", 100, MA_ERR)
c
c     == tally up over all the nodes ==

      nij = 0
      do i = 1, nbf_ao
        do j = 1, nbf_ao
          do iefg = 1,6
            nij = nij + 1
            call ga_dgop(msg_excrho,dbl_mb(ix2c +nij-1),1,'+')
          end do
        enddo
      enddo

c     == pack into GA and symmetrize ==

      do iefg = 1,6
        call ga_zero(g_w(iefg)) 
c        call ga_put(g_w(i,iefg),1,nbf_ao,1,nbf_ao,
c     &     dbl_mb(ix2c(i,iefg)),nbf_ao)
        call ga_put(g_w(iefg),1,nbf_ao,1,nbf_ao,
     &     dbl_mb(ix2c + ((iefg-1)*(nbf_ao**2)) ),nbf_ao)
        call ga_symmetrize(g_w(iefg))
      end do
      call ga_sync()
c
      if(.not.MA_chop_stack(lx2c))
     &     call errquit("x2c_getv_so: pop lx2c", 100, MA_ERR)

      if (dbg>0 .and. me) write (luout,*) 'leaving getw_efg'
c
      return
      end
