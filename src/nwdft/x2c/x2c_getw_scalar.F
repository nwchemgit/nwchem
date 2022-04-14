c
c     == compute the matrix elements for scalar x2c ==
      subroutine x2c_getw_scalar(rtdb, g_dens,  
     &                          g_w, g_v, 
     &                          nexc, skipnuc, skipcoul)
c
c     modified from zora_getv type routines

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
      integer rtdb      ! [inp] runtime database GA handle
      integer g_dens(2) ! [inp] model pot. density matrix
      integer g_w       ! [out] W-matrix (=pV.p)
      integer g_v(4)    ! scratch
      integer nexc      ! dummy needed for numerical integration code
      logical skipnuc   ! [inp] whether to skip nuclear potential
      logical skipcoul  ! [inp] whether to skip atomic Hartree pot
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
      integer lx2c,ix2c
      integer inntsize,ddblsize,ok

      integer g_tmp

      logical grid_file_rewind
      external grid_file_rewind,ga_antisymmetrize

c     ================================================================

      if (skipnuc .and. skipcoul) then
        if (ga_nodeid().eq.0) write (luout,*)
     &     'x2c_getw_scalar: noting to do. Input error?'
        call errquit ('X2C_GETW_SC skipnuc=T, skipcoul=T', 1,
     &     INPUT_ERR)
      end if

c     debug:
      if (ga_nodeid().eq.0) then
        write (luout,*) 'hello from x2c_getw_scalar'
        write (luout,*) 'skipnuc, skipcoul = ',skipnuc, skipcoul
      end if

c     == allocate memory ==
      if (.not.MA_Push_Get(mt_dbl,nbf_ao*nbf_ao,
     &   'x2ci',lx2c,ix2c)) 
     &   call errquit('x2c_getv_scalar: x2ci',666, MA_ERR)
      
c     == preliminaries ==
      do i= 1, nbf_ao*nbf_ao
        dbl_mb(ix2c+i-1)=0.d0
      enddo
c     == generate the grid ==
      dummy(1) = 0.d0
      dummy(2) = 0.d0

c     g_v is a scratch array in the following call, must
c     have dimension (4)

      call grid_quadv0(rtdb,g_dens,g_v,nexc,rho_n,dummy,tmat)
      if (ga_nodeid().eq.0) write (luout,*) 'rho_n =',rho_n
c
c     == ao basis set info used by xc_eval_basis ==
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('x2c_getv_sc:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('x2c_getv_sc: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('x2c_getv_sc: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
c
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('x2c_getv_sc: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("x2c_getv_sc:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
c
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('x2c_getv_sc: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('x2c_getv_sc: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('x2c_getv_sc: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('x2c_getv_sc: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('x2c_getv_sc: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('x2c_getv_sc: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('x2c_getv_sc: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('x2c_getv_sc: rewinding gridpts?', 0,
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
     &             call errquit("x2c_getv:rchi_atom",0, MA_ERR)
c
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &             l_rq,i_rq))
     &             call errquit("x2c_getv_sc:rq",0, MA_ERR)
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
c              == calculate spin-orbit x2c integrals (pVp) ==
               call calc_x2c_w_scalar(
     &            ao_bas_han, geom, ipol, g_dens, skipnuc, skipcoul,
     &            dbl_mb(idelchi_ao),
     &            dbl_mb(iqxyz),dbl_mb(iqwght),
     &            nbf_ao, nqpts, ncenters,
     &            dbl_mb(ix2c)) ! out

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
     &     call errquit("x2c_getv_sc: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("x2c_getv_sc: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("x2c_getv_sc: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("x2c_getv_sc: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("x2c_getv_sc: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("x2c_getv_sc: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("x2c_getv_sc: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("x2c_getv_sc: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("x2c_getv_sc: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("x2c_getv_sc: pop bas_cent_info", 100, MA_ERR)
c
c     == tally up over all the nodes ==
      nij = 0
      do i = 1, nbf_ao
        do j = 1, nbf_ao
         nij = nij + 1
          call ga_dgop(msg_excrho, dbl_mb(ix2c+nij-1), 1, '+')
        enddo
      enddo

c     == pack into GA and symmetrize ==
c     Note: we need to ADD to the g_w GA because it may alreday hold
c     analytic W integrals with the nuclear potential

      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'getw_tmp',0,0, g_tmp))
     &   call errquit('x2c_getw_sc: error cdreating tmp GA',2, GA_ERR)

      call ga_zero(g_tmp) 
      call ga_put(g_tmp,1,nbf_ao,1,nbf_ao,
     &   dbl_mb(ix2c),nbf_ao)
      call ga_symmetrize(g_tmp)
      call ga_add(1.0d0, g_w, 1.0d0, g_tmp, g_w)

      if (.not. ga_destroy(g_tmp)) call errquit
     &     ('x2c_getw_sc: could not destroy g_tmp', 2, GA_ERR)

      call ga_sync()
c
      if(.not.MA_chop_stack(lx2c))
     &     call errquit("x2c_getv_sc: pop lx2c", 100, MA_ERR)
c
      return
      end
