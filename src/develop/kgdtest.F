c $Id: kgdtest.F,v 1.8 2000-02-04 01:30:24 mg201 Exp $
      logical function kgdtest (rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "stdio.fh"
c::passed
      integer rtdb          ! rtdb handle
c::local
      integer me
      integer kgdtask, kgd_tmp
c
      call ga_sync()
      call ga_sync()
      kgdtask = 0
      if (rtdb_get(rtdb,'kgdtask',MT_INT,1,kgd_tmp))
     &    kgdtask = kgd_tmp
c
      call ga_sync()
      call ga_sync()
      me = ga_nodeid()
c
      if (kgdtask.eq.0) then    !...................................   0
        if (me.eq.0) then
          write(luout,*)' default kgdtest task '
          write(luout,*)' test use of kgdtest! '
        endif
        kgdtest = .true.
      else if (kgdtask.eq.1) then !.................................   1
        if (me.eq.0) write(luout,*)
     &      ' kgdtest task 1: relativistic integral test'
        call kgdtest_rel1e(rtdb)
        kgdtest = .true.
      end if
      end
      subroutine kgdtest_rel1e(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      logical int_normalize
      external int_normalize
c
c test relativistic one-electron integrals
c
      integer rtdb
      integer geom,basis, basis_id
      integer nshell, memscr, membuf
      integer h_scr, k_scr, h_buf, k_buf
      integer ish, jsh, ucont
      integer li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
C      integer nint_out
      integer ihi,jhi
      logical status
      character*255 mo_basis, geom_name
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      if (.not.context_rtdb_match(rtdb,'ao basis',mo_basis))
     &    mo_basis = 'ao basis'
      if (.not.context_rtdb_match(rtdb,'geometry',geom_name))
     &    geom_name = 'geometry'
c
      if(.not.geom_create(geom,geom_name))call errquit
     &    ('kgdtest_rel1e: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('kgdtest_rel1e: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('kgdtest_rel1e: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('kgdtest_rel1e: basis load ',911)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(basis,rtdb)) call errquit
     &    ('kgdtest_rel1e: error normalizing ',911)
c
      call int_init(rtdb,1,basis)
      memscr = 100 000
      membuf = 1000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911)
      if (.not.ma_push_get(mt_dbl,membuf*3,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911)
c
      do ish = 1,nshell
        do jsh = 1,ish
          write(6,*)' ============= shells <',ish,'|',jsh,'>',
     &        '==================== start =========='
          write(6,*)' '
          
          ucont = (sf_ibs_cn2ucn(ish,basis_id))
          Li      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          i_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          i_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          i_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          i_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          i_cent  = (sf_ibs_cn2ce(ish,basis_id))
          i_geom  = ibs_geom(basis_id)
          ihi = i_gen*(Li+1)*(Li+2)/2
c
          ucont = (sf_ibs_cn2ucn(jsh,basis_id))
          Lj      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          j_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          j_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          j_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          j_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          j_cent  = (sf_ibs_cn2ce(jsh,basis_id))
          j_geom  = ibs_geom(basis_id)
          jhi = j_gen*(Lj+1)*(Lj+2)/2
*
*   Calculate overlap and kinetic energy integrals
*
      call int_hf1sp(
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,i_cent,
     &        coords(1,j_cent,j_geom),
     &        dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,j_cent,
     &        coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
     &        dbl_mb(k_buf),dbl_mb(k_buf+1000),dbl_mb(k_buf+2000),
     &        .true.,.true.,.false.,.false.,
     &        .false.,dbl_mb(k_scr),memscr,'kgdtest_rel1e')
          call ecp_matpr (dbl_mb(k_buf),1,jhi,1,ihi,1,jhi,1,ihi,
     &        'overlap integrals','E',120,8)
          call ecp_matpr (dbl_mb(k_buf+1000),1,jhi,1,ihi,1,jhi,1,ihi,
     &        'kinetic integrals','E',120,8)
C     modified metric
          call rel_onel (
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,
     &        coords(1,j_cent,j_geom),
     &        dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
     &        dbl_mb(k_buf),dbl_mb(k_buf+1000),dbl_mb(k_buf+2000),
     &        membuf,.true.,.false.,.false.,.false.,.true.,.false.,
     &        .false.,.false.,dbl_mb(k_scr),memscr,0,1)
          call ecp_matpr (dbl_mb(k_buf),1,jhi,1,ihi,1,jhi,1,ihi,
     &        'modified overlap integrals','E',120,8)
C     modified kinetic energy
          call rel_onel (
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,
     &        coords(1,j_cent,j_geom),
     &        dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
     &        dbl_mb(k_buf),dbl_mb(k_buf+1000),dbl_mb(k_buf+2000),
     &        membuf,.false.,.true.,.false.,.false.,.true.,.false.,
     &        .false.,.false.,dbl_mb(k_scr),memscr,0,1)
          call ecp_matpr (dbl_mb(k_buf+1000),1,jhi,1,ihi,1,jhi,1,ihi,
     &        'modified kinetic integrals','E',120,8)
C     modified potential energy
          call rel_onel (
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,
     &        coords(1,j_cent,j_geom),
     &        dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
     &        dbl_mb(k_buf),dbl_mb(k_buf+1000),dbl_mb(k_buf+2000),
     &        membuf,.false.,.false.,.true.,.false.,.true.,.false.,
     &        .false.,.false.,dbl_mb(k_scr),memscr,0,1)
          call ecp_matpr (dbl_mb(k_buf+2000),1,jhi,1,ihi,1,jhi,1,ihi,
     &        'modified potential integrals','E',120,8)
          write(6,*)' ============= shells <',ish,'|',jsh,'>',
     &          '====================  end  =========='
          write(6,*)' '
        enddo
      enddo
c      
      call int_terminate()
      status = ma_pop_stack(h_buf)
      status = status.and.ma_pop_stack(h_scr)
      if (.not.status) call errquit('pop failed',911)
      status = bas_destroy(basis)
      status = status.and.geom_destroy(geom)
      if (.not.status) call errquit('b/g destroy failed',911)
      return
      end
