c $Id: kgdtest.F,v 1.10 2000-05-26 20:05:42 mg201 Exp $
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
      write (6,*) 'in kgdtest'
      call util_flush(6)
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
      else if (kgdtask.eq.2) then !.................................   2
        if (me.eq.0) write(luout,*)
     &      ' kgdtest task 2: general contraction test'
        call kgdtest_gencon(rtdb)
        kgdtest = .true.
      else if (kgdtask.eq.3) then !.................................   3
        if (me.eq.0) write(luout,*)
     &      ' kgdtest task 3: so-ecp integral test'
        call kgdtest_soecp(rtdb)
        kgdtest = .true.
      else if (kgdtask.eq.4) then !.................................   4
        if (me.eq.0) write(luout,*)
     &      ' kgdtest task 4: rel2e integral test'
        call kgdtest_rel2e(rtdb)
        kgdtest = .true.
      end if
      end
************************************************************************
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
      if (.not.int_normalize(rtdb,basis)) call errquit
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
************************************************************************
      subroutine kgdtest_gencon(rtdb)
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
#include "stdio.fh"
c
      logical int_normalize
      external int_normalize
      integer idamax
      external idamax
c
c test general contraction in McMD 2e integrals.
c
      integer rtdb
      integer geom,basis, basis_id
      integer nshell, memscr, membuf
      integer h_scr, k_scr, h_buf, k_buf
      integer ish, jsh, ksh,lsh, ucont
      integer li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer lk, k_prim, k_gen, k_iexp, k_icfp, k_cent, k_geom
      integer ll, l_prim, l_gen, l_iexp, l_icfp, l_cent, l_geom
      integer Nints
      integer i_eri,i_c,j_c,k_c,l_c,ic,jc,kc,lc,i_seg
      integer ihi,jhi,khi,lhi,i2,j2,k2,l2,n_cart,n_cont
      character*255 mo_basis, geom_name
      double precision difmax
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
     &    ('kgdtest_gencon: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('kgdtest_gencon: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('kgdtest_gencon: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('kgdtest_gencon: basis load ',911)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(rtdb,basis)) call errquit
     &    ('kgdtest_gencon: error normalizing ',911)
c
      call int_init(rtdb,1,basis)
      memscr = 1 000 000
      membuf = 10 000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911)
      if (.not.ma_push_get(mt_dbl,membuf,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911)
c
      difmax = 0.0d0
C      do ish = 1,nshell
      ish = 35
        ucont = (sf_ibs_cn2ucn(ish,basis_id))
        Li      = infbs_cont(CONT_TYPE ,ucont,basis_id)
        i_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
        i_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
        i_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
        i_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
        i_cent  = (sf_ibs_cn2ce(ish,basis_id))
        i_geom  = ibs_geom(basis_id)
        i2 = (Li+1)*(Li+2)/2
        ihi = i_gen*i2
c
C        do jsh = 1,ish
        jsh = 11
          ucont = (sf_ibs_cn2ucn(jsh,basis_id))
          Lj      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          j_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          j_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          j_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          j_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          j_cent  = (sf_ibs_cn2ce(jsh,basis_id))
          j_geom  = ibs_geom(basis_id)
          j2 = (Lj+1)*(Lj+2)/2
          jhi = j_gen*j2
c
C          do ksh = 1,ish
          ksh = 2
            ucont = (sf_ibs_cn2ucn(ksh,basis_id))
            Lk      = infbs_cont(CONT_TYPE ,ucont,basis_id)
            k_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
            k_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
            k_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
            k_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
            k_cent  = (sf_ibs_cn2ce(ish,basis_id))
            k_geom  = ibs_geom(basis_id)
            k2 = (Lk+1)*(Lk+2)/2
            khi = k_gen*k2
c
C            do lsh = 1,ksh
            lsh = 2
              ucont = (sf_ibs_cn2ucn(lsh,basis_id))
              Ll      = infbs_cont(CONT_TYPE ,ucont,basis_id)
              l_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
              l_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
              l_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
              l_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
              l_cent  = (sf_ibs_cn2ce(lsh,basis_id))
              l_geom  = ibs_geom(basis_id)
              l2 = (Ll+1)*(Ll+2)/2
              lhi = l_gen*l2
              
              n_cart = i2*j2*k2*l2
              n_cont = i_gen*j_gen*k_gen*l_gen
              Nints = ihi*jhi*khi*lhi
              write (LuOut,*) ish,Li,i_prim,i_gen,i_cent
              write (LuOut,*) jsh,Lj,j_prim,j_gen,j_cent
              write (LuOut,*) ksh,Lk,k_prim,k_gen,k_cent
              write (LuOut,*) lsh,Ll,l_prim,l_gen,l_cent
C              if (n_cont .gt. 1 .and. Li.ne.Lj .and. Lk.ne.Ll) then
C                write(luout,*)' '
C                write(luout,*)' ============= (',
C     &              ish,':',Li,',',
C     &              jsh,':',Lj,'|',
C     &              ksh,':',Lk,',',
C     &              lsh,':',Ll,')',
C     &              '===================='
C                write(luout,*)' '
                write (luout,*) 'Nints',Nints

                call hf2(
     &            coords(1,i_cent,i_geom),
     &              dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &              dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,
     &              coords(1,j_cent,j_geom),
     &              dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &              dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,
     &              coords(1,k_cent,k_geom),
     &              dbl_mb(mb_exndcf(k_iexp,basis_id)),
     &              dbl_mb(mb_exndcf(k_icfp,basis_id)),k_prim,k_gen,Lk,
     &              coords(1,l_cent,l_geom),
     &              dbl_mb(mb_exndcf(l_iexp,basis_id)),
     &              dbl_mb(mb_exndcf(l_icfp,basis_id)),l_prim,l_gen,Ll,
     &              dbl_mb(k_buf),Nints,.false.,.false.,.false.,
     &              .false.,dbl_mb(k_scr),memscr)
              
                i_eri = k_buf+Nints
                i_seg = i_eri
                i_c = i_icfp
                do ic = 1,i_gen
                  j_c = j_icfp
                  do jc = 1,j_gen
                    k_c = k_icfp
                    do kc = 1,k_gen
                      l_c = l_icfp
                      do lc = 1,l_gen
C                        write (luout,*) 'ic,jc,kc,lc',ic,jc,kc,lc
                        call hf2(
     &                      coords(1,i_cent,i_geom),
     &                      dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &                      dbl_mb(mb_exndcf(i_c,basis_id)),
     &                      i_prim,1,Li,
     &                      coords(1,j_cent,j_geom),
     &                      dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &                      dbl_mb(mb_exndcf(j_c,basis_id)),
     &                      j_prim,1,Lj,
     &                      coords(1,k_cent,k_geom),
     &                      dbl_mb(mb_exndcf(k_iexp,basis_id)),
     &                      dbl_mb(mb_exndcf(k_c,basis_id)),
     &                      k_prim,1,Lk,
     &                      coords(1,l_cent,l_geom),
     &                      dbl_mb(mb_exndcf(l_iexp,basis_id)),
     &                      dbl_mb(mb_exndcf(l_c,basis_id)),
     &                      l_prim,1,Ll,
     &                      dbl_mb(i_eri),Nints,.false.,.false.,.false.,
     &                      .false.,dbl_mb(k_scr),memscr)

                        i_eri = i_eri+i2*j2*k2*l2
                        l_c = l_c+l_prim
                      end do
                      k_c = k_c+k_prim
                    end do
                    j_c = j_c+j_prim
                  end do
                  i_c = i_c+i_prim
                end do

                call reorder_2eints(dbl_mb(i_eri),dbl_mb(i_seg),
     &              l2,l_gen,k2,k_gen,j2,j_gen,i2,i_gen)
                call dcopy (Nints,dbl_mb(i_eri),1,dbl_mb(i_seg),1)
                call daxpy(Nints,-1.0d0,dbl_mb(k_buf),1,dbl_mb(i_eri),1)
                ic = idamax(Nints,dbl_mb(i_eri),1)-1
                write (luout,*) 'Maximum difference',dbl_mb(i_eri+ic)
                difmax = max(difmax,abs(dbl_mb(i_eri+ic)))
                if (abs(dbl_mb(i_eri+ic)) .gt. 1.0d-12) then
                  call ecp_matpr (dbl_mb(k_buf),1,khi*lhi,1,ihi*jhi,
     &                1,khi*lhi,1,ihi*jhi,'General contraction',
     &                'E',120,6)
                  call ecp_matpr (dbl_mb(i_seg),1,khi*lhi,1,ihi*jhi,
     &                1,khi*lhi,1,ihi*jhi,'Segmented contraction',
     &                'E',120,6)
                end if
C              end if
C            end do
C          end do
C        end do
C      end do
      write (luout,*) 'Maximum difference of all blocks',difmax

      end
      subroutine reorder_2eints(bnew,bold,lc,lg,kc,kg,jc,jg,ic,ig)
      double precision bnew(lc,lg,kc,kg,jc,jg,ic*ig)
      double precision bold(lc,kc,jc,ic,lg,kg,jg*ig)
      integer lc,lg,kc,kg,jc,jg,ic,ig
      integer i,j,k,l,ii,jj,kk,ll,inew,iold
      iold = 0
      inew = 0
      do i = 1,ig
        do j = 1,jg
          iold = iold+1
          do k = 1,kg
            do l = 1,lg
              do ii = 1,ic
                do jj = 1,jc
                  do kk = 1,kc
                    do ll = 1,lc
                      bnew(ll,l,kk,k,jj,j,inew+ii) =
     &                    bold(ll,kk,jj,ii,l,k,iold)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
        inew = inew+ic
      end do
      end

************************************************************************
      subroutine kgdtest_soecp(rtdb)
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
#include "ecp_nwc.fh"
#include "stdio.fh"
c
      logical int_normalize
      external int_normalize
      integer idamax
      external idamax
c
c test spin-orbit integrals
c
      integer rtdb
      integer geom,basis, basis_id
      integer nshell, memscr, membuf
      integer h_scr, k_scr, h_buf, k_buf
      integer ish, jsh, ucont
      integer li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_geom
      integer lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_geom
      integer ihi,jhi,i2,j2
      integer ibug,n_blk
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
     &    ('kgdtest_gencon: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('kgdtest_gencon: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('kgdtest_gencon: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('kgdtest_gencon: basis load ',911)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(rtdb,basis)) call errquit
     &    ('kgdtest_gencon: error normalizing ',911)
c
      call int_init(rtdb,1,basis)
      memscr = 1 000 000
      membuf = 10 000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911)
      if (.not.ma_push_get(mt_dbl,membuf,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911)
c
      if (.not. rtdb_get(rtdb,'ecp:ibug',mt_int,1,ibug))
     &    ibug = 3
      if (.not. rtdb_get(rtdb,'ecp:n_blk',mt_int,1,n_blk))
     &    n_blk = 3
c
      do ish = 1,nshell
        ucont = (sf_ibs_cn2ucn(ish,basis_id))
        Li      = infbs_cont(CONT_TYPE ,ucont,basis_id)
        i_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
        i_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
        i_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
        i_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
        i_cent  = (sf_ibs_cn2ce(ish,basis_id))
        i_geom  = ibs_geom(basis_id)
        i2 = (Li+1)*(Li+2)/2
        ihi = i_gen*i2
c
        do jsh = 1,ish
          ucont = (sf_ibs_cn2ucn(jsh,basis_id))
          Lj      = infbs_cont(CONT_TYPE ,ucont,basis_id)
          j_prim  = infbs_cont(CONT_NPRIM,ucont,basis_id)
          j_gen   = infbs_cont(CONT_NGEN ,ucont,basis_id)
          j_iexp  = infbs_cont(CONT_IEXP ,ucont,basis_id)
          j_icfp  = infbs_cont(CONT_ICFP ,ucont,basis_id)
          j_cent  = (sf_ibs_cn2ce(jsh,basis_id))
          j_geom  = ibs_geom(basis_id)
          j2 = (Lj+1)*(Lj+2)/2
          jhi = j_gen*j2

          write(luout,*)' '
          write(luout,*)' ============= ',ish,jsh,
     &        '===================='
          write(luout,*)' '

          call ecp_integral (
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(i_icfp,basis_id)),i_prim,i_gen,Li,i_cent,
     &        coords(1,j_cent,j_geom),
     &        dbl_mb(mb_exndcf(j_iexp,basis_id)),
     &        dbl_mb(mb_exndcf(j_icfp,basis_id)),j_prim,j_gen,Lj,j_cent,
     &        dbl_mb(k_xyzecp),
     &        dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &        int_mb(k_ecp_nprim_c),
     &        int_mb(k_ecp_ncoef_c), ! new name is n_colc_C
     &        int_mb(k_ecp_ind_z),
     &        int_mb(k_ecp_ind_c),
     &        n_zeta_c,
     &        n_zeta_c,
     &        int_mb(k_ecp_l_c),
     &        int_mb(k_ecp_lip), 
     &        n_ecp,l_ecp,
     &        0,
     &        dbl_mb(k_ecp_c2s),mem_c2s,
     &        dbl_mb(k_buf),ihi*jhi,n_blk, 
     &        .false.,dbl_mb(k_scr),memscr,
     &        ibug)               

        end do
      end do
      call int_terminate()
      end
************************************************************************
      subroutine kgdtest_rel2e(rtdb)
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
#include "apiP.fh"
#include "rel_nwc.fh"
#include "stdio.fh"
c
      logical int_normalize
      external int_normalize
      integer idamax
      external idamax
c
c test relativistic 2e integrals for correct magnitude
c
      integer rtdb
      integer geom,basis, basis_id
      integer nshell, memscr, membuf
      integer h_scr, k_scr, h_buf, k_buf
      integer ish, jsh, ksh,lsh, ucont
      integer li, ipr, igc, iex, icf, icfs, ice, igm
      integer lj, jpr, jgc, jex, jcf, jcfs, jce, jgm
      integer lk, kpr, kgc, kex, kcf, kcfs, kce, kgm
      integer ll, lpr, lgc, lex, lcf, lcfs, lce, lgm
      integer Nints
      integer i_eri,i,abas,sbas
      integer ihi,jhi,khi,lhi,i2,j2,k2,l2,n_cart,n_cont
      character*255 mo_basis, geom_name
      double precision difmax,errmax,dif
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
     &    ('kgdtest_rel2e: geom create error',911)
      if(.not.bas_create(basis,mo_basis))call errquit
     &    ('kgdtest_rel2e: basis create error',911)
c
      if(.not.geom_rtdb_load(rtdb,geom,geom_name)) call errquit
     &    ('kgdtest_rel2e: geom load ',911)
      if(.not.bas_rtdb_load(rtdb,geom,basis,mo_basis)) call errquit
     &    ('kgdtest_rel2e: basis load ',911)
c
      basis_id = basis + BASIS_HANDLE_OFFSET
      nshell = ncont_tot_gb(basis_id)
      if (.not.int_normalize(rtdb,basis)) call errquit
     &    ('kgdtest_rel2e: error normalizing ',911)
c
      call int_init(rtdb,1,basis)
      memscr = 5 000 000
      membuf = 100 000
      if (.not.ma_push_get(mt_dbl,memscr,' scratch ',
     &    h_scr, k_scr)) call errquit
     &    (' ma error 1',911)
      if (.not.ma_push_get(mt_dbl,membuf,' buf ',
     &    h_buf, k_buf)) call errquit
     &    (' ma error 2',911)
c
      sbas = sc_bsh + BASIS_HANDLE_OFFSET
      abas = ao_bsh + BASIS_HANDLE_OFFSET
c
      difmax = 0.0d0
      errmax = 0.0d0
C      rel_dbg = 3
      do ish = 1,nshell
        ucont = (sf_ibs_cn2ucn(ish,abas))
        Li      = infbs_cont(CONT_TYPE ,ucont,abas)
        ipr  = infbs_cont(CONT_NPRIM,ucont,abas)
        igc   = infbs_cont(CONT_NGEN ,ucont,abas)
        iex  = infbs_cont(CONT_IEXP ,ucont,abas)
        icf  = infbs_cont(CONT_ICFP ,ucont,abas)
        ice  = (sf_ibs_cn2ce(ish,abas))
        igm  = ibs_geom(abas)
        i2 = (Li+1)*(Li+2)/2
        ihi = igc*i2
        ucont = ao_to_ls(ucont)
        icfs = infbs_cont(CONT_ICFP ,ucont,sbas)
c
        do jsh = 1,ish
          ucont = (sf_ibs_cn2ucn(jsh,abas))
          Lj      = infbs_cont(CONT_TYPE ,ucont,abas)
          jpr  = infbs_cont(CONT_NPRIM,ucont,abas)
          jgc   = infbs_cont(CONT_NGEN ,ucont,abas)
          jex  = infbs_cont(CONT_IEXP ,ucont,abas)
          jcf  = infbs_cont(CONT_ICFP ,ucont,abas)
          jce  = (sf_ibs_cn2ce(jsh,abas))
          jgm  = ibs_geom(abas)
          j2 = (Lj+1)*(Lj+2)/2
          jhi = jgc*j2
          ucont = ao_to_ls(ucont)
          jcfs = infbs_cont(CONT_ICFP ,ucont,sbas)
c
          do ksh = 1,ish
            ucont = (sf_ibs_cn2ucn(ksh,abas))
            Lk      = infbs_cont(CONT_TYPE ,ucont,abas)
            kpr  = infbs_cont(CONT_NPRIM,ucont,abas)
            kgc   = infbs_cont(CONT_NGEN ,ucont,abas)
            kex  = infbs_cont(CONT_IEXP ,ucont,abas)
            kcf  = infbs_cont(CONT_ICFP ,ucont,abas)
            kce  = (sf_ibs_cn2ce(ish,abas))
            kgm  = ibs_geom(abas)
            k2 = (Lk+1)*(Lk+2)/2
            khi = kgc*k2
            ucont = ao_to_ls(ucont)
            kcfs = infbs_cont(CONT_ICFP ,ucont,sbas)
c
            do lsh = 1,ksh
              ucont = (sf_ibs_cn2ucn(lsh,abas))
              Ll      = infbs_cont(CONT_TYPE ,ucont,abas)
              lpr  = infbs_cont(CONT_NPRIM,ucont,abas)
              lgc   = infbs_cont(CONT_NGEN ,ucont,abas)
              lex  = infbs_cont(CONT_IEXP ,ucont,abas)
              lcf  = infbs_cont(CONT_ICFP ,ucont,abas)
              lce  = (sf_ibs_cn2ce(lsh,abas))
              lgm  = ibs_geom(abas)
              l2 = (Ll+1)*(Ll+2)/2
              lhi = lgc*l2
              ucont = ao_to_ls(ucont)
              lcfs = infbs_cont(CONT_ICFP ,ucont,sbas)
              
              n_cart = i2*j2*k2*l2
              n_cont = igc*jgc*kgc*lgc
              Nints = ihi*jhi*khi*lhi
              i_eri = k_buf+Nints
C              if (n_cont .gt. 1 .and. Li.ne.Lj .and. Lk.ne.Ll) then
C                write(luout,*)' '
C                write(luout,*)' ============= (',
C     &              ish,':',Li,',',
C     &              jsh,':',Lj,'|',
C     &              ksh,':',Lk,',',
C     &              lsh,':',Ll,')',
C     &              '===================='
C                write(luout,*)' '
C                write (luout,*) 'Nints',Nints

              call hf2(
     &            coords(1,ice,igm),
     &            dbl_mb(mb_exndcf(iex,abas)),
     &            dbl_mb(mb_exndcf(icf,abas)),ipr,igc,Li,
     &            coords(1,jce,jgm),
     &            dbl_mb(mb_exndcf(jex,abas)),
     &            dbl_mb(mb_exndcf(jcf,abas)),jpr,jgc,Lj,
     &            coords(1,kce,kgm),
     &            dbl_mb(mb_exndcf(kex,abas)),
     &            dbl_mb(mb_exndcf(kcf,abas)),kpr,kgc,Lk,
     &            coords(1,lce,lgm),
     &            dbl_mb(mb_exndcf(lex,abas)),
     &            dbl_mb(mb_exndcf(lcf,abas)),lpr,lgc,Ll,
     &            dbl_mb(k_buf),Nints,.false.,.false.,.false.,
     &            .false.,dbl_mb(k_scr),memscr)
C              write (luout,*) 'hf2 completed'
C              call util_flush(6)
              call rel_2e4c_sf (
     &            coords(1,ice,igm),dbl_mb(mb_exndcf(iex,abas)),
     &            dbl_mb(mb_exndcf(icf,abas)),
     &            dbl_mb(mb_exndcf(icfs,sbas)),ipr,igc,Li,ice,
     &            coords(1,jce,jgm),dbl_mb(mb_exndcf(jex,abas)),
     &            dbl_mb(mb_exndcf(jcf,abas)),
     &            dbl_mb(mb_exndcf(jcfs,sbas)),jpr,jgc,Lj,jce,
     &            coords(1,kce,kgm),dbl_mb(mb_exndcf(kex,abas)),
     &            dbl_mb(mb_exndcf(kcf,abas)),
     &            dbl_mb(mb_exndcf(kcfs,sbas)),kpr,kgc,Lk,kce,
     &            coords(1,lce,lgm),dbl_mb(mb_exndcf(lex,abas)),
     &            dbl_mb(mb_exndcf(lcf,abas)),
     &            dbl_mb(mb_exndcf(lcfs,sbas)),lpr,lgc,Ll,lce,
c...................... canAB   canCD   canPQ   DryRun
     &            dbl_mb(i_eri),Nints,.false.,.false.,.false.,.false.,
     &            dbl_mb(k_scr),memscr,
     &            .true.,.true.,ss_one_cent,do_ssss,rel_dbg)

C              write (luout,*) 'rel_2e4c_sf completed'
C              call util_flush(6)
              call daxpy(Nints,-1.0d0,dbl_mb(k_buf),1,dbl_mb(i_eri),1)
              dif = 0.0d0
              do i = 0,Nints-1
                if (dbl_mb(k_buf+i) .ne. 0.0d0) then
                  dif = max(dif,abs(dbl_mb(i_eri+i)/dbl_mb(k_buf+i)))
                else if (dbl_mb(i_eri+i) .ne. 0.0d0) then
                  errmax = max(errmax, abs(dbl_mb(i_eri+i)))
                end if
              end do
              difmax = max(difmax,dif)
              if (dif .gt. 0.1d0) then
                write (LuOut,*) ish,Li,ipr,igc,ice
                write (LuOut,*) jsh,Lj,jpr,jgc,jce
                write (LuOut,*) ksh,Lk,kpr,kgc,kce
                write (LuOut,*) lsh,Ll,lpr,lgc,lce
                write (LuOut,*) dif
              end if
            end do
          end do
        end do
      end do
      write (luout,*) 'Maximum difference of all blocks',difmax
      write (luout,*) 'Maximum difference from zero integrals',errmax

      end
