c $Id$
c:tex-% part of the internal API routines.
c:tex-\subsection{int\_rel\_init}
c:tex-This routine initializes the relativistic integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function int_rel_init(rtdb,nbas,basis)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "global.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "rel_nwc.fh"
#include "rtdb.fh"
c::passed
      integer rtdb, nbas, basis
      dimension basis(nbas)
c::local
      integer geom
      integer i,lbas,sbas,nucontl,nuconts,iucontl,iuconts,iexpl,
     &    iexps,type,nprim,ngen,bastag
      logical odum, oauto
      logical any_large, any_small, auto_large, auto_small
      double precision xl,xs
c
      logical int_normalize
      external int_normalize
      logical  bas_name_exist_rtdb
      external bas_name_exist_rtdb
#include "bas_exndcf.fh"
c
      if (init_rel) call errquit('int_rel_init: already called',911,
     &       INT_ERR)
      init_rel = .true.
c
      if (.not.bas_get_ao_handle (basis,nbas,ao_bsh))
     &    call errquit('int_rel_init: could not find ao basis',911,
     &       BASIS_ERR)
      if (.not.bas_any_rel_shells (ao_bsh))
     &    call errquit('int_rel_init: ao basis has no rel shells',99,
     &       BASIS_ERR)
      if (.not. bas_geom(ao_bsh, geom)) 
     &    call errquit('int_rel_init: could not get ao basis geometry',
     &       BASIS_ERR,
     &    99)
      if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
     &    call errquit('int_rel_init: rtdb_get failed for atomic_nesc',
     &       RTDB_ERR,
     &    99)
      if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
     &    call errquit('int_rel_init: ' // '
     &    rtdb_get failed for nesc_1e_approx',99, RTDB_ERR)
c
      any_large = bas_name_exist_rtdb(rtdb,large_cmpt_name)
      any_small = bas_name_exist_rtdb(rtdb,small_cmpt_name)
      auto_large = bas_name_exist_rtdb(rtdb,auto_large_cmpt_name)
      auto_small = bas_name_exist_rtdb(rtdb,auto_small_cmpt_name)
      oauto = .false.
c
c   get small component basis set
c
      if (.not. bas_create(sc_bsh, 'small component'))
     &    call errquit('int_rel_init: ' // 
     &    'could not create basis set for small component', 99,
     &       BASIS_ERR)
      if (any_small) then
        if (.not. bas_rtdb_load(rtdb, geom, sc_bsh, small_cmpt_name))
     &      call errquit('int_rel_init: ' //
     &      'small component exists but could not load it',99,
     &       UNKNOWN_ERR)
      else if (auto_small) then
        oauto = bas_rtdb_load(rtdb, geom, sc_bsh, auto_small_cmpt_name)
        if (.not.oauto) call errquit('int_rel_init: ' //
     &      'could not load auto small component',99,
     &       UNKNOWN_ERR)
      else
        call errquit ('int_rel_init: ' //
     &      'small component does not exist!',99,
     &       UNKNOWN_ERR)
      end if
c 
c   Set up mapping array between ao basis and small component basis
c
      lbas = ao_bsh + BASIS_HANDLE_OFFSET
      nucontl = infbs_head(HEAD_NCONT,lbas)
    3 sbas = sc_bsh + BASIS_HANDLE_OFFSET
      nuconts = infbs_head(HEAD_NCONT,sbas)
      do iuconts = 1,nuconts
        ls_to_ao(iuconts) = 0
      end do
      do 11 iucontl = 1, nucontl
        ao_to_ls(iucontl) = 0
        if (infbs_cont(CONT_RELLS,iucontl,lbas) .ne. 0) then
          type  = infbs_cont(CONT_TYPE, iucontl,lbas)
          if (type .lt. 0) call errquit ('int_rel_init: ' //
     &        'sp/spd shells illegal for relativistic basis',99,
     &       BASIS_ERR)
          nprim = infbs_cont(CONT_NPRIM,iucontl,lbas) 
          ngen  = infbs_cont(CONT_NGEN,iucontl,lbas)
          iexpl  = infbs_cont(CONT_IEXP,iucontl,lbas)
          bastag  = infbs_cont(CONT_TAG,iucontl,lbas)
          do 10 iuconts = 1,nuconts
            if (ls_to_ao(iuconts) .ne. 0) go to 10
            if ((type.ne.infbs_cont(CONT_TYPE,iuconts,sbas)) .or.
     &          (nprim.ne.infbs_cont(CONT_NPRIM,iuconts,sbas)) .or.
     &          (ngen.ne.infbs_cont(CONT_NGEN, iuconts,sbas)) .or.
     &          (bastag.ne.infbs_cont(CONT_TAG, iuconts,sbas))) go to 10
            iexps  = infbs_cont(CONT_IEXP,iuconts,sbas)
            do i = 0,nprim-1
              xl = dbl_mb(mb_exndcf(iexpl+i,lbas))
              xs = dbl_mb(mb_exndcf(iexps+i,sbas))
              if (abs(xl-xs) .gt. abs(min(xl,xs))*1.0e-10) go to 10
            end do
            ao_to_ls(iucontl) = iuconts
            ls_to_ao(iuconts) = iucontl
            go to 11
   10     continue
c
c     Mismatch between ao and small component: shell from ao not found
c     in sc. Check if an auto-small basis has been loaded; if not, load it 
c     and try to match the ao and auto-small bases.
c
          if (ao_to_ls(iucontl) .eq. 0) then
            if (auto_small .and. .not.oauto) then
              oauto = bas_rtdb_load(rtdb, geom, sc_bsh, 
     &            auto_small_cmpt_name)
              if (.not.oauto) call errquit('int_rel_init: ' //
     &            'could not load auto small component',99,
     &       UNKNOWN_ERR)
              go to 3
            else
              write (luout,*) 'Error for ao basis unique contraction',
     &            iucontl
              call errquit ('int_rel_init: ' //
     &            'Could not match ao basis shell to small component',
     &            99, BASIS_ERR)
            end if
          end if
        end if
   11 continue
      odum = int_normalize(rtdb,sc_bsh)
c
c     get large component basis set
c
      if (nesc_1e_approx) then
        if (any_large .or. auto_large) then
          if (.not. bas_create(lc_bsh, large_cmpt_name))
     &        call errquit('int_rel_init: ' // 
     &        'could not create basis set for large component', 99,
     &       UNKNOWN_ERR)
          if (oauto) then
            if (.not. auto_large) call errquit('int_rel_init: ' //
     &          'could not find auto large component basis ',99,
     &       BASIS_ERR)
            if (bas_rtdb_load(rtdb, geom, lc_bsh, auto_large_cmpt_name))
     &          go to 12
          else
            if (.not. any_large) call errquit('int_rel_init: ' //
     &          'could not find large component basis ',99,
     &       BASIS_ERR)
            if (bas_rtdb_load(rtdb, geom, lc_bsh, large_cmpt_name)) 
     &          go to 12
          end if
          call errquit(
     &        'int_rel_init: could not load large component basis',99,
     &       BASIS_ERR)
        else
          call errquit('int_rel_init: ' //
     &        'large component basis expected but does not exist',99,
     &       BASIS_ERR)
        end if
c
c       check that large and small component basis sets match
c
   12   lbas = lc_bsh + BASIS_HANDLE_OFFSET
        nucontl = infbs_head(HEAD_NCONT,lbas)
        if (nucontl.ne.nuconts) call errquit('int_rel_init: ' //
     &      'Large and small component bases do not match',99,
     &       BASIS_ERR)
        do iucontl = 1, nucontl
          type  = infbs_cont(CONT_TYPE, iucontl,lbas)
          nprim = infbs_cont(CONT_NPRIM,iucontl,lbas) 
          ngen  = infbs_cont(CONT_NGEN, iucontl,lbas)
          iexpl  = infbs_cont(CONT_IEXP,iucontl,lbas)
          bastag  = infbs_cont(CONT_TAG, iucontl,lbas)
          if ((type.ne.infbs_cont(CONT_TYPE,iucontl,sbas)) .or.
     &        (nprim.ne.infbs_cont(CONT_NPRIM,iucontl,sbas)) .or.
     &        (ngen.ne.infbs_cont(CONT_NGEN, iucontl,sbas)) .or.
     &        (bastag.ne.infbs_cont(CONT_TAG, iucontl,sbas))) 
     &        call errquit ('int_rel_init: ' //
     &        'Large and small component bases do not match',99,
     &       BASIS_ERR)
          iexps  = infbs_cont(CONT_IEXP,iucontl,sbas)
          do i = 1,nprim
            xl = dbl_mb(mb_exndcf(iexpl,lbas))
            xs = dbl_mb(mb_exndcf(iexps,sbas))
            if (abs(xl-xs) .gt. abs(min(xl,xs))*1.0e-10) 
     &          call errquit ('int_rel_init: ' //
     &          'Large and small component bases do not match',99,
     &       BASIS_ERR)
          end do
        end do
        odum = int_normalize(rtdb,lc_bsh)
      end if
c
      level_SSSS = 0
      rel_dbg = 0
      rel_typ = 1
c
      int_rel_init = .true.
c
      end
*
      subroutine int_rel_term
      implicit none
#include "nwc_const.fh"
#include "bas.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "errquit.fh"
c
      if (nesc_1e_approx) then
        if (lc_bsh .ne. -1) then
          if (.not. bas_destroy(lc_bsh)) call errquit(
     &        'int_rel_terminate: failed to destroy' // 
     &        'large component basis',0, BASIS_ERR)
          lc_bsh = -1
        end if
      end if
      if (sc_bsh .ne. -1) then
        if (.not. bas_destroy(sc_bsh)) call errquit(
     &      'int_rel_terminate: failed to destroy small component basis'
     &      ,0, BASIS_ERR)
        sc_bsh = -1
      endif
      ao_bsh = -1
      init_rel = .false.
c
      end

      block data rel_init_bd
c
c   Block data for relativistic integrals
c
#include "nwc_const.fh"
#include "rel_nwc.fh"
      data init_rel /.false./
      data small_cmpt_name/'small component'/
      data large_cmpt_name/'large component'/
      data auto_small_cmpt_name/'auto small component'/
      data auto_large_cmpt_name/'auto large component'/
      end
