c $Id: int_rel_init.F,v 1.1 1999-10-26 19:14:27 mg201 Exp $
c:tex-% part of the internal API routines.
c:tex-\subsection{int\_rel\_init}
c:tex-This routine initializes the relativistic integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function int_rel_init(rtdb,nbas,basis)
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "apiP.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "bas.fh"
#include "rel_nwc.fh"
c::passed
      integer rtdb, nbas, basis
      dimension basis(nbas)
c::local
      integer geom
      integer ibas,lbas,sbas,nucont,iucont,type,nprim,ngen,bastag
      logical odum
      character*255 basis_name,trans_name
c
      logical int_normalize
      external int_normalize
      logical  bas_name_exist_rtdb
      external bas_name_exist_rtdb
c
      any_large = bas_name_exist_rtdb(rtdb,large_cmpt_name)
      any_small = bas_name_exist_rtdb(rtdb,small_cmpt_name)
      if (init_rel) then
        write(6,*)' already called int_rel_init' 
        call errquit('int_rel_init error',911)
      end if
      init_rel = .true.
      do ibas = 1,nbas
        odum = bas_name(basis(ibas),basis_name,trans_name)
        if (basis_name(1:8) .eq. 'ao basis') go to 1
      end do
      call errquit('int_rel_init: could not find ao basis',911)
    1 if (.not. bas_geom(basis(ibas), geom)) 
     &    call errquit('int_rel_init: bas_geom?',0)
c
c   get small component basis set
c
      if (.not. bas_create(sc_bsh, 'small component'))
     &    call errquit('int_rel_init: bas_create sc_bsh?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, sc_bsh, small_cmpt_name))
     &    call errquit('int_rel_init: could not load small component',
     &    0)
      odum = int_normalize(sc_bsh)
c
c     check that small component basis matches beginning of ao basis
c
      lbas = basis(ibas) + BASIS_HANDLE_OFFSET
      sbas = sc_bsh + BASIS_HANDLE_OFFSET
      nucont = infbs_head(HEAD_NCONT,sbas)
      do iucont = 1, nucont
        type  = infbs_cont(CONT_TYPE, iucont,lbas)
        if (type .lt. 0) call errquit (
     &      'sp/spd shells illegal for relativistic basis',911)
        nprim = infbs_cont(CONT_NPRIM,iucont,lbas) 
        ngen  = infbs_cont(CONT_NGEN, iucont,lbas)
        bastag  = infbs_cont(CONT_TAG, iucont,lbas)
        if ((type.ne.infbs_cont(CONT_TYPE,iucont,sbas)) .or.
     &      (nprim.ne.infbs_cont(CONT_NPRIM,iucont,sbas)) .or.
     &      (ngen.ne.infbs_cont(CONT_NGEN, iucont,sbas)) .or.
     &      (bastag.ne.infbs_cont(CONT_TAG, iucont,sbas))) 
     &      call errquit (
     &      'Small component basis does not match ao basis',911) 
      end do
c
c     get large component basis set
c
      if (any_large) then
        if (.not. bas_create(lc_bsh, large_cmpt_name))
     $      call errquit('int_rel_init: bas_create lc_bsh?', 0)
        if (.not.bas_rtdb_load(rtdb, geom, lc_bsh, large_cmpt_name))
     &      call errquit(
     &      'int_rel_init: could not load large component basis',0)
        odum = int_normalize(lc_bsh)
c
c       check that large and small component basis sets match
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        nucont = infbs_head(HEAD_NCONT,lbas)
        if (nucont.ne.infbs_head(HEAD_NCONT,sbas)) call errquit
     &      ('Large and small component bases do not match',911)
        do iucont = 1, nucont
          type  = infbs_cont(CONT_TYPE, iucont,lbas)
          nprim = infbs_cont(CONT_NPRIM,iucont,lbas) 
          ngen  = infbs_cont(CONT_NGEN, iucont,lbas)
          bastag  = infbs_cont(CONT_TAG, iucont,lbas)
          if ((type.ne.infbs_cont(CONT_TYPE,iucont,sbas)) .or.
     &        (nprim.ne.infbs_cont(CONT_NPRIM,iucont,sbas)) .or.
     &        (ngen.ne.infbs_cont(CONT_NGEN, iucont,sbas)) .or.
     &        (bastag.ne.infbs_cont(CONT_TAG, iucont,sbas))) 
     &        call errquit (
     &        'Large and small component bases do not match',911) 
        end do
c
c  this next line is only for the one-electron approximation. It must
c  be removed when the two-electron code is added.
c
      else
        call errquit('Must have both large and small component',911)
      end if
c
      call rel_init_con
      int_rel_init = .true.
c
      end
*
      subroutine int_rel_term
      implicit none
#include "bas.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
c
      if (lc_bsh .ne. -1) then
        if (.not. bas_destroy(lc_bsh)) call errquit(
     &      'int_rel_terminate: failed to destroy large component basis'
     &      ,0)
        lc_bsh = -1
      end if
      if (sc_bsh .ne. -1) then
        if (.not. bas_destroy(sc_bsh)) call errquit(
     &      'int_rel_terminate: failed to destroy small component basis'
     &      ,0)
        sc_bsh = -1
      endif
      any_rel = .false.
      any_large = .false.
      any_small = .false.
      init_rel = .false.
c
      end

      block data rel_init_bd
c
c   Block data for relativistic integrals
c
#include "rel_nwc.fh"
      data init_rel /.false./
      data small_cmpt_name/'small component'/
      data large_cmpt_name/'large component'/
      end

