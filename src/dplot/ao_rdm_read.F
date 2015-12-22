      logical function ao_rdm_read(filename,basis, g_dens)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "bas.fh"
c
c
c      integer basis
c      integer g_dens(3)
c
c      integer nbf               ! No. of functions in basis
c      integer l_ao_rdm, k_ao_rdm
c      integer unitno,i
c      parameter (unitno = 67)
c
c
c      if (.not. bas_numbf(basis,nbf)) call
c     &    errquit('ao_rdm_read: could not get nbf',0, BASIS_ERR)
c      call ga_sync()
c      if (ga_nodeid() .eq. 0) then
c         open(unitno, status='old', form='unformatted', file='ao_rdm',
c     $        err=1000)
cc
c         if (.not. ma_push_get(mt_dbl,nbf,'ao_rdm_read',
c     1       l_ao_rdm,k_ao_rdm))
c     2       call errquit('ao_rdm_read: ma failed', nbf, MA_ERR)
c         do i = 1, nbf
c            call sread(unitno, dbl_mb(k_ao_rdm), nbf)
c            call ga_put(g_dens(1), 1, nbf, i, i, dbl_mb(k_ao_rdm), 1)
c         enddo
c        close(unitno,err=1002)
c      endif
c      call ga_sync()
c 1000 write(6,*) ' ao_rdm_read: failed to open '
c 1002 write(6,*) 'ao_rdm_read: failed to close'
c
c      end


c
      character*(*) filename
      integer basis
      integer g_dens(3)
c
      integer nbf
      integer l_ao_rdm,k_ao_rdm,ok
      integer unitno,i
      parameter (unitno = 67)

c      integer iset              ! No. (1,2) of set of vectors to read
c      double precision occ(*)   ! Must be at least nbf long (not nmo)
c      double precision evals(*) ! Must be at least nbf long (not nmo)
c      integer g_vecs
c
c      integer nsets             ! No. of sets of vectors
c      integer nbf               ! No. of functions in basis
c      integer nmo(2)            ! No. of vectors in each set
c      integer ok, jset, i, j
c      integer l_vecs, k_vecs
c      integer unitno
c      parameter (unitno = 67)
      integer inntsize,ddblsize
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
      if (.not.bas_numbf(basis,nbf)) call
     &    errquit('ao_rdm_read: could not get nbf',0, BASIS_ERR)
         if (.not. ma_push_get(mt_dbl,nbf,'ao_rdm_read',
     1       l_ao_rdm,k_ao_rdm))
     2   call errquit('ao_rdm_read: ma problem',0,ma_err)
         do i = 1, nbf
            call sread(unitno, dbl_mb(k_ao_rdm), nbf)
            call ga_put(g_dens, 1, nbf, i, i, dbl_mb(k_ao_rdm), nbf)
         enddo
         if(.not.ma_pop_stack(l_ao_rdm))
     1   call errquit('ao_rdm_read: pop failed',l_ao_rdm,ma_err)
c
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
c      if (l_vecs .ne. -1) then
c	 if (.not. ma_pop_stack(l_vecs)) call errquit
c     $      ('movecs_read: pop failed', l_vecs, MA_ERR)
c      endif
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
c      if (ok .eq. 1) then
c	 call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, inntsize, 0)
c	 call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, ddblsize*nbf, 0)
c	 call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ, ddblsize*nbf, 0)
c      endif
c
      ao_rdm_read = ok .eq. 1
c      if (ga_nodeid() .eq. 0 .and. ao_rdm_read .and.
c     $     util_print('vectors i/o', print_high)) then
c         write(6,22) filename(1:inp_strlen(filename))
c 22      format(/' Read ao 1-particle density matrix from ',a/)
c         call util_flush(luout)
c      endif
c      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' ao_rdm_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      return
      goto 10
c
 1001 write(6,*) ' ao_rdm_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      return
      goto 10
c
 1002 write(6,*) ' ao_rdm_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      return
      goto 10
c
      end
c $Id$
