      logical function ao_1prdm_read(nbf,g_vecs)
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
c
      integer g_vecs
c
      integer nbf               ! No. of functions in basis
      integer ok
      integer l_vecs, k_vecs
      integer i
      integer unitno
      parameter (unitno = 66)
      integer inntsize
c
      l_vecs = -1               ! An invalid MA handle
c
      inntsize = MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file='ao_1prdm',
     $        err=1000)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_read',l_vecs,k_vecs))
     $        call errquit('movecs_read: ma failed', nbf, MA_ERR)
         do i = 1, nbf
            call sread(unitno, dbl_mb(k_vecs), nbf)
            call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), 1)
         enddo
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
         if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('ao_1prdm_read: pop failed', l_vecs, MA_ERR)
      endif
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
c
      ao_1prdm_read = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. ao_1prdm_read .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) 'ao_1prdm'
 22      format(/' Read one-particle reduced density matrix from ',a/)
         call util_flush(luout)
      endif
      return
c
 1000 write(6,*) ' ao_1prdm_read: failed to open ',
     $     'ao_1prdm'
      call util_flush(luout)
      ok = 0
      goto 10
c
 1002 write(6,*) ' ao_1prdm_read: failed to close',
     $     'ao_1prdm'
      call util_flush(luout)
      ok = 0
      goto 10
c
      end

c $Id$
