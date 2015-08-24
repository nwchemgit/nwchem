      logical function ao_1prdm_write(nbf,g_vecs)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
c     Temporary routine
c
      integer nbf               ! [input] No. of functions in basis
      integer g_vecs     ! Global array with eigen-vectors
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 66)   ! These need to be managed !!!
      integer l_vecs, k_vecs
      integer i
      integer ok
      integer inntsize
c
      inntsize = MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
      ok = 0
c
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='unknown', form='unformatted',
     $        file='ao_1prdm', err=1000)
         if (.not. ma_push_get(mt_dbl,nbf,'ao_1prdm_write',
     1      l_vecs,k_vecs))
     $        call errquit('ao_1prdm_write: ma failed', nbf, MA_ERR)
         do i = 1, nbf
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs),1)
            call swrite(unitno, dbl_mb(k_vecs), nbf)
         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('ao_1prdm_write: ma pop failed', l_vecs,
     &       MA_ERR)
         close(unitno,err=1002)
         ok = 1
      endif
c
   10 call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)
c
      ao_1prdm_write = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) 'ao_1prdm'
   22    format(/' Wrote one-particle reduced density matrix to ',a/)
         call util_flush(luout)
      endif
      call ga_sync()
      return
c
 1000 write(6,*) ' ao_1prdm_write: failed to open ',
     $     'ao_1prdm'
      call util_flush(luout)
      ok = 0
      goto 10
c
 1002 write(6,*) ' ao_1prdm_write: failed to close',
     $     'ao_1prdm'
      call util_flush(luout)
      ok = 0
      goto 10
c
      end

c $Id$
