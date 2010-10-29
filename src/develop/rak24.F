*
* $Id$
*
      logical function raktask_ecppe(rtdb)
*
* compute pe matrix with given basis set
*
      implicit none
#include "errquit.fh"
*::includes
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
*::functions
      logical int_normalize
      external int_normalize
*::passed
      integer rtdb
*::local      
      integer geom, basis, nbf, nshells, max1e, mscr, h_buff, k_buff
      integer h_scr, k_scr, ish, jsh, ilo, ihi, jlo, jhi
      integer nbfi, nbfj, ish_start, jsh_start
      raktask_ecppe = .false.

* create geometry handle
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('task_ecppe: geom_create failed?',911,
     &       GEOM_ERR)
* load geometry from rtdb
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &      call errquit
     &      ('task_ecppe: geom_rtdb_load failed?',911, RTDB_ERR)
* create basis set handle
      if (.not.bas_create(basis,'ao basis'))
     &      call errquit('task_ecppe: bas_create failed?',911,
     &       BASIS_ERR)
* load basis set
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &      call errquit
     &      ('task_ecppe: bas_rtdb_load failed?',911, RTDB_ERR)
c
* query basis set object for the number of functions and shells      
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &      ('task_ecppe: bas_numbf failed?',911, BASIS_ERR)
      if (.not.bas_numcont(basis,nshells)) call errquit
     &      ('task_ecppe: bas_numcont failed?',911, BASIS_ERR)
*
*normalize basis set
      if (.not.int_normalize(rtdb,basis)) call errquit
     &      ('task_ecppe: int_normalize failed?',911, INT_ERR)
      if (.not.bas_print(basis)) stop ' err bas_print'
      if (.not.gbs_map_print(basis)) stop ' err bas_print'
*
* initialize the integral API
      call int_init(rtdb,1,basis)
*
* query integral API for memory requirements.
      call int_mem_1e(max1e, mscr)
*
* allocate local memory buffers for integrals and scratch space for 
* the integral API
      if (.not. ma_push_get(mt_dbl,max1e,'integral buffer',
     &      h_buff,k_buff)) call errquit
     &      ('task_ecppe: could not allocat integral buffer',911,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,mscr,'integral scratch',
     &      h_scr,k_scr)) call errquit
     &      ('task_ecppe: could not allocat integral buffer',911,
     &       MA_ERR)
      if (.not.rtdb_get(rtdb,'rak24:ish_start',mt_int,1,ish_start))
     &    ish_start = 1
      if (.not.rtdb_get(rtdb,'rak24:jsh_start',mt_int,1,jsh_start))
     &    jsh_start = 1
      do ish = ish_start,nshells
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) call errquit
     &      ('task_ecppe: bas_cn2bfr failed',911, BASIS_ERR)
        nbfi = ihi - ilo + 1
        do jsh = jsh_start,nshells
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     &        ('task_ecppe: bas_cn2bfr failed',911, BASIS_ERR)
          nbfj = jhi - jlo + 1
*...compute it: pe(ilo:ihi,jlo:jhi)
          call util_flush(6)
          call dfill(mscr,0.0d00,dbl_mb(k_scr),1)
          call dfill(max1e,0.0d00,dbl_mb(k_buff),1)
          write(6,*)'<ish,jsh> = <',ish,',',jsh,'>...'
          call util_flush(6)
          call int_1epe(basis,jsh,basis,ish,
     &        mscr,dbl_mb(k_scr),
     &        max1e,dbl_mb(k_buff))
          write(6,*)'<ish,jsh> = <',ish,',',jsh,'>...done'
          call util_flush(6)
        enddo
      enddo
      
*
      raktask_ecppe = .true.
*
      end
