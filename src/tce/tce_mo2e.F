      subroutine tce_mo2e(rtdb,d_ao2e,d_mo2e,k_2e_offset)
c
c $Id$
c
c     Returns a SF handle for the MO 2-e integral matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb                 ! Run-time database
      logical nodezero             ! True if node zero
      integer d_ao2e               ! Half MO integrals
      integer d_mo2e               ! MO integrals
      integer k_2e_offset          ! MO integrals offset
      integer l_kl,k_kl            ! Half MO integrals
      integer l_sl,k_sl            ! Quarter transformed
      integer l_pqkl_packed,k_pqkl_packed! Packed
      integer pq_packed_off
      integer kl
      integer l_pqrs,k_pqrs        ! MO integrals
      integer pb,qb,rb,sb
      integer p,q
      integer spin_pq,spin_rs
      integer sym_pq,sym_rs
      integer sf_size_kl
      integer sf_size_rs
      integer sf_size_pq
      integer sf_size_pqrs
      integer pqoff
      integer l_pqboff,k_pqboff
      integer pq
      integer r_range,s_range
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      integer next,count,nprocs
      integer Nb
c     integer pqrs,r,s
c     integer l_debug,k_debug
c     double precision debugdbl
      logical transpose
c
c     ===========
c     Preliminary
c     ===========
c
      nodezero = (ga_nodeid().eq.0)
      Nb = noa + nob + nva + nvb
c
c     ===========
c     Allocate SF
c     ===========
c
      sf_size_pqrs = 0
      pq_packed_off = 0
      if (.not.ma_push_get(mt_int,Nb*Nb,'pqboff',l_pqboff,k_pqboff))
     1  call errquit('tce_mo2e: MA problem',11,MA_ERR)
      do pb = 1,Nb
        do qb = pb,Nb
          do rb = 1,Nb
            do sb = rb,Nb
              if (int_mb(k_spin+pb-1)+int_mb(k_spin+qb-1) .eq.
     1            int_mb(k_spin+rb-1)+int_mb(k_spin+sb-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+pb-1)+int_mb(k_spin+qb-1)+
     2           int_mb(k_spin+rb-1)+int_mb(k_spin+sb-1).ne.8)) then
              if (ieor(int_mb(k_sym+pb-1),ieor(int_mb(k_sym+qb-1),
     1            ieor(int_mb(k_sym+rb-1),int_mb(k_sym+sb-1))))
     2            .eq. 0) then
                sf_size_pqrs = sf_size_pqrs
     1            + int_mb(k_range+pb-1) * int_mb(k_range+qb-1)
     2            * int_mb(k_range+rb-1) * int_mb(k_range+sb-1)
              endif
              endif
              endif
            enddo
          enddo
          int_mb(k_pqboff+(pb-1)*Nb+qb-1) = pq_packed_off
          pq_packed_off = pq_packed_off + int_mb(k_range+pb-1)
     1                                  * int_mb(k_range+qb-1)
        enddo
      enddo
      sf_size_kl = nbf**2
c
c     ====================
c     Allocate work arrays
c     ====================
c
      if (.not.ma_push_get(mt_dbl,sf_size_kl,'kl',l_kl,k_kl))
     1  call errquit('tce_mo2e: MA problem',1,MA_ERR)
c
c     ============
c     Loop over pq
c     ============
c
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do pb = 1,Nb
        do qb = pb,Nb
          if (next.eq.count) then
            if (util_print('mo2e',print_debug))
     1        write(LuOut,9190) pb,qb,ga_nodeid()
            spin_pq = int_mb(k_spin+pb-1)+int_mb(k_spin+qb-1)
            sym_pq = ieor(int_mb(k_sym+pb-1),int_mb(k_sym+qb-1))
            sf_size_pq = int_mb(k_range+pb-1)*int_mb(k_range+qb-1)
            if (.not.ma_push_get(mt_dbl,sf_size_pq*sf_size_kl,
     1        'pqkl packed',l_pqkl_packed,k_pqkl_packed))
     2        call errquit('tce_mo2e: MA problem',2,MA_ERR)
            do kl = 1,sf_size_kl
              call get_block(d_ao2e,dbl_mb(k_pqkl_packed
     1          +(kl-1)*sf_size_pq),sf_size_pq,(kl-1)*pq_packed_off
     2          +int_mb(k_pqboff+(pb-1)*Nb+qb-1))
            enddo
            do rb = 1,Nb
              do sb = rb,Nb
                transpose = (int_mb(k_spin+pb-1).eq.int_mb(k_spin+rb-1))
                spin_rs = int_mb(k_spin+rb-1)+int_mb(k_spin+sb-1)
                sym_rs = ieor(int_mb(k_sym+rb-1),int_mb(k_sym+sb-1))
                if (spin_pq .eq. spin_rs) then
                if ((.not.restricted).or.(spin_pq+spin_rs.ne.8)) then
                if (ieor(sym_pq,sym_rs) .eq. 0) then
                  r_range = int_mb(k_range+rb-1)
                  s_range = int_mb(k_range+sb-1)
                  sf_size_rs = r_range * s_range
                  if (.not.ma_push_get(mt_dbl,s_range*nbf,
     1              'sl',l_sl,k_sl))
     2              call errquit('tce_mo2e: MA problem',3,MA_ERR)
                  if (.not.ma_push_get(mt_dbl,sf_size_pq*sf_size_rs,
     1              'rs',l_pqrs,k_pqrs))
     2              call errquit('tce_mo2e: MA problem',4,MA_ERR)
                  pq = 0
                  pqoff = 0
                  do p = 1,int_mb(k_range+pb-1)
                    do q = 1,int_mb(k_range+qb-1)
                      pq = pq + 1
                      call dcopy(sf_size_kl,
     1                  dbl_mb(k_pqkl_packed+pq-1),
     2                  sf_size_pq,dbl_mb(k_kl),1)
                      if (transpose) then
#ifdef MA_TRANS_BLOCKED
                        call ma_transpose_blocked(dbl_mb(k_kl),nbf)
#else
                        call ma_transpose(dbl_mb(k_kl),nbf)
#endif
                      else
                        call ma_scale(dbl_mb(k_kl),nbf**2,-1.0d0)
                      endif
                      call dgemm('T','N',s_range,nbf,nbf,
     1                  1.0d0,dbl_mb(k_movecs_sorted
     2                  +int_mb(k_offset+sb-1)*nbf),nbf,
     3                  dbl_mb(k_kl),nbf,
     4                  0.0d0,dbl_mb(k_sl),s_range)
                      call dgemm('N','N',s_range,r_range,nbf,
     1                  1.0d0,dbl_mb(k_sl),s_range,
     2                  dbl_mb(k_movecs_sorted
     3                  +int_mb(k_offset+rb-1)*nbf),nbf,
     4                  0.0d0,dbl_mb(k_pqrs+pqoff),s_range)
                      pqoff = pqoff + r_range * s_range
                    enddo
                  enddo
                  if (util_print('mo2e',print_debug))
     1              then
                    write(LuOut,9200) pb,qb,rb,sb
                    call ma_print_compact(dbl_mb(k_pqrs),
     1              sf_size_pq*sf_size_rs,1,'MO 2e')
                  endif
                  call put_hash_block(d_mo2e,dbl_mb(k_pqrs),
     1              sf_size_pq*sf_size_rs,int_mb(k_2e_offset),
     2              ((((pb-1)*Nb+qb-1)*Nb+rb-1)*Nb+sb-1))
                  if (.not.ma_pop_stack(l_pqrs))
     1              call errquit('tce_mo2e: MA problem',5,MA_ERR)
                  if (.not.ma_pop_stack(l_sl))
     1              call errquit('tce_mo2e: MA problem',6,MA_ERR)
                endif
                endif
                endif
              enddo
            enddo
            if (.not.ma_pop_stack(l_pqkl_packed))
     1        call errquit('tce_mo2e: MA problem',3,MA_ERR)
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ma_pop_stack(l_kl))
     1  call errquit('tce_mo2e: MA problem',10,MA_ERR)
      if (.not.ma_pop_stack(l_pqboff))
     1  call errquit('tce_mo2e: MA problem',12,MA_ERR)
      if (nodezero.and.util_print('time',print_debug)) then
        write(LuOut,8000) 'createfile   ', cpusecs(1), cpusecs(51)
        write(LuOut,8000) 'deletefile   ', cpusecs(2), cpusecs(52)
        write(LuOut,8000) 'put_block    ', cpusecs(3), cpusecs(53)
        write(LuOut,8000) 'add_block    ', cpusecs(4), cpusecs(54)
        write(LuOut,8000) 'get_block    ', cpusecs(5), cpusecs(55)
        write(LuOut,8000) 'reconcilefile', cpusecs(6), cpusecs(56)
      endif
 8000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
c
c     Debug only
c
c     if (.not.ma_push_get(mt_dbl,sf_size_pqrs,
c    1  'mo2e',l_debug,k_debug))
c    2  call errquit('tce_mo2e: MA problem',-1,MA_ERR)
c     call get_block(d_mo2e,dbl_mb(k_debug),sf_size_pqrs,0)
c     pqrs = 0
c     do pb = 1,Nb
c       do qb = pb,Nb
c         do rb = 1,Nb
c           do sb = rb,Nb
c             if (int_mb(k_spin+pb-1)+int_mb(k_spin+qb-1) .eq.
c    1            int_mb(k_spin+rb-1)+int_mb(k_spin+sb-1)) then
c               if (ieor(int_mb(k_sym+pb-1),ieor(int_mb(k_sym+qb-1),
c    1              ieor(int_mb(k_sym+rb-1),int_mb(k_sym+sb-1))))
c    2              .eq. 0) then
c                 do p = 1,int_mb(k_range+pb-1)
c                   do q = 1,int_mb(k_range+qb-1)
c                     do r = 1,int_mb(k_range+rb-1)
c                       do s = 1,int_mb(k_range+sb-1)
c                         pqrs = pqrs + 1
c                         debugdbl = dbl_mb(k_debug+pqrs-1)
c                         write(LuOut,9100) p+int_mb(k_offset+pb-1),
c    1                                      q+int_mb(k_offset+qb-1),
c    2                                      r+int_mb(k_offset+rb-1),
c    3                                      s+int_mb(k_offset+sb-1),
c    4                                      debugdbl
c                       enddo
c                     enddo
c                   enddo
c                 enddo
c               endif
c             endif
c           enddo
c         enddo
c       enddo
c     enddo
c     if (.not.ma_pop_stack(l_debug))
c    1  call errquit('tce_mo2e: MA problem',7,MA_ERR)
c9100 format(4i3,f20.10)
c
c     ======
c     Return
c     ======
c
 9190 format(/,' p q blocks ',2i3,' has been assigned to node ',i3)
 9200 format(/,' p q r s blocks ',4i3)
      call deletefile(d_ao2e)
      return
      end
