      SUBROUTINE ccsd_t_restart(rtdb,d_t1,k_t1_offset,
     1                               d_t2,k_t2_offset,
     2                               d_v2,k_v2_offset,
     3                               t_energy)
C
C     $Id$
C
      IMPLICIT NONE
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      integer rtdb,node
      integer d_t1
      integer k_t1_offset
      integer d_t2
      integer k_t2_offset
      integer d_v2
      integer k_v2_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_singles,l_singles
      integer k_doubles,l_doubles
      integer size,i
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer size_table
      integer ccsd_t_restart_begin
      integer outer_virtual_index
      integer k_triples_table,l_triples_table
      double precision energy,t_energy
      double precision factor
#ifdef DEBUG_PRINT
      double precision energy_debug
#endif
      external nxtask
      logical nodezero
c
      node = ga_nodeid()
      nodezero=(node.eq.0)
c
c     Restarting for the occupied triplets
c
      size_table = nvab
c
c     Initialize the starting index
c
      if (.not.rtdb_get(rtdb,'tce:ccsd_t_restart_begin',mt_int,1,
     1                  ccsd_t_restart_begin)) then
        ccsd_t_restart_begin = 1
        if (.not.rtdb_put(rtdb,'tce:ccsd_t_restart_begin',mt_dbl,1,
     1                    ccsd_t_restart_begin)) then
          call errquit('ccsd_t_restart: failed to write to RTDB',0,
     1                 RTDB_ERR)
        endif
      endif
#ifdef DEBUG_PRINT
      write(6,*) 'size_table = ',size_table
      write(6,*) 'ccsd_t_restart_begin = ',ccsd_t_restart_begin
      call util_flush(6)
#endif
c
c     Allocate the triples table
c
      if (.not.ma_push_get(mt_dbl,size_table,'(T) restart table',
     1                     l_triples_table,k_triples_table)) 
     2    call errquit('(T) restart table',size_table,MA_ERR) 
c
c     Initialize the triples table
c
      if (.not.rtdb_get(rtdb,'tce:restart_triples_table',mt_dbl,
     1    size_table,dbl_mb(k_triples_table))) then
        do i = 1,size_table
          dbl_mb(k_triples_table+i-1) = 0.0d0
        enddo
        if (.not.rtdb_put(rtdb,'tce:restart_triples_table',mt_dbl,
     1      size_table,dbl_mb(k_triples_table))) then
            call errquit('ccsd_t_restart: triples table RTDB put',
     1                   0,RTDB_ERR)
        endif
      endif
#ifdef DEBUG_PRINT
      write(6,9000)
      do i = 1,size_table
        write(6,9001) node,i,dbl_mb(k_triples_table+i-1)
      enddo
      write(6,9000)
      call util_flush(6)
#endif
c
c     Evaluating the (T) energy
c
c      nprocs = GA_NNODES()
c      count = 0
c      next = nxtask(nprocs,1)
#ifdef DEBUG_PRINT
      energy_debug = 0.0d0
#endif
      do t_p4b = noab+ccsd_t_restart_begin,noab+nvab
       energy=0.0d0
       outer_virtual_index = t_p4b - noab

       nprocs = GA_NNODES()
       count = 0
       next = nxtask(nprocs,1)
c
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
            if (next.eq.count) then                        
            if (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3      .eq.int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1)) then
            if ((.not.restricted).or.
     1         (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3         +int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1).le.8)) then
            if (ieor(int_mb(k_sym+t_p4b-1),
     1          ieor(int_mb(k_sym+t_p5b-1),
     2          ieor(int_mb(k_sym+t_p6b-1),
     3          ieor(int_mb(k_sym+t_h1b-1),
     4          ieor(int_mb(k_sym+t_h2b-1),
     5               int_mb(k_sym+t_h3b-1)))))).eq.0) then
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
            do i = 1, size
             dbl_mb(k_singles+i-1) = 0.0d0
             dbl_mb(k_doubles+i-1) = 0.0d0
            enddo
            call ccsd_t_singles(dbl_mb(k_singles),d_t1,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            if (restricted) then
              factor = 2.0d0
            else
              factor = 1.0d0
            endif
            if ((t_p4b.eq.t_p5b).and.(t_p5b.eq.t_p6b)) then
              factor = factor / 6.0d0
            else if ((t_p4b.eq.t_p5b).or.(t_p5b.eq.t_p6b)) then
              factor = factor / 2.0d0
            endif
            if ((t_h1b.eq.t_h2b).and.(t_h2b.eq.t_h3b)) then
              factor = factor / 6.0d0
            else if ((t_h1b.eq.t_h2b).or.(t_h2b.eq.t_h3b)) then
              factor = factor / 2.0d0
            endif
            i = 0
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
                  energy = energy + factor * dbl_mb(k_doubles+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
            endif
            count = count + 1
           enddo
          enddo
         enddo
        enddo
       enddo
c
c      Make sure all processess have finished these loops
c
       call ga_sync()
       next = nxtask(-nprocs,1)
c
c      Write new starting index to RTDB
c
       ccsd_t_restart_begin = outer_virtual_index+1
       if (.not.rtdb_put(rtdb,'tce:ccsd_t_restart_begin',mt_int,1,
     1                   ccsd_t_restart_begin))
     2   call errquit('ccsd_t_restart: failed writing to RTDB',1,
     3   RTDB_ERR)
c
c      Global reduction on (T) energy
c
#ifdef DEBUG_PRINT
       if (nodezero) then
         write(6,*) '====================================='
         write(6,*) 'outer_virtual_index = ',outer_virtual_index
         write(6,*) 'energy before reduction',energy
       endif
#endif
       call ga_dgop(mt_dbl,energy,1,'+')
#ifdef DEBUG_PRINT
       energy_debug = energy_debug + energy
       if (nodezero) then
         write(6,*) 'energy after reduction',energy
         write(6,*) 'energy_debug',energy_debug
         write(6,*) '====================================='
       endif
       call util_flush(6)
#endif
c
c      Update the triples table in memory
c
       dbl_mb(k_triples_table+outer_virtual_index-1) = energy
c
c      Write updated triples table to RTDB
c
       if (.not.rtdb_put(rtdb,'tce:restart_triples_table',mt_dbl,
     1     size_table,dbl_mb(k_triples_table)))
     2     call errquit('ccsd_t_restart: triples table RTDB put',
     3                  0,RTDB_ERR)
c
      enddo ! outer_virtual_index
c      next = nxtask(-nprocs,1)
c
c     Sum triples correction when finished
c
      t_energy = 0.0d0
      do i = 1,size_table
        t_energy = t_energy + dbl_mb(k_triples_table+i-1)
      enddo
#ifdef DEBUG_PRINT
      if (nodezero) write(6,*) 't_energy = ',t_energy
c
      write(6,9000)
      do i = 1,size_table
        write(6,9001) node,i,dbl_mb(k_triples_table+i-1)
      enddo
      write(6,9000)
      call util_flush(6)
#endif
c
      if (.not.ma_pop_stack(l_triples_table))
     1  call errquit("ccsd_t_restart: MA problem",0,MA_ERR)
c
c     Reset starting index once complete?
c
c      ccsd_t_restart_begin = 1
c      if (.not.rtdb_put(rtdb,'tce:ccsd_t_restart_begin',mt_int,1,
c     1                  ccsd_t_restart_begin))
c     2  call errquit('ccsd_t_restart: failed writing to RTDB',1,
c     3  RTDB_ERR)
c
      return
 9000 format(1x,'---------------------------------------------------')
 9001 format(1x,'node ',i4,' table(',i4,') = ',f25.14)
      end
