       subroutine moints_stats_print( algo )
       implicit none
#include "global.fh"
       character*(*) algo
#ifdef MOINTS_STATS
#include "moints_stats.fh"       
       integer NSTATS
       integer MSG_MOINTS_STATS
       parameter(NSTATS=15)
       parameter(MSG_MOINTS_STATS=1963)
       double precision mv(NSTATS),amv(NSTATS)
       integer i,npass,ntasks,num_nodes,nflop

       num_nodes = ga_nnodes()
       npass = mi_npass

       mv(1) = mi_ttotal
       mv(2) = mi_thalf
       mv(3) = mi_tint
       mv(4) = mi_t1
       mv(5) = mi_t1j
       mv(6) = mi_t1k
       mv(7) = mi_t2k
       mv(8) = mi_t3k
       mv(9) = mi_t4k
       mv(10) = mi_t2j
       mv(11) = mi_t34j
       mv(12) = mi_t34k
       mv(13) = mi_flop1
       mv(14) = mi_synch
       mv(15) = mi_ntasks

       call ga_sync()
       call ga_dgop(MSG_MOINTS_STATS,mv,NSTATS,'+')
       do i=1,NSTATS
         amv(i) = mv(i)/num_nodes
       enddo
       ntasks = mv(15)
       nflop = mi_flop1

       if (algo.eq.'sixfold') then
         if (ga_nodeid().eq.0) write(6,333) npass, ntasks,
     $                mi_tint, mv(3), amv(3),
     $                mi_t1j, mv(5), amv(5),
     $                mi_t2j, mv(10), amv(10),
     $                mi_t34j, mv(11), amv(11),
     $                mi_t1k, mv(6), amv(6),
     $                mi_t2k, mv(7), amv(7),
     $                mi_t34k, mv(12), amv(12),
     $                mi_synch, mv(14), amv(14),
     $                mi_thalf, mv(2), amv(2),
     $                mi_ttotal, mv(1), amv(1)
 333     format(//,
     $          /,20x,'Statistics + Timings',
     $          /,20x,20('-'),
     $          /,5x,'Number of passes:',14x,i10,
     $          /,5x,'Number of tasks:',15x,i10,
     $          /,5x,57('-'),
     $          /,30x,'Node 0',5x,'Aggregate',5x,'Average',
     $          /,5x,57('-'),
     $          /,5x,'Integrals:',8x,3(3x,f10.2),
     $          /,5x,'1st index J:',6x,3(3x,f10.2),
     $          /,5x,'2nd index J:',6x,3(3x,f10.2),
     $          /,5x,'3rd & 4th index J:',3(3x,f10.2),
     $          /,5x,'1st index K:',6x,3(3x,f10.2),
     $          /,5x,'2nd index K:',6x,3(3x,f10.2),
     $          /,5x,'3rd & 4th index K:',3(3x,f10.2),
     $          /,5x,'Synchronization:',2x,3(3x,f10.2),
     $          /,5x,'Four-fold shell:',2x,3(3x,f10.2),
     $          /,5x,'Total:',12x,3(3x,f10.2),
     $          /,5x,57('-'))

       else
         if (ga_nodeid().eq.0) write(6,334) npass, ntasks,
     $                mi_flop1, mv(13), amv(13),
     $                mi_tint, mv(3), amv(3),
     $                mi_t1, mv(4), amv(4),
     $                mi_t2k, mv(7), amv(7),
     $                mi_t3k, mv(8), amv(8),
     $                mi_t4k, mv(9), amv(9),
     $                mi_t2j, mv(10), amv(10),
     $                mi_t34j, mv(11), amv(11),
     $                mi_synch, mv(14), amv(14),
     $                mi_thalf, mv(2), amv(2),
     $                mi_ttotal, mv(1), amv(1)
 334     format(//,
     $          /,20x,'Statistics + Timings',
     $          /,20x,20('-'),
     $          /,5x,'Number of passes:',14x,i10,
     $          /,5x,'Number of tasks:',15x,i10,
     $          /,5x,57('-'),
     $          /,30x,'Node 0',5x,'Aggregate',5x,'Average',
     $          /,5x,57('-'),
     $          /,5x,'MFlop count: 1st idx',1x,f10.2,2(3x,f10.2),
     $          /,5x,'Integrals:',8x,3(3x,f10.2),
     $          /,5x,'1st index:',8x,3(3x,f10.2),
     $          /,5x,'2nd index K:',6x,3(3x,f10.2),
     $          /,5x,'3rd index K:',6x,3(3x,f10.2),
     $          /,5x,'4th index K:',6x,3(3x,f10.2),
     $          /,5x,'2nd index J:',6x,3(3x,f10.2),
     $          /,5x,'3rd & 4th index J:',3(3x,f10.2),
     $          /,5x,'Synchronization:',2x,3(3x,f10.2),
     $          /,5x,'Four-fold shell:',2x,3(3x,f10.2),
     $          /,5x,'Total:',12x,3(3x,f10.2),
     $          /,5x,57('-'))
       endif
       if (ga_nodeid().eq.0) call util_flush(6)

#endif /* MOINTS_STATS */
       return
       end



       subroutine moints_stats_clear()
       implicit none
#include "global.fh"
#ifdef MOINTS_STATS
#include "moints_stats.fh"
c
       mi_npass = 0.d0
       mi_ttotal = 0.d0
       mi_thalf = 0.d0
       mi_tint = 0.d0
       mi_t1 = 0.d0 
       mi_t1j = 0.d0 
       mi_t1k = 0.d0
       mi_t2k = 0.d0
       mi_t3k = 0.d0
       mi_t4k = 0.d0
       mi_t2j = 0.d0
       mi_t34j = 0.d0
       mi_t34k = 0.d0
       mi_flop1 = 0.d0
       mi_synch = 0.d0
       mi_ntasks = 0.d0

#endif /* MOINTS_STATS */
       return
       end









