       subroutine moints_stats_print( algo )
C$Id$
       implicit none
#include "msgids.fh"
#include "global.fh"
       character*(*) algo
#include "moints_stats.fh"       
       integer NSTATS
       parameter(NSTATS=21)
       double precision mv(NSTATS),amv(NSTATS)
       integer i,npass,num_nodes,nflop
       double precision taskavg, synchavg
       double precision aggtasks, aggsynch
       integer nsynchs, ntasks, nsynch0

       num_nodes = ga_nnodes()
       call dfill( NSTATS, 0.d0, mv, 1 )
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
       mv(17) = mi_aggtask
       mv(19) = mi_nsynchs
       mv(20) = mi_aggsynch
       mv(21) = mi_tio

       call ga_sync()

       call ga_dgop(msg_moints_stats_1,mv,nstats,'+')
       ntasks = mv(15)
       aggtasks = mv(17)
       nsynchs = mv(19)
       aggsynch = mv(20)

       if (nsynchs.eq.0) then
         nsynchs = 1
         aggsynch = 0.d0
         mi_minsynch = 0.d0
         mi_maxsynch = 0.d0
       endif

       call ga_dgop(msg_moints_stats_2,mi_maxsynch,1,'max')
       call ga_dgop(msg_moints_stats_3,mi_minsynch,1,'min')
       call ga_dgop(msg_moints_stats_4,mi_maxtask,1,'max')
       call ga_dgop(msg_moints_stats_5,mi_mintask,1,'min')
         
       do i=1,nstats
         amv(i) = mv(i)/num_nodes
       enddo

       taskavg = aggtasks/ntasks
       synchavg = aggsynch/nsynchs
       nflop = mi_flop1
       nsynch0 = mi_nsynchs
       
       if (algo.eq.'semi') then
         if (ga_nodeid().eq.0) write(6,337) npass, ntasks, nsynch0,
     $                mi_flop1, mv(13), amv(13),
     $                mi_tint, mv(3), amv(3),
     $                mi_t1, mv(4), amv(4),
     $                mi_t2k, mv(7), amv(7),
     $                mi_tio, mv(21), amv(21),
     $                mi_synch, mv(14), amv(14),
     $                mi_thalf, mv(2), amv(2),
     $                mi_ttotal, mv(1), amv(1)
 337     format(//,
     $          /,20x,'Statistics + Timings',
     $          /,20x,23('-'),
     $          /,5x,'Number of passes:',14x,i10,
     $          /,5x,'Number of tasks:',15x,i10,
     $          /,5x,'Number of synchs:',14x,i10,
     $          /,5x,60('-'),
     $          /,30x,'Node 0',6x,'Aggregate',6x,'Average',
     $          /,5x,60('-'),
     $          /,5x,'MFlop count: 1st idx',1x,f11.2,2(3x,f11.2),
     $          /,5x,'Integrals:',8x,3(3x,f11.2),
     $          /,5x,'1st index:',8x,3(3x,f11.2),
     $          /,5x,'2nd index K:',6x,3(3x,f11.2),
     $          /,5x,'Disk I/O:',9x,3(3x,f11.2),
     $          /,5x,'Synchronization:',2x,3(3x,f11.2),
     $          /,5x,'Four-fold shell:',2x,3(3x,f11.2),
     $          /,5x,'Total:',12x,3(3x,f11.2),
     $          /,5x,60('-'))

       else if (algo.eq.'sixfold') then
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
     $          /,20x,23('-'),
     $          /,5x,'Number of passes:',14x,i10,
     $          /,5x,'Number of tasks:',15x,i10,
     $          /,5x,60('-'),
     $          /,30x,'Node 0',6x,'Aggregate',6x,'Average',
     $          /,5x,60('-'),
     $          /,5x,'Integrals:',8x,3(3x,f11.2),
     $          /,5x,'1st index J:',6x,3(3x,f11.2),
     $          /,5x,'2nd index J:',6x,3(3x,f11.2),
     $          /,5x,'3rd & 4th index J:',3(3x,f11.2),
     $          /,5x,'1st index K:',6x,3(3x,f11.2),
     $          /,5x,'2nd index K:',6x,3(3x,f11.2),
     $          /,5x,'3rd & 4th index K:',3(3x,f11.2),
     $          /,5x,'Synchronization:',2x,3(3x,f11.2),
     $          /,5x,'Four-fold shell:',2x,3(3x,f11.2),
     $          /,5x,'Total:',12x,3(3x,f11.2),
     $          /,5x,60('-'))

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
     $          /,20x,23('-'),
     $          /,5x,'Number of passes:',14x,i10,
     $          /,5x,'Number of tasks:',15x,i10,
     $          /,5x,60('-'),
     $          /,30x,'Node 0',6x,'Aggregate',6x,'Average',
     $          /,5x,60('-'),
     $          /,5x,'MFlop count: 1st idx',1x,f11.2,2(3x,f11.2),
     $          /,5x,'Integrals:',8x,3(3x,f11.2),
     $          /,5x,'1st index:',8x,3(3x,f11.2),
     $          /,5x,'2nd index K:',6x,3(3x,f11.2),
     $          /,5x,'3rd index K:',6x,3(3x,f11.2),
     $          /,5x,'4th index K:',6x,3(3x,f11.2),
     $          /,5x,'2nd index J:',6x,3(3x,f11.2),
     $          /,5x,'3rd & 4th index J:',3(3x,f11.2),
     $          /,5x,'Synchronization:',2x,3(3x,f11.2),
     $          /,5x,'Four-fold shell:',2x,3(3x,f11.2),
     $          /,5x,'Total:',12x,3(3x,f11.2),
     $          /,5x,60('-'))
       endif
       if (ga_nodeid().eq.0) then
         write(6,933)  mi_maxsynch, mi_minsynch, synchavg,
     $                 mi_maxtask, mi_mintask, taskavg
 933     format(//,5x,57('-'),
     $          /,21x,'max',12x,'min',12x,'avg',
     $          /,5x,57('-'),
     $          /,5x,'Synch:',5x,3(f10.4,5x),
     $          /,5x,'Task:',6x,3(f10.4,5x),
     $          /,5x,57('-'))
       endif
       if (ga_nodeid().eq.0) call util_flush(6)

       return
       end



       subroutine moints_stats_clear()
       implicit none
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
       mi_maxsynch = 0.d0
       mi_minsynch = 1.d24
       mi_aggsynch = 0.d0
       mi_nsynchs = 0.d0
       mi_maxtask = 0.d0
       mi_mintask = 1.d24
       mi_aggtask = 0.d0
       mi_ntasks = 0.d0
       mi_tio = 0.d0

       return
       end


