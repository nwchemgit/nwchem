      subroutine argos_timing
c
      implicit none
c
#include "argos.fh"
c
      integer timer_calls
      real*8 timer_cpu_average,timer_wall_average
      external timer_calls
      external timer_cpu_average,timer_wall_average
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' TIMING ANALYSIS',//,t20,' Number of calls',
     + t44,' Average times per call',
     + //,t40,'        Wall',t55,'         CPU',/)
      write(lfnout,1001) timer_calls(201),timer_wall_average(201),
     + timer_cpu_average(201)
 1001 format(' MD step',t20,i10,t40,f12.3,t55,f12.3)
c
      if(timer_calls(206).gt.0) then
      write(lfnout,2001) timer_calls(206),timer_wall_average(206),
     + timer_cpu_average(206)
 2001 format(' FFT',t20,i10,t40,f12.3,t55,f12.3)
      endif
      close(unit=lfnout,status='keep')
      endif
c
      return
      end
c $Id$
