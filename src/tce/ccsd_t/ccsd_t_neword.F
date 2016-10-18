      integer function ccsd_t_6tasks(
     R     restricted,noab,nvab,kspin,ksym)
C
C     $Id$
C
      implicit none
      logical restricted
      integer noab,nvab
      integer kspin(*),ksym(*)
c
      integer t_p4b,t_p5b,t_p6b
      integer t_h1b,t_h2b,t_h3b
c
      ccsd_t_6tasks=0
      do t_p4b = noab+1,noab+nvab
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
            if (kspin(t_p4b)+kspin(t_p5b)+kspin(t_p6b).eq.
     3             kspin(t_h1b)+kspin(t_h2b)+kspin(t_h3b)) then
            if ((.not.restricted).or.
     1              (kspin(t_p4b)+kspin(t_p5b)+kspin(t_p6b)+
     3         kspin(t_h1b)+kspin(t_h2b)+kspin(t_h3b).le.8)) then
            if (ieor(ksym(t_p4b),ieor(ksym(t_p5b),ieor(ksym(t_p6b),
     3          ieor(ksym(t_h1b),ieor(ksym(t_h2b),ksym(t_h3b)
     P              ))))).eq.0) then
               ccsd_t_6tasks=ccsd_t_6tasks + 1
            endif
            endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      return
      end
c
      subroutine ccsd_t_neword(tot_task, 
     w     w_thresh,task_thresh,
     R     restricted,noab,nvab,kspin,ksym,
     K     krange,kaux,klist)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      logical restricted
      integer w_thresh,task_thresh
      integer tot_task
      integer noab,nvab
      integer kspin(*),ksym(*)
      integer krange(*)
      integer kaux(7,*),klist(7,*)
c
      integer t_p4b,t_p5b,t_p6b
      integer t_h1b,t_h2b,t_h3b
      integer i,j,m
      integer wl_max,wl_min
      integer range_p4,range_p5,range_p6,
     R     range_h1,range_h2,range_h3
      integer nsplits,ii,w_in
      integer highindx
      parameter (highindx=7)
      double precision rr456,rr4561,rr45612
c
      task_thresh=tot_task
      m=0
      do t_p4b = noab+1,noab+nvab
       range_p4 = krange(t_p4b)
       do t_p5b = t_p4b,noab+nvab
        range_p5 = krange(t_p5b)
        do t_p6b = t_p5b,noab+nvab
         range_p6 = krange(t_p6b)
         rr456=range_p4*range_p5*range_p6
         do t_h1b = 1,noab
          range_h1 = krange(t_h1b)
            rr4561=rr456*range_h1
          do t_h2b = t_h1b,noab
             range_h2 = krange(t_h2b)
             rr45612 = rr4561*range_h2
           do t_h3b = t_h2b,noab
            if (kspin(t_p4b)+kspin(t_p5b)+kspin(t_p6b).eq.
     3             kspin(t_h1b)+kspin(t_h2b)+kspin(t_h3b)) then
            if ((.not.restricted).or.
     1              (kspin(t_p4b)+kspin(t_p5b)+kspin(t_p6b)+
     3         kspin(t_h1b)+kspin(t_h2b)+kspin(t_h3b).le.8)) then
            if (ieor(ksym(t_p4b),ieor(ksym(t_p5b),ieor(ksym(t_p6b),
     3          ieor(ksym(t_h1b),ieor(ksym(t_h2b),ksym(t_h3b)
     P              ))))).eq.0) then
            m=m+1
            kaux(1,m)=t_p4b
            kaux(2,m)=t_p5b
            kaux(3,m)=t_p6b
            kaux(4,m)=t_h1b
            kaux(5,m)=t_h2b
            kaux(6,m)=t_h3b
            range_h3 = krange(t_h3b)
            kaux(highindx,m)=rr45612*range_h3
            endif
            endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
#ifdef DEBUG
      do i=1,tot_task
         write(6,12) i,(kaux(j,i),j=1,highindx)
 12      format(I7,6I8,I20)
      enddo
#endif
c
      if(m.ne.tot_task)  then
         write(6,*) ' m ',m,' tot_task ',tot_task
         call errquit("incompatibility_task",6,MA_ERR)
      endif
c
      wl_max=0
      do m=1,tot_task
        wl_max=max(kaux(highindx,m),wl_max)
      enddo
c
      wl_min=wl_max
      do m=1,tot_task
        wl_min=min(kaux(highindx,m),wl_min)
      enddo
c     all tiles the same. no ordering needed
#if 1
      if(((wl_max-wl_min)*100d0)/wl_max.lt.1d0) then 
      do m=1,tot_task
         do j=1,highindx
            klist(j,m)=kaux(j,m)
         enddo
      enddo
      return
      endif
#endif

c
c simple sorting
c
c
      do m=1,tot_task
         do j=1,highindx
            klist(j,m)=-1
         enddo
      enddo
      m=0
      if(ga_nodeid().eq.0) then
         write(6,33) ' total no. of tasks ',tot_task
         write(6,33) ' total no. of tasks / no. procs ',
     c        tot_task/ga_nnodes()
         write(6,33) ' wl_min ',wl_min,(wl_min)**(1d0/6d0)
         write(6,33) ' wl_max ',wl_max,(wl_max)**(1d0/6d0)
      endif
   33 format(a32,i12,f14.1)

#ifdef DEBUG
      write(6,*) ' sort1 '
      do i=1,tot_task
         write(6,11) i,(klist(j,i),j=1,6)
      enddo
#endif
      nsplits=16
      do ii=nsplits,1,-1
         w_in=wl_min+((wl_max-wl_min)*(ii-1))/nsplits
          call ccsd_sillysort(w_in,
     R     kaux,klist, tot_task,m,highindx)
       enddo
      call ccsd_sillysort(0,
     R     kaux,klist, tot_task,m,highindx)
c     check that everybody has been invalidated
      do i=1,tot_task
         if(kaux(highindx,i).ge.0) then
            write(6,*) ' kaux element ',i,
     R           'not sorted. val = ',kaux(highindx,i)
            call errquit('neworder error ',0,0)
            endif
      enddo
c
c     figure out task with w lt w_thresh
c
c      w_thresh=w_thresh**6

      do i=1,tot_task
         if(klist(highindx,i).lt.w_thresh) then
            task_thresh=i
            goto 1973
         endif
      enddo
c
c      
 1973 continue
cbaad
c     task_thresh must be gt ga_nnodes because of the static nxtask for the first task
cbaad      task_thresh=max(task_thresh,ga_nnodes()+1) 
      if(ga_nodeid().eq.0) then
         if(task_thresh.lt.ga_nnodes()+1) then
            write(6,*) ' warning: CCSD(T) oversubscribed '
            write(6,'(a,i8)') '  warning: tip: use nproc <',task_thresh
         endif
         write(6,33) ' thresh for no. of tasks ',task_thresh
      endif
#ifdef DEBUG
      write(6,*) ' sorting ok '
      do m=1,tot_task
         write(6,11) m,(klist(j,m),j=1,6)
 11      format(I7,6I8)
      enddo
#endif
      return
      end
      subroutine ccsd_sillysort(value,
     R     kaux,klist,n,found,highindx )
      implicit none
c
c very simple sorting
c
      integer value
      integer n
      integer found
      integer highindx
      integer kaux(highindx,*)
      integer klist(highindx,*)
c
      integer i,j,m
c
      do m=1,n
         if(kaux(highindx,m).gt.value) then
            found=found+1
            do j=1,highindx
               klist(j,found)=kaux(j,m)
            enddo
c     invalidate
            kaux(highindx,m)=-99
         endif
      enddo
      return
      end
