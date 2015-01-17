      integer function icsd_t2_8_tasks(irrep_v,irrep_t,
     R     restricted,noab,nvab,kspin,ksym)
C
C     $Id$
C
      implicit none
      logical restricted
      integer irrep_v,irrep_t
      integer noab,nvab
      integer kspin(*),ksym(*)
c
      integer p3b,p4b,p5b,p6b
      integer h1b,h2b
c
      icsd_t2_8_tasks=0
      DO p3b = noab+1,noab+nvab
         DO p4b = p3b,noab+nvab
            DO h1b = 1,noab
               DO h2b = h1b,noab
                  IF ((.not.restricted).or.(kspin(p3b)+kspin(p4b)+
     &                 kspin(h1b)+kspin(h2b).ne.8)) THEN
                     IF (kspin(p3b)+kspin(p4b) .eq. 
     E                    kspin(h1b)+kspin(h2b)) THEN
                        IF(ieor(ksym(p3b),ieor(ksym(p4b),ieor(ksym(h1b),
     K                       ksym(h2b)))) .eq. 
     I                       ieor(irrep_v,irrep_t)) then
                           DO p5b = noab+1,noab+nvab
                              DO p6b = p5b,noab+nvab
                                 IF (kspin(p5b)+kspin(p6b) .eq. 
     K                                kspin(h1b)+kspin(h2b)) THEN
                                    IF (ieor(ksym(p5b),ieor(ksym(p6b),
     I                                   ieor(ksym(h1b),ksym(h2b)))).
     E                                   eq.irrep_t) THEN
                                       icsd_t2_8_tasks=icsd_t2_8_tasks+1
                                    END IF
                                 END IF
                              END DO ! p6b
                           END DO ! p5b
                        endif
                     endif
                  endif
               END DO
            END DO
         END DO
      END DO
      return
      end
      subroutine icsd_t2_8_ord(dryrun,tot_task,
     I     irrep_v,irrep_t,
     R     restricted,noab,nvab,kspin,ksym,
     K     krange,kaux,klist)

C
C     $Id$
C
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer tot_task
      logical dryrun
      logical restricted
      integer irrep_v,irrep_t
      integer noab,nvab
      integer kspin(*),ksym(*)
      integer krange(*)
      integer kaux(5,*),klist(4,*)

c
      integer p3b,p4b,p5b,p6b
      integer h1b,h2b
      integer i,j,m
      integer ii,nsplits
      integer wl_max,wl_min,w_in
      integer range_p3,range_p4,range_p5,
     R     range_p6,range_h1,range_h2
      integer rr12,rr34,rr56,rr1234
c     A h1*h2*      p5*p6
c     B       p3*p4*p5*p6
c     C h1*h2*p3*p4
c
      m=0
      DO p3b = noab+1,noab+nvab
       range_p3 = krange(p3b)
         DO p4b = p3b,noab+nvab
            range_p4 = krange(p4b)
            rr34=range_p4*range_p3
            DO h1b = 1,noab
               range_h1 = krange(h1b)
               DO h2b = h1b,noab
                  IF ((.not.restricted).or.(kspin(p3b)+kspin(p4b)+
     &                 kspin(h1b)+kspin(h2b).ne.8)) THEN
                     IF (kspin(p3b)+kspin(p4b) .eq. 
     E                    kspin(h1b)+kspin(h2b)) THEN
                        IF(ieor(ksym(p3b),ieor(ksym(p4b),ieor(ksym(h1b),
     K                       ksym(h2b)))) .eq. 
     I                       ieor(irrep_v,irrep_t)) then
                           m=m+1
                           if(.not.dryrun) then
                              kaux(1,m)=p3b
                              kaux(2,m)=p4b
                              kaux(3,m)=h1b
                              kaux(4,m)=h2b
                              range_h2 = krange(h2b)
                              rr12=range_h1*range_h2
                              kaux(5,m)=rr12*rr34
                           endif
                        endif
                     endif
                  endif
               END DO
            END DO
         END DO
      END DO
      tot_task=m
      if(dryrun) return
      wl_max=0
      do m=1,tot_task
        wl_max=max(kaux(5,m),wl_max)
      enddo
c
      wl_min=wl_max
      do m=1,tot_task
        wl_min=min(kaux(5,m),wl_min)
      enddo
c
c simple sorting
c
c
      do m=1,tot_task
         do j=1,4
            klist(j,m)=-1
         enddo
      enddo
      m=0
      if(ga_nodeid().eq.0) then
         write(6,33) ' total no. of tasks ',tot_task
         write(6,33) ' total no. of tasks / no. procs ',
     c        tot_task/ga_nnodes()
         write(6,33) ' wl_min ',wl_min,nint((wl_min)**(1d0/4d0))
         write(6,33) ' wl_max ',wl_max,nint((wl_max)**(1d0/4d0))
      endif
 33   format(a32,2i12)
      nsplits=16
      do ii=nsplits,1,-1
         w_in=wl_min+((wl_max-wl_min)*(ii-1))/nsplits
          call ccsd_sillysort(w_in,
     R     kaux,klist, tot_task,m,5)
       enddo
      call ccsd_sillysort(0,
     R     kaux,klist, tot_task,m,5)
c     check that everybody has been invalidated
      do i=1,tot_task
         if(kaux(5,i).ge.0) then
            write(6,*) ' kaux element ',i,
     R           'not sorted. val = ',kaux(7,i)
            call errquit('neworder error ',0,0)
            endif
      enddo

      return
      end
