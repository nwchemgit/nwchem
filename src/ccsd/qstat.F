      subroutine qenter(namein,iprint)
C     $Id$
      Implicit NONE
#include "errquit.fh"
      Character*(*) NameIn
      Integer IPrint
C
      Integer N, I
      Double Precision Tcp1, Tio1
C
C
#include "qdata.fh"
***      logical ma_verify_allocator_stuff
***      external ma_verify_allocator_stuff
c
      double precision zip
      data zip/0.0d00/
*

***      write(6,*) 'qenter:', namein
***      call util_flush(6)
***      if (.not. ma_verify_allocator_stuff()) call errquit(namein,0,
***     &       MA_ERR)

      call cputime(tcp1)
      tio1=zip
*
* update caller
*
      if(level.gt.0) then
         n=no(level)
         tcp(n)=tcp(n)+tcp1-tcp0
         tio(n)=tio(n)+tio1-tio0
         if(iprint.gt.0) then
            write(*,*)
            write(*,*) '****** ',name(n),' calling ',namein,' ******'
            write(*,*)
         end if
      else if(iprint.gt.0) then
         write(*,*)
         write(*,*) '****** main calling ',namein,' ******'
         write(*,*)
      end if
*
* update called
*
      level=level+1
      if(level.gt.mxlev) then
         write(*,*) 'overflow in subroutine stack',
     &              ' in routine qentry'
         stop 20
      end if
      n=0
      do 100 i=1,nrout
         if(namein.eq.name(i)) n=i
100   continue
      if(n.eq.0) then
         nrout=nrout+1
         if(nrout.gt.mxrout) then
            write(*,*) 'overflow in subroutine list',
     &                 ' in routine qentry'
            stop 20
         end if
         n=nrout
         name(n)=namein
         tcp(n)=zip
         tio(n)=zip
         ncall(n)=0
      end if
      ncall(n)=ncall(n)+1
      no(level)=n
      call cputime(tcp0)
      tio0=zip
      return
      end
      subroutine qexit(namein,iprint)
      Implicit NONE
      Character*(*) NameIn
      Integer IPrint
C
      Double Precision Tcp1, Tio1
      Integer M, N
C
C
#include "qdata.fh"
*
      double precision zip
      data zip/0.0d00/
c
      call cputime(tcp1)
      tio1=zip
      if(level.eq.0) then
         write(*,*) 'qexit called with empty subroutine stack'
         call errquit(' qexit failure ',0,0)
      end if
*
* update exited routine
*
      n=no(level)
      if(name(n).ne.namein) then
         write(*,*) 'attempt to exit ',namein,' when in ',name(n)
         stop 20
      end if
      tcp(n)=tcp(n)+tcp1-tcp0
      tio(n)=tio(n)+tio1-tio0
      if(iprint.gt.0) then
         if(level.gt.1) then
            m=no(level-1)
            write(*,*)
            write(*,*) '****** ',namein,' returning to ',name(m),
     &                 ' ******'
            write(*,*)
         else
            write(*,*)
            write(*,*) '****** ',name(n),' returning to main ******'
            write(*,*)
         end if
      end if
*
*
*
      level=level-1
      call cputime(tcp0)
      tio0=zip
      return
      end
      subroutine qstat
      Implicit NONE
*
      Double Precision Tcp2, Tcp3, Tio2, Tio3, TcpTot, TioTot
      Double Precision tcp_min, tcp_ave, tcp_max
      Integer Nodes, IAm, I
      Logical oprint
*old:      Integer NDe
C
C
#include "qdata.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
*
      double precision zip
      data zip/0.0d00/
*
      oprint = util_print('information', print_low)
      call cputime(tcp2)
      tio2=zip
      tcptot=zip
      tiotot=zip
      nodes = ga_nnodes()
      iam = ga_nodeid()
      call ga_sync()
      do 100 i=1,nrout
      tcptot=tcptot+tcp(i)
      tiotot=tiotot+tio(i)
 100  continue
      call ga_sync()
     
*old:      do 310 nde=0,min(3,nodes-1)
*old:      if (iam.eq.nde)then
*old:      write(*,*)
*old:c     write(*,'(5x,a)') 'runtime statistics'
*old:c     write(*,'(5x,a)') '------------------'
*old:c     write(*,*)
*old:      write(*,'(a)') ' node routine      calls     cpu     i/o'
*old:      do 300 i=1,nrout
*old:         write(*,'(i4,2x,a8,i9,f9.2,f8.2)')
*old:     &      iam,name(i),ncall(i),tcp(i),tio(i)
*old:300   continue
*old:      write(*,'(6x,a8,9x,f9.2,f8.2)')
*old:     &   'Total   ',tcptot,tiotot
*old:      endif
      if (iam.eq.0.and.oprint)write(*,*)
      if (iam.eq.0.and.oprint)
     &    write(*,'(a,a)')
     &    ' routine      calls  cpu(0)   ',
     &    'cpu-min  cpu-ave  cpu-max   i/o '
cedo  fix for hangs when nrout is not constant among procs
      call ga_dgop(msg_qstat_3,nrout,1,'min')
      do 300 i=1,nrout
         tcp_min = tcp(i)
         call ga_dgop(msg_qstat_1,tcp_min,1,'min')
         tcp_ave = tcp(i)
         call ga_dgop(msg_qstat_2,tcp_ave,1,'+')
         tcp_ave = tcp_ave/dble(nodes)
         tcp_max = tcp(i)
         call ga_dgop(msg_qstat_3,tcp_max,1,'max')
         if (iam.eq.0.and.oprint) then
           write(*,'(1x,a8,i9,f9.2,f9.2,f9.2,f9.2,f8.2)')
     &         name(i),ncall(i),tcp(i),
     &         tcp_min,tcp_ave,tcp_max,tio(i)
           call util_flush(6)
         endif
300   continue
      tcp_min = tcptot
      call ga_dgop(msg_qstat_1,tcp_min,1,'min')
      tcp_ave = tcptot
      call ga_dgop(msg_qstat_2,tcp_ave,1,'+')
      tcp_ave = tcp_ave/dble(nodes)
      tcp_max = tcptot
      call ga_dgop(msg_qstat_3,tcp_max,1,'max')
      if (iam.eq.0.and.oprint) then
        write(*,'(1x,a8,9x,f9.2,f9.2,f9.2,f9.2,f8.2)')
     &      'Total   ',tcptot,tcp_min,tcp_ave,tcp_max,tiotot
        call util_flush(6)
      endif
      call ga_sync()
c
      call cputime(tcp3)
      tio3=zip
      tcp0=tcp0+tcp3-tcp2
      tio0=tio0+tio3-tio2
      return
      end
      block data qblock
      Implicit NONE
c$$$      implicit real*8 (a-h,o-z)
*
#include "qdata.fh"      
*     
      data level,nrout/0,0/
*
      end
      subroutine qstat_init
      implicit none
#include "qdata.fh"
      integer i
*
      call ifill(mxlev,0,no,1)
      Tcp0 = 0.0d00
      Tio0 = 0.0d00
      call dfill(mxrout,0.0d00,tcp   ,1)
      call dfill(mxrout,0.0d00,tio   ,1)
      call ifill(mxrout,0     ,ncall ,1)
      do i = 1,mxrout
        name(i) = ' '
      enddo
      level = 0
      nrout = 0
      end
