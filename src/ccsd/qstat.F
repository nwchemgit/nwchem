      subroutine qenter(namein,iprint)
C     $Id: qstat.F,v 2.4 1997-10-13 06:34:52 d3e129 Exp $
      Implicit NONE
      Character*(*) NameIn
      Integer IPrint
C
      Integer N, I
      Double Precision Tcp1, Tio1
C
C
#include "qdata.fh"
c
      double precision zip
      data zip/0.0d00/
*
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
         stop 20
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
      Integer Nodes, IAm, I, NDe
C
C
#include "qdata.fh"
#include "global.fh"
*
      double precision zip
      data zip/0.0d00/
*
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
     
      do 310 nde=0,min(3,nodes-1)
      if (iam.eq.nde)then
      write(*,*)
c     write(*,'(5x,a)') 'runtime statistics'
c     write(*,'(5x,a)') '------------------'
c     write(*,*)
      write(*,'(a)') ' node routine      calls     cpu     i/o'
      do 300 i=1,nrout
         write(*,'(i4,2x,a8,i9,f9.2,f8.2)')
     &      iam,name(i),ncall(i),tcp(i),tio(i)
300   continue
      write(*,'(6x,a8,9x,f9.2,f8.2)')
     &   'Total   ',tcptot,tiotot
      endif
      call ga_sync()
 310  continue
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
