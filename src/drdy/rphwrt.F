C $Id$
#if defined(LINUX) && !defined(XLFLINUX) && !defined(CRAYFORTRAN) && !defined(__flang__)
#define APPEND_OK
#endif
      Subroutine drdy_rphwrt (iop)
C
C Write out RPH interpolation information.
C
C    IOP = 0, first call to RPHSET
C    IOP = 1,2, write out reactant information
C    IOP = 3,4, write out product information
C    IOP = 5, write out saddle point information
C    IOP = 6, set up to write out RPH information for points along
C             the MEP
C    IOP = 7, write out RPH information for points along the MEP
C
      Implicit None
#include "errquit.fh"
C
#include "drdyP.fh"
#include "stdio.fh"
C
*::local:
*::passed:
      Integer iop
      Character*255 filename ! we want to give the file 30 a name
*::local:
#if !(defined(APPEND_OK) || defined(F90_IO))
      character*80 dummy_card
#endif
      Logical isopen
      Integer i,j,llopt(40),lifreq,nfreq,nend
      Double Precision swrt
C
      nfreq  = 0
      nend   = 0
      lifreq = 0
c
c get a unique file name for fu30
c
      call util_file_name('file30',.false.,.false.,filename)
C
      if (iop.eq.0) then
c
C First call; write comment card and options
C
         isopen = .false.
         inquire(unit=fu30,opened=isopen)
         if (isopen) then
            write(fu6,*)' deleting old unit fu30 data'
            close (unit=fu30, status='delete')
         endif
         open (unit=fu30,status='new',form='formatted',file=filename)
         write (fu30,3000) ctitle
C
C Options:
C    LLOPT(1)  = 2, print out titles, options, summary information, 
C                   and echo input as it is read
C    LLOPT(2)  = maximum number of points used in Lagrange interpolation
C    LLOPT(3)  = 1, write coordinates and gradients in mass-weighted
C    LLOPT(4)  =-2, packed F matrix written in mass-weighted coordinates
C    LLOPT(5)  = 0, no anharmonicity data
C    LLOPT(6)  = 2, use gradient on save grid and fit to quadratic
C
         do i = 1,40
            llopt(i) = 0
         enddo
         llopt(1) = 2
         llopt(2) = 3
         llopt(3) = 1
         llopt(4) = -2
         llopt(5) = 0
         llopt(6) = 2
         llopt(7) = 0
         llopt(8) = 0
         write (fu30,3005) llopt
      elseif (iop.lt.5) then
c
c reactants or products
c
         isopen = .false.
         inquire(unit=fu30,opened=isopen)
         if (isopen) then
            write(fu6,*)' deleting old unit fu30 data'
            close (unit=fu30, status='keep')
         endif
#if defined(APPEND_OK)
         open (unit=fu30,access='append',form='formatted',
     &         file=filename)
#endif
#if defined(F90_IO)
         open (unit=fu30,position='append',form='formatted',
     $         file=filename)
#endif
#if !(defined(APPEND_OK) || defined(F90_IO))
         open (unit=fu30,form='formatted',file=filename)
         rewind fu30
00001    continue
         read (fu30,'(a)',end=00002) dummy_card
         goto 00001
00002    backspace fu30
#endif
         write (fu30,3010) iop
         swrt = 0.0d00
         if (iop.eq.1) then
            vzero = v
            v = v - vzero
            if (lgs(39).eq.0) then
               write (fu30,3020) swrt,v,vzero
            else
               vzerosp = vspc
               vspc = vspc - vzerosp
               write (fu30,3021) swrt,v,vzero,vspc,vzerosp
            endif
         else if (iop.eq.3) then
            v = v - vzero
            if (lgs(39).eq.0) then
               write (fu30,3020) swrt,v
            else
               vspc = vspc - vzerosp
               write (fu30,3020) swrt,v,vspc
            endif
         else
            write (fu30,3020) swrt,0.0d00,0.0d00
         endif
         lifreq = 0
         if (iop.eq.2 .or. iop.eq.4) lifreq = nf(iop-1)
         if (icode(iop).eq.1) then
            write (fu30,*) ' '
         elseif (icode(iop).eq.2) then
            write (fu6,6000) iop
            call errquit('drdy_rphwrt: fatal error ',911, INPUT_ERR)
         else
            nfreq = nf(iop)
            nend = 3*nratom(iop)
            write (fu30,3005) (lifreq+i,i=1,nfreq)
C  write out packed hessian in mass-weighted coordinates
            write (fu30,3040)
            do i = 1, nend
               write (fu30,3020) (fsv(i,j),j=1,i)
            enddo
         endif
      elseif (iop.eq.5) then
c
c saddle point
c
        isopen = .false.
        inquire(unit=fu30,opened=isopen)
        if (isopen) then
          write(fu6,*)' deleting old unit fu30 data'
          close (unit=fu30, status='keep')
        endif
#if defined(APPEND_OK)
        open (unit=fu30,access='append',form='formatted',
     &        file=filename)
#endif
#if defined(F90_IO)
        open (unit=fu30,position='append',form='formatted',
     &        file=filename)
#endif
#if !(defined(APPEND_OK) || defined(F90_IO))
        open (unit=fu30,form='formatted',file=filename)
        rewind fu30
00003   continue
        read (fu30,'(a)',end=00004) dummy_card
        goto 00003
00004   backspace fu30
#endif
        write (fu30,3050)
         nfreq = nf(5)
         if(lgs(39).eq.0) then
            write (fu30,3020) s,v-vzero
         else
            write (fu30,3020) s,v-vzero,vspc-vzerosp
         endif
c  write out mass-weighted x
C         write (fu30,3030) (x(i),i=1,n3)
         write (fu30,3005) (i,i=1,nfreq)
C  write out packed hessian in mass-weighted coordinates
         write (fu30,3040)
         do i = 1, n3
            write (fu30,3020) (f(i,j),j=1,i)
         enddo
      elseif (iop.eq.6) then
*
         isopen = .false.
         inquire(unit=fu30,opened=isopen)
         if (isopen) then
            write(fu6,*)' deleting old unit fu30 data'
            close (unit=fu30, status='keep')
         endif
#if defined(APPEND_OK)
         open (unit=fu30,access='append',form='formatted',
     &         file=filename)
#endif
#if defined(F90_IO)
         open (unit=fu30,position='append',form='formatted',
     &         file=filename)
#endif
#if !(defined(APPEND_OK) || defined(F90_IO))
         open (unit=fu30,form='formatted',file=filename)
         rewind fu30
00005    continue
         read (fu30,'(a)',end=00006) dummy_card
         goto 00005
00006    backspace fu30
#endif
C
         write (fu30,3060)
      else
c
c points along mep
c
         isopen = .false.
         inquire(unit=fu30,opened=isopen)
         if (isopen) then
            write(fu6,*)' deleting old unit fu30 data'
            close (unit=fu30, status='keep')
         endif
#if defined(APPEND_OK)
         open (unit=fu30,access='append',form='formatted',
     &         file=filename)
#endif
#if defined(F90_IO)
         open (unit=fu30,position='append',form='formatted',
     &         file=filename)
#endif
#if !(defined(APPEND_OK) || defined(F90_IO))
         open (unit=fu30,form='formatted',file=filename)
         rewind fu30
00007    continue
         read (fu30,'(a)',end=00008) dummy_card
         goto 00007
00008    backspace fu30
#endif
         if (icode(5).eq.3) then
            nfreq = n3 - 6
         else
            nfreq = n3 - 7
         endif
         if(lgs(39).eq.0) then
            write (fu30,3070) s,v-vzero
         else
            write (fu30,3070) s,v-vzero,vspc-vzerosp
         endif
c  write out mass-weighted x
         write (fu30,3030) (x(i),i=1,n3)
c  write out mass-weighted dx
         write (fu30,3080) (dx(i),i=1,n3)
         write (fu30,3005) (i,i=1,nfreq)
c  write out packed hessian in mass-weighted coordinates
         write (fu30,3040)
         do i = 1, n3
            write (fu30,3020) (f(i,j),j=1,i)
         enddo
      endif
      close (fu30)
      return
C
3000  Format (a80)
3005  Format (1x,i2,39i3)
3010  Format (' iop =', i5, ', s,v=')
3020  Format (1x, 1pe23.15, 3e24.15)
3021  Format (1x, 1pe23.15, 4e24.15)
3030  Format (' x='/ (1x, 1pe23.15, 3e24.15))
3040  Format (' f matrix (packed)=')
3050  Format (' saddle point,  s,v=')
3060  Format (' begin points along reaction coordinate'/ ' number of',
     *   ' points goes here')
3070  Format (' s,v='/ (1x, 1pe23.15, 3e24.15))
3080  Format (' dx='/ (1x, 1pe23.15, 3e24.15))
6000  Format (' icode=2 was chosen for iop=',i5,', but it is not',
     *   ' supported')
      END
