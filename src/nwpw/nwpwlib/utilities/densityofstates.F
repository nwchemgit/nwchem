c
c     $Id$
c

*     *********************************
*     *                               *
*     *        densityofstates        *
*     *                               *
*     *********************************
*
      subroutine densityofstates(filename,
     >                           eigs,weight,neigs,
     >                           sign,alpha,npoints,emin,emax)
      implicit none
      character*(*) filename
      real*8 eigs(*),weight(*)
      integer neigs
      real*8 sign
      real*8 alpha
      integer npoints
      real*8 emin,emax
      integer ne_in(2)

#include "bafdecls.fh"
#include "msgtypesf.h"
#include "btdb.fh"
#include "errquit.fh"
      
*     **** local variables ***
      integer taskid,MASTER,rtdb
      parameter (MASTER=0)

      integer j,n,l
      real*8  f,e,de,sf,ff
      character*255 full_filename,prefilename

      integer  control_rtdb
      external control_rtdb


      call Parallel_taskid(taskid)
      rtdb = control_rtdb()

      if (.not.btdb_cget(rtdb,'nwpw:dos:filename',1,prefilename)) then
         call util_file_name(filename,.true.,.false.,full_filename)
      else
        l = index(prefilename,' ') - 1
        full_filename = prefilename(1:l)//'.'//filename
        call util_file_name_resolve(full_filename, .false.)
      end if

      if (taskid.eq.MASTER) then
      l = index(full_filename,' ') - 1
      write(*,1200)
      write(*,1210) full_filename(1:l),sign,alpha,emin,emax,npoints

      open(unit=19,file=full_filename,form='formatted')

       
      ff = 0.0d0
      sf = alpha/(4.0d0*datan(1.0d0))
      de = (emax-emin)/dble(npoints-1)
      do j=1,npoints
        e = emin + (j-1)*de

        f=0.0d0
        do n=1,neigs
          f = f + weight(n)/((e-eigs(n))**2 + alpha**2)
        end do
        f = sf*f

        if ((j.eq.1).or.(j.eq.npoints)) then
           ff = ff + 0.5d0*f*de
        else
           ff = ff + f*de
        end if

        write(19,1310) e,f*sign,ff*sign
      end do
      close(19)
      end if

      return

 1200 FORMAT(//'== DOS generation =='/)
 1210 FORMAT(' filename:',A,
     >       /' sign=',F8.2,
     >       ', alpha=',F8.4,', (emin:emax)= (',F10.4,':',F10.4,
     >       '), npoints=',I5)
 1310 FORMAT(3E15.5)

      end


*     *********************************
*     *                               *
*     *    writedensityofstates       *
*     *                               *
*     *********************************
*
      subroutine writedensityofstates(filename,
     >                           sign,alpha,npoints,emin,emax,dstates)
      implicit none
      character*(*) filename
      real*8 sign
      real*8 alpha
      integer npoints
      real*8 emin,emax
      real*8 dstates(*)

#include "bafdecls.fh"
#include "msgtypesf.h"
#include "btdb.fh"
#include "errquit.fh"
      
*     **** local variables ***
      integer taskid,MASTER,rtdb
      parameter (MASTER=0)

      integer j,n,l
      real*8  f,e,de,ff
      character*255 full_filename,prefilename

      integer  control_rtdb
      external control_rtdb


      call Parallel_taskid(taskid)
      rtdb = control_rtdb()

      if (.not.btdb_cget(rtdb,'nwpw:dos:filename',1,prefilename)) then
         call util_file_name(filename,.true.,.false.,full_filename)
      else
        l = index(prefilename,' ') - 1
        full_filename = prefilename(1:l)//'.'//filename
        call util_file_name_resolve(full_filename, .false.)
      end if

      if (taskid.eq.MASTER) then
      l = index(full_filename,' ') - 1
      write(*,1200)
      write(*,1210) full_filename(1:l),sign,alpha,emin,emax,npoints

      open(unit=19,file=full_filename,form='formatted')

       
      ff = 0.0d0
      de = (emax-emin)/dble(npoints-1)
      do j=1,npoints
        e = emin + (j-1)*de
       
        if ((j.eq.1).or.(j.eq.npoints)) then
           ff = ff + 0.5d0*dstates(j)*de
        else
           ff = ff + dstates(j)*de
        end if
        write(19,1310) e,dstates(j),ff
      end do
      close(19)
      end if

      return

 1200 FORMAT(//'== DOS generation =='/)
 1210 FORMAT(' filename:',A,
     >       /' sign=',F8.2,
     >       ', alpha=',F8.4,', (emin:emax)= (',F10.4,':',F10.4,
     >       '), npoints=',I5)
 1310 FORMAT(3E15.5)

      end



*     *********************************
*     *                               *
*     *        adddensityofstates     *
*     *                               *
*     *********************************
*
      subroutine adddensityofstates(eigs,weight,neigs,
     >                              sign,alpha,npoints,emin,emax,
     >                              dstates)
      implicit none
      real*8 eigs(*),weight(*)
      integer neigs
      real*8 sign
      real*8 alpha
      integer npoints
      real*8 emin,emax
      real*8 dstates(*) 

*     **** local variables ***
      integer j,n,l
      real*8  f,e,de,sf
       
      sf = alpha/(4.0d0*datan(1.0d0))
      de = (emax-emin)/dble(npoints-1)
      do j=1,npoints
        e = emin + (j-1)*de

        f=0.0d0
        do n=1,neigs
          f = f + weight(n)/((e-eigs(n))**2 + alpha**2)
        end do
        f = sf*f
        dstates(j) = dstates(j) + f*sign

      end do

      return
      end

