c
c     $Id: densityofstates.F,v 1.1 2003-12-02 23:55:34 bylaska Exp $
c

*     *********************************
*     *                               *
*     *        densityofstates        *
*     *                               *
*     *********************************
*
*

      subroutine densityofstates(filename,
     >                           eigs,neigs,alpha,npoints,emin,emax)
      implicit none
      character*(*) filename
      real*8 eigs(*)
      integer neigs
      real*8 alpha
      integer npoints
      real*8 emin,emax
      
*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)

      integer j,n,l
      real*8  f,e,de,sf
      character*255 full_filename

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
      l = index(filename,' ') - 1
      write(*,1200)
      write(*,1210) filename(1:l),alpha,emin,emax,npoints

      full_filename = filename
      call util_file_name_resolve(full_filename, .false.)
      open(unit=19,file=full_filename,form='formatted')

       
      sf = alpha/(4.0d0*datan(1.0d0))/dble(neigs)
      de = (emax-emin)/dble(npoints-1)
      do j=1,npoints
        e = emin + (j-1)*de

        f=0.0d0
        do n=1,neigs
          f = f + 1.0d0/((e-eigs(n))**2 + alpha**2)
        end do
        f = sf*f

        write(19,1310) e,f
      end do
      close(19)
      end if

      return

 1200 FORMAT(//'== DOS generation =='/)
 1210 FORMAT('  filename=',A,
     >       ', alpha=',F10.4,', (emin:emax)= (',F10.4,':',F10.4,
     >       '),  npoints=',I5)
 1310 FORMAT(2E15.5)

      end

