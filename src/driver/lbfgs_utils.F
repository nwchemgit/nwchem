*
* $Id$
      subroutine lbfgs_compute_info(n,x,x0,g,xmax,xrms,gmax,grms)
      implicit none
#include "rtdb.fh"
      integer n
      double precision x(n),x0(n),g(n)
c
c     Compute stuff used for printing and convergence tests
c
      integer i
      double precision grms,gmax,xrms,xmax
      double precision dx
c
      double precision gmax_tol,grms_tol
      double precision xmax_tol,xrms_tol

      grms = 0d0
      gmax  = 0d0
      xmax  = 0d0
      xrms  = 0d0
c
      do i = 1, n
         grms = grms + g(i)*g(i)
         gmax  = max(gmax, abs(g(i)))
         dx = x(i) - x0(i)
         xrms = xrms + dx*dx
         xmax  = max(xmax, abs(dx))
      enddo
      grms = sqrt(grms/dble(n))
      xrms = sqrt(xrms/dble(n))

c
      end

       function lbfgs_converged(xmax,xrms,gmax,grms,tol)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
c
      logical lbfgs_converged
      double precision grms,gmax,xrms,xmax,tol(4)
c
      double precision gmax_tol,grms_tol
      double precision xmax_tol,xrms_tol
c
      gmax_tol = tol(1)
      grms_tol = tol(2)
      xmax_tol = tol(3)
      xrms_tol = tol(4)
c
      lbfgs_converged = 
     $     ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol) .and.
     $     (xrms .lt. xrms_tol)  .and. (xmax .lt. xmax_tol))
     $     .or. 
     $     ((gmax.lt.0.01d0*gmax_tol) .and. (grms.lt.0.01d0*grms_tol))


c
      end

      subroutine lbfgs_print_info(istep,
     +                             e,de,xmax,xrms,gmax,grms,
     +                             tol)
      implicit none
#include "util.fh"
#include "geom.fh"
#include "inp.fh"
#include "errquit.fh"
      integer n
c
c     Compute stuff used for printing and convergence tests
c
      integer istep
      double precision e,de,grms,gmax,xrms,xmax,tol(4)
      character*1 mark
      character*9 cvg1, cvg2, cvg3, cvg4

      double precision gmax_tol,grms_tol
      double precision xmax_tol,xrms_tol
      character*255 filename, dir
      character*30 pname
c
      pname = "lbfgs_print_info"
c
      gmax_tol = tol(1)
      grms_tol = tol(2)
      xmax_tol = tol(3)
      xrms_tol = tol(4)

      cvg1 = ' '
      cvg2 = ' '
      cvg3 = ' '
      cvg4 = ' '
      if (gmax .lt. gmax_tol) cvg1 = '     ok  '
      if (grms .lt. grms_tol) cvg2 = '     ok  '
      if (xrms .lt. xrms_tol) cvg3 = '     ok  '
      if (xmax .lt. xmax_tol) cvg4 = '     ok  '

      mark = '@'
      if (istep-1 .gt. 1) mark = ' '
      write(6,1) mark, mark
      mark = '@'
      write(6,2) mark, istep-1, e, de,
     $  gmax, grms, xrms, xmax, util_wallsec(),
     $  cvg1, cvg2, cvg3, cvg4
 1    format(
     $     /,a1,' Step       Energy      Delta E   Gmax',
     $     '     Grms     Xrms     Xmax   Walltime',
     $     /,a1,' ---- ---------------- -------- --------',
     $     ' -------- -------- -------- --------')
 2    format(
     $     a1,i5,f17.8,1p,d9.1,0p,4f9.5,f9.1,/,
     $     1x,5x,17x,9x,4a9,/)
      return
      end

      subroutine lbfgs_print_xyz(istep,igeom,xyz)
      implicit none
#include "util.fh"
#include "geom.fh"
#include "inp.fh"
#include "errquit.fh"
      integer n
c
c     Compute stuff used for printing and convergence tests
c
      integer istep,igeom
      character*(*) xyz

      character*255 filename, dir
      character*30 pname
c
      pname = "lbfgs_print_xyz"
c
      if (xyz.ne.' ') then
         dir      = ' '
         filename = ' ' 
         call util_directory_name(dir, .false., 0)
         write(filename,13) dir(1:inp_strlen(dir)), 
     $        xyz(1:inp_strlen(xyz)),istep-1
 13      format(a,'/',a,'-',i3.3,'.xyz')
         write(*,*) filename
         open(88,file=filename,form='formatted',status='unknown',
     $        err=133)
         if (.not. geom_print_xyz(igeom, 88))
     $        call errquit('driver: geom_print_xyz?',0, GEOM_ERR)
         close(88,status='keep',err=133)
      endif
c
      return
c
 133  call errquit(pname//'error opening/closing xyz file',0, GEOM_ERR)
c
      end

