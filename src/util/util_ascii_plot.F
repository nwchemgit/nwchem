*     **************************************************
*     *                                                *
*     *             util_ascii_plotter                 *
*     *                                                *
*     **************************************************
      subroutine util_ascii_plotter(mark,lunit,n,x,y,Symbol,Title,
     >                              Xlabel,Ylabel)
      implicit none
      character*(*) mark
      integer lunit,n
      real*8 x(n),y(n)
      character*(*) Symbol,Title,Xlabel,Ylabel

      integer MaxRow,MaxCol
      !parameter (MaxRow=20,MaxCol=75)
      parameter (MaxRow=20,MaxCol=86)

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      character    Point(MaxRow,MaxCol)
      character*10 MyFormat,MyformatT,MyformatX,nstru,nstrd
      character*20 ttt
      integer i,j,l

      integer  ascii_xscale,ascii_yscale
      external ascii_xscale,ascii_yscale

      l = index(Title,' ')-1
      l = 43 - (l/2)
      write(MyformatT,'("(A,"I2, "x,A)")' ) l
      l = index(Xlabel,' ') -1
      l = 43 - (l/2)
      write(MyformatX,'("(A,"I2, "x,A)")' ) l

      write(nstru,'(E10.3)') YU
      write(nstrd,'(E10.3)') YD

      write(MyFormat,'("(A,"I2, "A1)")' ) MaxCol 
      Point = " "

      !*** set y-axis ***
      do j=1,MaxRow
         Point(j,min(maxcol,11+ascii_xscale(0.0d0))) = ":" ! y-axis
      end do

      !*** set x-axis ***
      do i=12,MaxCol
         j=min(maxrow,ascii_yscale(0.0d0))
         if (j.ge.1) Point(j,i) = "-"  ! x-axis
c         j=max(1,j)
c         Point(j,i) = "-"  ! x-axis
      end do

      !*** set ylabels ***
      do i=1,10
         Point(1, i) = nstru(i:i)
         Point(20,i) = nstrd(i:i)
      end do
      if (ascii_yscale(0.0d0).lt.20) then
         i=min(maxrow,ascii_yscale(0.0d0))
c         i=max(1,i)
         if (i.ge.1) then
            Point(i,1)  = " "
            Point(i,2)  = "0"
            Point(i,3)  = "."
            Point(i,4)  = "0"
            Point(i,5)  = "0"
            Point(i,6)  = "0"
            Point(i,7)  = "E"
            Point(i,8)  = "+"
            Point(i,9)  = "0"
            Point(i,10) = "0"
         end if
      end if

      !**** plot points ****
      do i=1,n
         Point(min(ascii_yscale(y(i)),maxrow),
     *        min(11+ascii_xscale(x(I)),maxcol)) = Symbol
      end do

      !**** write graph to lunit ****
      write(lunit,MyformatT) mark,Title
      write(lunit,'(A,A)') mark,Ylabel
      write(lunit,MyFormat) (mark,(Point(i,j), j=1,MaxCol), i=1,MaxRow)

      !*** set xlabels ***
      write(lunit,'(A,A,A,A)') mark,"           ",
     >               "|....................................|",
     >               "....................................|"
      write(lunit,900) mark,XL,(XL+XR)/2.0d0,XR
      write(lunit,MyformatX) mark,Xlabel
  900 format(A,5x,E10.3,27x,E10.3,27x,E10.3)

      return
      end

*     **************************************************
*     *                                                *
*     *             util_ascii_setwindow               *
*     *                                                *
*     **************************************************
      subroutine util_ascii_setwindow(xmin,xmax,ymin,ymax)
      implicit none
      real*8 xmin,xmax,ymin,ymax

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      XL=xmin
      XR=xmax
      YD=ymin
      YU=ymax
      return
      end

*     **************************************************
*     *                                                *
*     *             util_ascii_xscale                  *
*     *                                                *
*     **************************************************
      integer function ascii_xscale(x)
      implicit none
      real*8 x

      integer MaxRow,MaxCol
!      parameter (MaxRow=20,MaxCol=75)
      parameter (MaxRow=20,MaxCol=86)

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      ascii_xscale = NINT((x-XL)*(MaxCol-1)/(XR-XL)+1)
      return
      end

*     **************************************************
*     *                                                *
*     *             util_ascii_yscale                  *
*     *                                                *
*     **************************************************
      integer function ascii_yscale(y)
      implicit none
      real*8 y

      integer MaxRow,MaxCol
!      parameter (MaxRow=20,MaxCol=75)
      parameter (MaxRow=20,MaxCol=86)

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      ascii_yscale = NINT((y-YU)*(MaxRow-1)/(YD-YU)+1)
      return
      end
c $Id$
