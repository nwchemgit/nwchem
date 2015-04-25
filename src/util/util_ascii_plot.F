*     **************************************************
*     *                                                *
*     *             util_ascii_plotter                 *
*     *                                                *
*     **************************************************
      subroutine util_ascii_plotter(lunit,n,x,y,Symbol,Title,
     >                              Xlabel,Ylabel)
      implicit none
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
      write(MyFormatT,'("("I2, "x,A)")' ) l
      l = index(Xlabel,' ') -1
      l = 43 - (l/2)
      write(MyFormatX,'("("I2, "x,A)")' ) l

      write(nstru,'(E10.3)') YU
      write(nstrd,'(E10.3)') YD

      write(MyFormat,'("("I2, "A1)")' ) MaxCol 
      Point = " "

      !*** set y-axis ***
      do j=1,MaxRow
         Point(j,11+ascii_xscale(0.0d0) ) = ":" ! y-axis
      end do

      !*** set x-axis ***
      do i=12,MaxCol
         Point(ascii_yscale(0.0d0),i) = "-"  ! x-axis
      end do

      !*** set ylabels ***
      do i=1,10
         Point(1, i) = nstru(i:i)
         Point(20,i) = nstrd(i:i)
      end do
      if (ascii_yscale(0.0d0).lt.20) then
         Point(ascii_yscale(0.0d0),1)  = " "
         Point(ascii_yscale(0.0d0),2)  = "0"
         Point(ascii_yscale(0.0d0),3)  = "."
         Point(ascii_yscale(0.0d0),4)  = "0"
         Point(ascii_yscale(0.0d0),5)  = "0"
         Point(ascii_yscale(0.0d0),6)  = "0"
         Point(ascii_yscale(0.0d0),7)  = "E"
         Point(ascii_yscale(0.0d0),8)  = "+"
         Point(ascii_yscale(0.0d0),9)  = "0"
         Point(ascii_yscale(0.0d0),10) = "0"
      end if

      !**** plot points ****
      do i=1,n
        Point(ascii_yscale(y(i)),11+ascii_xscale(x(I))) = Symbol
      end do

      !**** write graph to lunit ****
      write(lunit,MyformatT) Title
      write(lunit,'(A)') Ylabel
      write(lunit,MyFormat) ((Point(i,j), j=1,MaxCol), i=1,MaxRow)

      !*** set xlabels ***
      write(lunit,*) "          ",
     >               "|....................................|",
     >               "....................................|"
      write(lunit,900) ,XL,(XL+XR)/2.0d0,XR
      write(lunit,MyformatX) Xlabel
  900 format(7x,f6.1,31x,f6.1,31x,f6.1)

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
      parameter (MaxRow=20,MaxCol=75)

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
      parameter (MaxRow=20,MaxCol=75)

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      ascii_yscale = NINT((y-YU)*(MaxRow-1)/(YD-YU)+1)
      return
      end
c $Id$
