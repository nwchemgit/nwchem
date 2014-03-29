
      subroutine util_ascii_plotter(lunit,n,x,y,Symbol,Title)
      implicit none
      integer lunit,n
      real*8 x(n),y(n)
      character*(*) Symbol,Title

      integer MaxRow,MaxCol
      parameter (MaxRow=20,MaxCol=75)

      real*8  XL,XR,YD,YU
      common /ascii_pcommon/ XL,XR,YD,YU

      character    Point(MaxRow,MaxCol)
      character*10 MyFormat
      integer i,j

      integer  ascii_xscale,ascii_yscale
      external ascii_xscale,ascii_yscale

      write(MyFormat,'("("I2, "A1)")' ) MaxCol 
      Point = " "
      do j=1,MaxRow
         Point(j,ascii_xscale(0.0d0) ) = ":" ! y-axis
      end do
      do i=1,MaxCol
         Point(ascii_yscale(0.0d0),i) = "-"  ! x-axis
      end do

      do i=1,n
        Point(ascii_yscale(y(i)),ascii_xscale(x(I))) = Symbol
      end do
      write(lunit,'(A80 /)') Title
      write(lunit,MyFormat) ((Point(i,j), j=1,MaxCol), i=1,MaxRow)

      return
      end

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
