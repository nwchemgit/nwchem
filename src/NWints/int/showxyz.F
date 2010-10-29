      program showxyz
C$Id$
      implicit none
      integer Lvalm,Lval
      integer nxyz(3)
      integer icrng, ic
c
      Lvalm = 5
c
      call defnxyz(Lvalm)
c
      do 00100 Lval = 0,Lvalm*4
        icrng = ((Lval+1)*(Lval+2))/2
        do 00200 ic = 1,icrng
          call getNxyz(lval,ic,nxyz)
          call prntnxyz(lval,ic,nxyz)
00200   continue
        write(6,*)'  '
00100 continue
      end
      subroutine prntnxyz(lval,ic,nxyz)
      implicit none
      integer lval, ic, nxyz(3)
      character*30 cxyz
      character*1 pxyz(3)
      integer i,j, icnt
c
      pxyz(1) = 'x'
      pxyz(2) = 'y'
      pxyz(3) = 'z'
      cxyz=' '
      icnt = 1
      do 00100 i=1,3
        do 00200 j=1,nxyz(i)
          cxyz(icnt:icnt) = pxyz(i)
          icnt = icnt + 1
00200   continue
00100 continue
      write(6,*)'<',lval,',',ic,'>= ',cxyz
      end
