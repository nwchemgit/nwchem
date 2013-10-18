      SubRoutine Get_Grid(Grid,nGrid,CooFrom,CooTo,nGridXYZ,
     &                    Volume,Out_Unit,iProc,Where,
     &                    nAtom,XYZ,What,lgaussian,charge)
*
* $Id$
*
      Implicit None
#include "errquit.fh"
*
************************************************************************
#include "stdio.fh"
#include "util.fh"
************************************************************************
*
*---- Arguments of a Get_Grid call
      Integer nGrid, nGridXYZ(3), Out_Unit, iProc, nAtom
      Double Precision Volume
      Double Precision Grid(3,nGrid),XYZ(3,nAtom)
      Double Precision CooFrom(3), CooTo(3),charge(*)
      Character*8 Where,What
*---- Define local variables
      Double Precision Step(3),Spread(3),Angle(3)
      Double Precision Dist,Thr
      Integer i,iGrid,iX,iY,iZ,Const_M(3),Const_P(3),iAtom
      Integer InRow,nTimes,iTimes,LeftNo,nLoop,Length,j,iOff,nRejected
      Logical DeBug,lgaussian
*
      Parameter (Thr = 1.d-8)
*
************************************************************************
*
      DeBug = util_print('dplot', print_high)
*
      nRejected = 0
      call  dplot_defgrid(volume,angle,
     ,     spread,step,const_m,const_p,
     ,     cooto,coofrom,ngridxyz)
*
      If (What.eq.'CHARGES') GoTo 2345   ! ***** temporary *****
      if(lgaussian) then
      iGrid = 0
      Do ix = 0, nGridXYZ(1)
         Do iY = 0, nGridXYZ(2)
            Do iz = 0, nGridXYZ(3)
*
               iGrid = iGrid + 1
               Grid(1,iGrid) = CooFrom(1) + dble(iX)*Step(1)
               Grid(2,iGrid) = CooFrom(2) + dble(iY)*Step(2)
               Grid(3,iGrid) = CooFrom(3) + dble(iZ)*Step(3)
               If (What.eq.'CHARGES') Then
                  Do iAtom = 1, nAtom
                     Dist = 0.d0
                     Do i = 1, 3
                        Dist = Dist + ( Grid(i,iGrid)-XYZ(i,iAtom) )**2
                     End Do
                     Dist = Sqrt(Dist)
                     If (Dist.le.Thr) Then
                        If (iProc.eq.0) Write(*,*)' Grid point No',
     &                     iGrid,' was rejected.'
                        iGrid = iGrid - 1
                        nRejected = nRejected + 1
                        GoTo 12340
                     End If
                  End Do
12340             Continue
               End If
*
            End Do
         End Do
      End Do
      else
      iGrid = 0
      Do iz = 0, nGridXYZ(3)
         Do iY = 0, nGridXYZ(2)
            Do ix = 0, nGridXYZ(1)
*
               iGrid = iGrid + 1
               Grid(1,iGrid) = CooFrom(1) + dble(iX)*Step(1)
               Grid(2,iGrid) = CooFrom(2) + dble(iY)*Step(2)
               Grid(3,iGrid) = CooFrom(3) + dble(iZ)*Step(3)
               If (What.eq.'CHARGES') Then
                  Do iAtom = 1, nAtom
                     Dist = 0.d0
                     Do i = 1, 3
                        Dist = Dist + ( Grid(i,iGrid)-XYZ(i,iAtom) )**2
                     End Do
                     Dist = Sqrt(Dist)
                     If (Dist.le.Thr) Then
                        If (iProc.eq.0) Write(*,*)' Grid point No',
     &                     iGrid,' was rejected.'
                        iGrid = iGrid - 1
                        nRejected = nRejected + 1
                        GoTo 1234
                     End If
                  End Do
 1234             Continue
               End If
*
            End Do
         End Do
      End Do
      endif
      If (iGrid+nRejected.ne.nGrid)
     &   Call ErrQuit('Get_Grid: iGrid .ne. nGrid', 0, INPUT_ERR)
      GoTo 3456
*
*---- In the case of calculating the atomic charges the grid is chosen
*---- to consist of points on the surface of the specified cube.
*---- This will be modified in the future: the points will lie on the
*---- 'molecular surface'
 2345 Continue
      iGrid = 0
      Do iZ = 0, nGridXYZ(3), nGridXYZ(3)
         Do iY = 0, nGridXYZ(2)
            Do iX = 0, nGridXYZ(1)
                  iGrid = iGrid + 1
                  Grid(1,iGrid) = CooFrom(1) + dble(iX)*Step(1)
                  Grid(2,iGrid) = CooFrom(2) + dble(iY)*Step(2)
                  Grid(3,iGrid) = CooFrom(3) + dble(iZ)*Step(3)
                  Do iAtom = 1, nAtom
                     Dist = 0.d0
                     Do i = 1, 3
                        Dist = Dist + ( Grid(i,iGrid)-XYZ(i,iAtom) )**2
                     End Do
                     Dist = Sqrt(Dist)
                     If (Dist.le.Thr) Then
                        If (iProc.eq.0) Write(*,*)' Grid point No',
     &                     iGrid,' was rejected.'
                        iGrid = iGrid - 1
                        nRejected = nRejected + 1
                        GoTo 1235
                     End If
                  End Do
 1235             Continue
            End Do
         End Do
      End Do
      Do iZ = 1, nGridXYZ(3) - 1
         Do iY = 0, nGridXYZ(2)
            Do iX = 0, nGridXYZ(1)
               If (iX.eq.0 . or . iX.eq.nGridXYZ(1) . or .
     &             iY.eq.0 . or . iY.eq.nGridXYZ(2)) Then
                  iGrid = iGrid + 1
                  Grid(1,iGrid) = CooFrom(1) + dble(iX)*Step(1)
                  Grid(2,iGrid) = CooFrom(2) + dble(iY)*Step(2)
                  Grid(3,iGrid) = CooFrom(3) + dble(iZ)*Step(3)
                  Do iAtom = 1, nAtom
                     Dist = 0.d0
                     Do i = 1, 3
                        Dist = Dist + ( Grid(i,iGrid)-XYZ(i,iAtom) )**2
                     End Do
                     Dist = Sqrt(Dist)
                     If (Dist.le.Thr) Then
                        If (iProc.eq.0) Write(*,*)' Grid point No',
     &                     iGrid,' was rejected.'
                        iGrid = iGrid - 1
                        nRejected = nRejected + 1
                        GoTo 1236
                     End If
                  End Do
 1236             Continue
               End If
            End Do
         End Do
      End Do
*
 3456 Continue
      nGrid = iGrid
*     If (iProc.eq.0) Write(*,*)' nGrid = ', nGrid
*
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)' Grid Points:'
         Do iTimes = 1, ngrid
               Write(LuOut,'(3f10.6)')(Grid(j,itimes),j=1,3)
         End Do
      End If
*
      Return
      End
