      SubRoutine Get_Grid(Grid,nGrid,CooFrom,CooTo,nGridXYZ,
     &                    Volume,Out_Unit,iProc,Where,
     &                    nAtom,XYZ,What)
*
* $Id: get_grid.F,v 1.5 1998-05-06 08:08:21 d3e129 Exp $
*
      Implicit None
*
************************************************************************
#include "stdio.fh"
************************************************************************
*
*---- Arguments of a Get_Grid call
      Integer nGrid, nGridXYZ(3), Out_Unit, iProc, nAtom
      Double Precision Volume
      Double Precision Grid(3,nGrid),XYZ(3,nAtom)
      Double Precision CooFrom(3), CooTo(3)
      Character*8 Where,What
*---- Define local variables
      Double Precision Step(3),Spread(3),Angle(3),Scale,Frac_P
      Double Precision Dist,Thr
      Integer i,iGrid,iX,iY,iZ,Const_M(3),Const_P(3),iFast,iAtom
      Integer InRow,nTimes,iTimes,LeftNo,nLoop,Length,j,iOff,nRejected
      Logical DeBug
*
      Parameter (Thr = 1.d-8)
*
************************************************************************
*
      DeBug = .False.
*
      nRejected = 0
      Volume = 1.0d0
      Scale  = 1.8897265d0
      Do i = 1, 3
         Angle(i)  = 90.0d0
         If (nGridXYZ(i).gt.0) Then
            Spread(i) = CooTo(i)-CooFrom(i)
            Step(i)   = Spread(i)/dble( nGridXYZ(i) )
            Volume    = Volume * Step(i)
            If (CooFrom(i).ge.0.0d0) Then
               Const_M(i) = 0
               Const_P(i) = nGridXYZ(i)
            Else If (CooTo(i).le.0.0d0) Then
               Const_M(i) =-nGridXYZ(i)
               Const_P(i) = 0
            Else
               Frac_P     = (CooTo(i)-0.5d0*Step(i))/Spread(i)
               Const_P(i) = Int(Frac_P*nGridXYZ(i)) + 1
               Const_M(i) = Const_P(i) - nGridXYZ(i)
            End If
            Spread(i) = Spread(i)/Scale
         Else If (nGridXYZ(i).eq.0) Then
            Volume     =0.d0
            CooFrom(i) = (CooFrom(i)+CooTo(i))/2.d0
            CooTo(i)   = CooFrom(i)
            Spread(i)  = 0.d0
            Const_M(i) = 0
            Const_P(i) = 0
         Else
            Call ErrQuit('Get_Grid: Negative no_of_spacings', 0)
         End If
      End Do
      If (iProc.eq.0 . and . Where.ne.'NUCLEI') Then
         Write(Out_Unit,'(3(1x,f10.5),3(3x,f5.2))')
     &        (Spread(i),i=1,3),(Angle(i),i=1,3)
         Write(Out_Unit,'(3(2x,i4))')(nGridXYZ(i),i=1,3)
         iFast = 1
         Write(Out_Unit,'(5x,i1,6(2x,i4))')
     &         iFast,(Const_M(i),Const_P(i),i=1,3)
      End If
*
      If (What.eq.'CHARGES') GoTo 2345   ! ***** temporary *****
      iGrid = 0
      Do iZ = 0, nGridXYZ(3)
         Do iY = 0, nGridXYZ(2)
            Do iX = 0, nGridXYZ(1)
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
      If (iGrid+nRejected.ne.nGrid)
     &   Call ErrQuit('Get_Grid: iGrid .ne. nGrid', 0)
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
         InRow  = 10
         nTimes = nGrid/InRow
         LeftNo = nGrid - nTimes*InRow
         iOff   = 1
         nLoop  = nTimes
         If (LeftNo.gt.0) nLoop = nTimes + 1
         Write(LuOut,*)' Grid Points:'
         Do iTimes = 1, nLoop
            Length = InRow
            If (iTimes.gt.nTimes) Length = LeftNo
            Do j = 1, 3
               Write(LuOut,'(20f8.2)')(Grid(j,iOff+i),i=0,Length-1)
            End Do
            Write(LuOut,*)
            iOff = iOff + InRow
         End Do
      End If
*
      Return
      End
