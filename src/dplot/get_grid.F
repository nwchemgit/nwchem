      SubRoutine Get_Grid(Grid,nGrid,CooFrom,CooTo,nGridXYZ,
     &                    Volume,Out_Unit,iProc,Where)
      Implicit None
*
************************************************************************
#include "stdio.fh"
************************************************************************
*
*---- Arguments of a Get_Grid call
      Integer nGrid, nGridXYZ(3), Out_Unit, iProc
      Double Precision Volume
      Double Precision Grid(3,nGrid)
      Double Precision CooFrom(3), CooTo(3)
      Character*8 Where
*---- Define local variables
      Double Precision Step(3),Spread(3),Angle(3),Scale,Frac_P
      Integer i,iGrid,iX,iY,iZ,Const_M(3),Const_P(3),iFast
      integer InRow, nTimes, iTimes, LeftNo, nLoop, Length, j, iOff
      Logical DeBug
*
************************************************************************
*
      DeBug = .False.
*
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
      iGrid = 0
      Do iZ = 0, nGridXYZ(3)
         Do iY = 0, nGridXYZ(2)
            Do iX = 0, nGridXYZ(1)
               iGrid = iGrid + 1
               Grid(1,iGrid) = CooFrom(1) + dble(iX)*Step(1)
               Grid(2,iGrid) = CooFrom(2) + dble(iY)*Step(2)
               Grid(3,iGrid) = CooFrom(3) + dble(iZ)*Step(3)
            End Do
         End Do
      End Do
      If (iGrid.ne.nGrid)
     &   Call ErrQuit('Get_Grid: iGrid .ne. nGrid', 0)
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
