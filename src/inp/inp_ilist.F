C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     inp_ilist -- read an integer list, possibly including ranges
C
C REVISION
C     $Id: inp_ilist.F,v 1.2 1995-11-03 20:42:47 d3g681 Exp $
C
C SYNOPSIS
      Logical Function Inp_IList( MaxList, List, N)
      Implicit NONE
      Integer MaxList           ! [in]    Size of List
      Integer List( MaxList )   ! [inout] Contents of list
      Integer N                 ! [out]   Number of elements in list
C
C DESCRIPTION
C     Reads the remainder of the line for a list of integers and puts
c     the results in List.  Ranges of integers may be compactly using
c     F90-style triplet notation. The number of elements set from the
c     input is returned in N. 
C
C     Inp_IList returns .TRUE. if the input is a valid integer list, 
C     and .FALSE. otherwise, also setting an appropriate error message.
C     If N > MaxList, it indicates that there is too much data on
C     the line to fit in List.
C
C INCLUDED FILES
#include "inpP.fh"
C
C EXTERNAL ROUTINES
      Integer  Inp_Cur_Field, Inp_N_Field
      Logical  Inp_IRange
      External Inp_Cur_Field, Inp_N_Field, Inp_IRange
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Field, First, Last, Stride, I, Start
C
      N = 0
C
C     Fields are numbered from 0..n-1
C
      Start = JRec
      Do Field = JRec, Jump-1
         If ( Inp_IRange(First, Last, Stride) ) then
C
C           Expand the triplet.  Notice that if there are too many
C           elements for the List, we count them but do not overwrite
C           the array bounds.
C
            Do I = First, Last, Stride
               N = N + 1
               If ( N .le. MaxList) List(N) = I
               If ( N .eq. MaxList+1) IErrPos = IStrt(JRec)
            EndDo
         Else
C
C           Inp_IRange will set IErrPos and ErrMsg appropriately.
C           We need to shove the pointer back to the beginning of
C           the whole list though.
C
            JRec = Start
            Inp_IList = .FALSE.
            Return
         EndIf
      EndDo
C
C     It is possible we have exceeded the array capacity.
C
      If ( N .le. MaxList) then
         Inp_IList = .TRUE.
      Else
C
C        IErrPos is set at the field which causes the overwrite
C
         ErrMsg = 'insufficient space for integer list'
         JRec = Start
         Inp_IList = .FALSE.
      EndIf
      Return
      End
