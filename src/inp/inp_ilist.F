C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     inp_ilist -- read an integer list, possibly including ranges
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Inp_IList( MaxList, List, N)
      Implicit NONE
      Integer MaxList           ! [in]    Size of List
      Integer List( MaxList )   ! [inout] Contents of list
      Integer N                 ! [out]   Number of elements in list
C
C DESCRIPTION
C     Reads the line for a list of integers and puts
c     the results in List.  Ranges of integers may be compactly using
c     F90-style triplet notation. The number of elements set from the
c     input is returned in N. 
c
c     It USED to read the entire remainder of the input line. BUT now
c     reads until it finds something that is not an integer range.
c
c     IF (the remainder of the input line was read as a valid list)
c     .  inp_list returns .true. with N as the number of elements read
c     ELSE 
c     .  inp_list returns .false. with N as the number of elements 
c     .  and also sets the INP internal error message appropriately
c     .  which can be cleared with inp_clear_err
c
c     N may be returned in either case as zero.
c
C     If N > MaxList, it indicates that there is too much data on
C     the line to fit in List. NOTE that inp_ilist now returns
c     true in this scenario so the value of N must be checked 
c     for this condition.  A recovery action would be to reposition
c     the input to the beginning of the list, allocate a sufficiently
c     large array and then call inp_ilist again.  The input position may
c     be recorded with inp_cur_field() and reset with inp_set_field().
c     The INP error message is also set in this circumstance.
C
C INCLUDED FILES
c     none.
C
C EXTERNAL ROUTINES
      Integer  Inp_Cur_Field, Inp_N_Field
      Logical  Inp_IRange
      External Inp_Cur_Field, Inp_N_Field, Inp_IRange
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer First, Last, Stride, I
C
      N = 0
C
c     Read as many fields as possible as integer ranges
c
 10   If ( Inp_IRange(First, Last, Stride) ) then
C
C     Expand the triplet.  If there are too many elements for the List, 
C     we count them but do not overwrite the array bounds.
C
         Do I = First, Last, Stride
            N = N + 1
            If ( N .le. MaxList) List(N) = I
            If ( N .eq. MaxList+1) call inp_mark_err
     $          ('insufficient space for integer list')
         EndDo
         goto 10
      EndIf
C
C     OK if we are we are the end of the input line 
C
      inp_ilist = (inp_cur_field() .eq. inp_n_field()) 
C
      End

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c name
c     inp_ilist_size -- returns number of elements in list
c
c revision
c     $id: inp_ilist.f,v 1.3 1997/03/05 23:20:31 d3g681 exp $
c
c synopsis
      logical function inp_ilist_size(n)
      implicit none
      integer n                 ! [out]   number of elements in list

c
c description
c     reads the line of a list of integers and 
c     and returns total number of elements
c     designed to be use a helper routine for inp_ilist
c     on exit from this routine the input is repositioned back
c     to the original value
c
c external routines
      integer  inp_cur_field, inp_n_field
      external inp_cur_field, inp_n_field
      logical  inp_irange
      external inp_irange
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c local variables
      integer first, last, stride, i
      integer nc
c
c     
      nc = inp_cur_field()
      n = 0
c
c     read as many fields as possible as integer ranges
c
 10   if ( Inp_IRange(first, last, stride) ) then
c
c     expand the triplet.  if there are too many elements for the list, 
c     we count them but do not overwrite the array bounds.
c
         do i = first, last, stride
            n = n + 1
         enddo
         goto 10
      endif
c
c
      inp_ilist_size = (inp_cur_field() .eq. inp_n_field()) 
      call inp_set_field(nc)
      end
