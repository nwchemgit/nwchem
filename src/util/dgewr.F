      Subroutine DGeWr (N, M, A, LDA, N0, M0, LUN, Label, LinWid,
     $   EFmt)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id$
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     dgewr -- Write a general matrix nicely formatted 
C
C SYNOPSIS
      Implicit NONE
#include "errquit.fh"
      Integer N, M, LDA, N0, M0, LinWid, LUN
      Double precision A(LDA,*)
      Character*(*) Label, EFmt
C
C ARGUMENTS
C     N       Row dimension of A (input)
C     M       Column dimension of A (input)
C     A       Matrix to be symmetrized (modified)
C     LDA     Leading dimension of A (input)
C     N0      Offset for numering of rows (input)
C     M0      Offset for numbering of columns (input)
C     LUN     Unit number to write to (input)
C     Label   Text label for matrix (input)
C     LinWid  Maximum line width to use (input)
C     EFmt    (Fortran) format of individual elements (input)
C             Allowed formats are {Ew.d,Fw.d,Dw.d,Gw.d}
C
C DESCRIPTION
C     "Pretty prints" a general matrix to the specified unit, which
C     must be opened for formatted writing, and positioned as desired.
C
C     The given element format string is analyzed to determine the
C     width required by each element, and two spaces are allowed
C     between elements.  This width is used to determine the number of
C     elements that can be output on each line.
C
C     The matrix is printed headed by the given label, which is
C     surrounded by lines of '-'s above and below to set it off.
C     If the label is passed in blank, no heading is printed.
C
C     Non-zero row/column offsets allow displayed row/column numbers
C     to be different from the numbering of the matrix itself.
C
C LIMITATIONS
C     Supported formats are: [DEFG]w.d
C     Long labels are wrapped without regard for their content.
C     The dashed lines that are used to surround the label are at
C        most MAXSTR in length.
C     Internally generated format strings rely on the fact that
C        blanks are not significant.  If this changes, it will be
C        necessary to be very careful in constructing the parts of
C        the format strings that pertain to integers (Iw specifiers).
C
C WARNING
C     The string PosFmt is case-sensitive.  If this source is
C     forced to all upper- or lowercase without some amount of care,
C     the routine may only accept format specifiers in a specific
C     case.
C
C TO BE ADDED
C     If row dimension is 1, don't print row label
C
C INCLUDED FILES
#include "inp.fh"
C
C ROUTINES REQUIRED
C     IndInt  Locate the left-most integer in a string (FStrings)
C     AtoI    Convert string to an integer (FStrings)
C
      Integer AtoI
      External AtoI
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C FIXED PARAMETERS
C     PosFmt  Allowed format speficiers
C
      Character*(*) PosFmt
      Parameter (PosFmt = 'DdEeFfGg')
C
C LOCAL VARIABLES
C     I       General use
C     J       General use
C     J2      General use
C     ii      General use
C     JJ      general use
C     LLabNr  Length of the row/column labeling number
C     NPerL   Number of elements per line
C     IW      Width of element, extracted from EFmt
C     ID      Decimal places of element, extracted from EFmt
C     LElem   Larger of LLabNr and IW, used as actual element width
C     LLen    Label length, used to pretty print it
C     MaxStr  Maximum length allocated to internal strings
C     CLFmt   Format string for column label
C     DFmt    Format string for lines of data
C     Dash    Holds the lines used in pretty printing Label
C
      Integer i, j, J2, ii, jj, LLabNr, NPerL, IW, ID, LElem, LLen
      Integer MaxStr
      Parameter (MaxStr = 255)
      Character*(MaxStr) CLFmt, DFmt, Dash
C
C     *****************************
C     * Check matrix dimensioning *
C     *****************************
C
      If (N .lt. 0 .OR. M .lt. 0) Call ErrQuit(
     $   'DGEWR: matrix dimensions cannot be negative', Min(M, N) ,
     &       INPUT_ERR)
      If (LDA .lt. Max(1, N)) Call ErrQuit(
     $   'DGEWR: leading dimension too small', LDA, INPUT_ERR)
C
C     ************************
C     * Check element format *
C     ************************
C     Make sure the format specifier is reasonable
C
      If (Index(PosFmt, EFmt(1:1)) .eq. 0) Call ErrQuit(
     $   'DGEWR: Unknown format specifer', 0, INPUT_ERR)
C
C     Get the total width
C
      Call IndInt (EFmt, I, J)
      IW = AtoI(EFmt(i:j))
      If (IW .le. 0) Call ErrQuit(
     $   'DGEWR: line width must be positive', IW, INPUT_ERR)
C
C     Next should be a '.'
C
      If (EFmt(j+1:j+1) .ne. '.') Call ErrQuit(
     $   'DGEWR: bad format specifier', 0, INPUT_ERR)
C
C     Get the number of decimal places
C
      Call IndInt (EFmt(j+1:), I, J2)
      ID = AtoI(EFmt(j+i:j+j2))
      If (ID .le. 0) Call ErrQuit('DGEWR: bad format specifier', ID,
     &       INPUT_ERR)
C
C     *********************************
C     * See if we can get away easily *
C     *********************************
C
      If (N .eq. 0 .OR. M .eq. 0) Return
C
C     ********************************
C     * Determine space requirements *
C     ********************************
C     Width of the labeling row/column number
C     We have to account for the offsets and what they might do to the
C     signs of the labels (which require an extra space).
C
      I = Min( 1+N0, 1+M0 )
      J = Max( N+N0, M+M0 )
      II = Int( Log10( dble( Abs(I)) ) +1 )
      JJ = Int( Log10( dble( Abs(J)) ) +1 )
      LLabNr = Max( II, JJ)
      If (II .ge. JJ .AND. I .lt. 0) LLabNr = LLabNr + 1
C
C     Determine now much we can fit on a line under the following
C     guidelines:
C        Column 1 blank
C        Column LinWid blank (to avoid extra newlines getting put in)
C        Elements separated by two spaces
C        Column width is larger of the element and the label
C
      LElem = Max (IW, LLabNr)
      NPerL = Int( (LinWid-LLabNr-2) / (LElem+2) )
C
C     Make sure we can fit something!!!
C
      If (NPerL .eq. 0) Call ErrQuit(
     $   'DGEWR: Specified format will not fit in line width', 0,
     &       INPUT_ERR)
C
C     *****************************
C     * Create the format strings *
C     *****************************
C
      Write (CLFmt, 9000) 1+LLabNr, NPerL, LElem - LLabNr + 2, LLabNr
 9000 Format( '(/', I10,'X,', I10,'(',I10, 'X,I', I10, '))' )
C
C     Number of blanks between elements must be increased if
C     column labels are larger than the element width
C
      Write (DFmt, 9010) LLabNr, NPerL, LElem - IW + 2,
     $   EFmt(:Inp_StrLen(EFmt))
 9010 Format( '(1X,I', I10,',', I10,'(', I10,'X,',A,'))' )
C
C     *******************
C     * Write the label *
C     *******************
C
      LLen = Inp_StrLen(Label)
      If (LLen .ne. 0) then
         Do 100 I = 1, Min(LinWid-2, LLen, MaxStr)
            Dash(i:i) = '-'
 100     Continue
         Write (LUN, 9100) Dash(:Min(LinWid-2, LLen, MaxStr))
         Do 200 i = 1, LLen, LinWid-2
            Write (LUN, 9100) Label(i:Min( LLen, i+LinWid-3))
 200     Continue
         Write (LUN, 9100) Dash(:Min(LinWid-2, LLen, MaxStr))
      EndIf
 9100 Format (1X, A)
C
C     ********************
C     * Write the matrix *
C     ********************
C
      Do 1000 j = 1, M, NPerL
         J2 = Min( J+NPerL-1, M)
         If (M .gt. 1) Write (LUN, CLFmt) (jj+M0, jj = J, J2)
         Write (LUN, '(A)') ' '
         Do 1010 i = 1, N
            Write (LUN, DFmt) i+N0, (A(i,jj), jj = J, J2)
 1010    Continue
 1000 Continue
C
      Return
      End
