C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_stats -- compute statistics on t amplitudes
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Stats(LowOnly, ILo, IHi, NVirI, JLo, JHi,
     $     NVirJ, G_T,
     $     TopN, TopN_Ind, TopN_Val)
      Implicit NONE
      Logical LowOnly
      Integer ILo
      Integer IHi
      Integer NVirI
      Integer JLo
      Integer JHi
      Integer NVirJ
      Integer G_T
      Integer TopN
      Integer TopN_Ind(4, TopN)
      Double Precision TopN_Val(TopN)
C
C DESCRIPTION
C EXTERNAL ROUTINES
      Integer  IDAMin
      External IDAMin
C
C INCLUDED FILES
#include "mafdecls.fh"
#include "global.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer MinPtr
      Double Precision MinVal
      Integer LD, I_T
      Integer AILo, AIHi, BJLo, BJHi
      Integer ALo, AHi, BLo, BHi, OffSet, ABot, ATop, BBot, BTop
      Integer MyILo, MyIHi, MyJLo, MyJHi
      Integer I, J, A, B
C
      If ( TopN .lt. 1) Return
C
C     Find out what smallest value in current TopN list is.
C     This is what an element exceed to displace a current TopN value.
C
      MinPtr = IDAMin(TopN, TopN_Val, 1)
      MinVal = TopN_Val(MinPtr)
C
C     Get access to local region of t
C
      Call GA_Distribution( G_T, GA_NodeID(), AIlo, AIhi, BJlo, BJhi)
C     
C     GA_Access doesn't like it if we request a bogus
C     patch (as might be returned if this node has no data
C     on it).
C     
      If ( (AIhi-AIlo+1) * (BJHi-BJLo+1) .ne. 0) then
C        
         Call GA_Access( G_T, AIlo, AIhi, BJlo, BJhi, I_T, LD)
C        
C        Figure out what these correspond to in more 
C        basic terms.  Note that everything here is relative
C        to the base I and J for this matrix, which this routine
C        doesn't really care about.
C        
         ALo = Mod( AILo, NVirI )
         If ( ALo .eq. 0) ALo = NVirI
         AHi = Mod( AIHi, NVirI )
         If ( AHi .eq. 0) AHi = NVirI
C        
         BLo = Mod( BJLo, NVirJ )
         If ( BLo .eq. 0) BLo = NVirJ
         BHi = Mod( BJHi, NVirJ )
         If ( BHi .eq. 0) BHi = NVirJ
C        
         MyILo = Int( (AILo-1) / NVirI ) + 1
         MyIHi = Int( (AIHi-1) / NVirI ) + 1
C        
         MyJLo = Int( (BJLo-1) / NVirJ ) + 1
         MyJHi = Int( (BJHi-1) / NVirJ ) + 1
C        
         Offset = -1
         Do J = MyJLo, MyJHi
C           
C           First and last J might not cover the whole range
C           
            BBot = 1
            BTop = NVirJ
            If ( J .eq. MyJLo) BBot = BLo
            If ( J .eq. MyJHi) BTop = BHi
C           
            Do B = BBot, BTop
C              
C              These two loops walk down a column of our patch
C
               Do I = MyILo, MyIHi
C                 
C                 As with B, A doesn't always cover the whole range
C                 
                  ABot = 1
                  ATop = NVirI
                  If ( I .eq. MyILo) ABot = ALo
                  If ( I .eq. MyIHi) ATop = AHi
C
C                 Do this section only if necessary
C
                  If ( .NOT. LowOnly .OR. (LowOnly .AND. J .le. I)) then
C                 
                     Do A = ABot, ATop
                        Offset = Offset + 1
C                    
C                       Does this value displace the smallest TopN?
C
                        If ( Abs( Dbl_MB( I_T + Offset ) ) .gt.
     $                       MinVal ) then
C
C                          Record value & indices
C
                           TopN_Val(MinPtr) = Dbl_MB(I_T + Offset)
                           TopN_Ind(1, MinPtr) = I
                           TopN_Ind(2, MinPtr) = J
                           TopN_Ind(3, MinPtr) = A
                           TopN_Ind(4, MinPtr) = B
C
C                          Find new minimum
C
                           MinPtr = IDAMin(TopN, TopN_Val, 1)
                           MinVal = TopN_Val(MinPtr)

                        EndIf
C                    
                     EndDo      ! A
                  Else
                     OffSet = Offset + ATop - ABot + 1
                  EndIf
               EndDo            ! I
C
C              Jump to the proper offset for the head of the next
C              column.
C
               Offset = Offset + LD - (AIhi - AIlo + 1)
C              
            EndDo               ! B
         EndDo                  ! J
C        
         Call GA_Release( G_T, AIlo, AIhi, BJlo, BJhi)
      EndIf
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_stats_merge -- merge statistics on t amplitudes
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Stats_Merge( TopN, TopN_Ind, TopN_Val)
      Implicit NONE
      Integer TopN
      Integer TopN_Ind(4, TopN)
      Double precision TopN_Val(TopN)
C
C DESCRIPTION
C     Merge statistics across all processors
C
C EXTERNAL ROUTINES
      Integer  IDAMin
      External IDAMin
C
C INCLUDED FILES
#include "global.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Info
C
      If ( TopN .lt. 1) Return
C
C     First, sort our local idea of TopN to make merging easier
C
      Call RIMP2_t_Sort('D', TopN, TopN_Ind, TopN_Val)
C
      If ( GA_NNodes() .gt. 1) then
         Write (LuOut, *) '!! WARNING !! t statistics not correct ',
     $        'for more than one processor !! WARNING !!'
      EndIf
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_stats_print -- print statistics on t amplitudes
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Stats_Print( TopN, TopN_Ind, TopN_Val,
     $     NCorI, NCorJ, NOccI, NOccJ, SpinLabel)
      Implicit NONE
      Integer TopN
      Integer TopN_Ind(4, TopN)
      Double precision TopN_Val(TopN)
      Integer NCorI
      Integer NCorJ
      Integer NOccI
      Integer NOccJ
      Character*(*) SpinLabel
C
C DESCRIPTION
C     Print t amplitude statistics
C
C EXTERNAL ROUTINES
      Integer  IDAMin
      External IDAMin
C
C INCLUDED FILES
#include "global.fh"
#include "stdio.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I
C
      If (TopN .lt. 1) Return
C
C     Headline and matching data format.  Note: SpinLabel assumed to hold
C     'Spin adapted' or something shorter.
C     
      Write (LuOut, 9000)
 9000 Format(/4X, 'Amplitude', 7X, 'I', 4X, 'J', 4X, 'A', 4X, 'B', 3X,
     $     'Spin')
 9010 Format(1X, E15.6, 4I5, 3X, A)
C
      Do I = 1, TopN
C
C        Make sure this is a real value
C
         If ( TopN_Ind(1, I) .ne. 0 .AND. TopN_Ind(2, I) .ne. 0 .AND.
     $        TopN_Ind(3, I) .ne. 0 .AND. TopN_Ind(4, I) .ne. 0 ) then
C
            Write (LuOut, 9010) TopN_Val(I),
     $           TopN_Ind(1, I) + NCorI,
     $           TopN_Ind(2, I) + NCorJ,
     $           TopN_Ind(3, I) + NCorI + NOccI,
     $           TopN_Ind(4, I) + NCorJ + NOccJ, SpinLabel
         EndIf
      EndDo
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_stats_init -- initialize statistics on t amplitudes 
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Stats_Init( TopN, H_TopN_Ind, H_TopN_Val)
      Implicit NONE
#include "errquit.fh"
      Integer TopN ![in] Number of largest amplitudes desired
      Integer H_TopN_Ind ![out] MA handle of array for top N indices
      Integer H_TopN_Val ![out] MA handle of array for top N values
C
C DESCRIPTION
C     Convenience routine to allocate and initialize arrays required
C     for accumulation of statistics on the t amplitudes.
C
C     Memory is allocated from the heap. See also RIMP2_t_Stats_Term.
C
C INCLUDED FILES
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I_TopN_Ind, I_TopN_Val
C
      If (TopN .lt. 1) Return
C
      If ( .NOT. MA_Alloc_Get(MT_Int, TopN*4,
     $     'Indices of top N amplitudes', H_TopN_Ind, I_TopN_Ind) )
     $     Call ErrQuit(
     $     'RIMP2_t_Stats_Init: can''t allocate local array', TopN*4,
     &       MA_ERR)
      Call IFill(TopN*4, 0, Int_MB(I_TopN_Ind), 1)
C
      If ( .NOT. MA_Alloc_Get(MT_Dbl, TopN,
     $     'Values of top N amplitudes', H_TopN_Val, I_TopN_Val) )
     $     Call ErrQuit(
     $     'RIMP2_t_Stats_Init: can''t allocate local array', TopN,
     &       MA_ERR)
      Call DFill(TopN, 0.0d0, Dbl_MB(I_TopN_Val), 1)
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_stats_term -- terminate statistics on t amplitudes 
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Stats_Term( TopN, H_TopN_Ind, H_TopN_Val)
      Implicit NONE
#include "errquit.fh"
      Integer TopN       ![in]  Number of largest amplitudes desired
      Integer H_TopN_Ind ![out] MA handle of array for top N indices
      Integer H_TopN_Val ![out] MA handle of array for top N values
C
C DESCRIPTION
C     Convenience routine to free arrays created by RIMP2_t_Stats_Init
C
C INCLUDED FILES
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      If ( TopN .lt. 1) Return
C
      If ( .NOT. MA_Free_Heap(H_TopN_Val) ) Call ErrQuit(
     $     'RIMP2_t_Stats_Term: can''t free local array', 1, MA_ERR)
C
      If ( .NOT. MA_Free_Heap(H_TopN_Ind) ) Call ErrQuit(
     $     'RIMP2_t_Stats_Term: can''t free local array', 2, MA_ERR)
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_sort -- specialized sort for t amplitudes
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_Sort( ID, N, Indices, D)
      Implicit NONE
#include "errquit.fh"
      Character*1 ID ![in] 'I' to sort in increasing order, 'D' for decreasing
      Integer N      ![in]
      Integer Indices(4, N) ![inout] t amplitude indices
      Double precision D(N) ![inout] t amplitude values
C
C DESCRIPTION
C     Sort a list of amplitudes & indices by ascending or decending 
C     absolute value.
C
C     Derived from LAPACK's DLASRT routine.
C
C EXTERNAL ROUTINES
      LOGICAL            LSAME
      EXTERNAL           LSAME
      EXTERNAL           XERBLA
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Info
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
C
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         Call ErrQuit('RIMP2_t_Sort: Error in argument', Info,
     &       INPUT_ERR)
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( Abs(D( J )).GT.Abs(D( J-1 )) ) THEN
                     DMNMX = D( J )
                     Call RIMP2_t_Swap( D(j), D(j-1), Indices(1,j),
     $                    Indices(1, j-1) )
c$$$                     D( J ) = D( J-1 )
c$$$                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( Abs(D( J )).LT.Abs(D( J-1 )) ) THEN
                     DMNMX = D( J )
                     Call RIMP2_t_Swap( D(j), D(j-1), Indices(1,j),
     $                    Indices(1, j-1) )
c$$$                     D( J ) = D( J-1 )
c$$$                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( Abs(D1).LT.Abs(D2) ) THEN
            IF( Abs(D3).LT.Abs(D1) ) THEN
               DMNMX = D1
            ELSE IF( Abs(D3).LT.Abs(D2) ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( Abs(D3).LT.Abs(D2) ) THEN
               DMNMX = D2
            ELSE IF( Abs(D3).LT.Abs(D1) ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( Abs( D( J ) ).LT.Abs(DMNMX) )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( Abs( D( I ) ).GT.Abs(DMNMX) )
     $         GO TO 80
            IF( I.LT.J ) THEN
               Call RIMP2_t_Swap( D(i), D(j), Indices(1,i),
     $              Indices(1, j) )
c$$$               TMP = D( I )
c$$$               D( I ) = D( J )
c$$$               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( Abs( D( J ) ).GT.Abs(DMNMX) )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( Abs( D( I ) ).LT.Abs(DMNMX) )
     $         GO TO 110
            IF( I.LT.J ) THEN
               Call RIMP2_t_Swap( D(i), D(j), Indices(1,i),
     $              Indices(1, j) )
c$$$               TMP = D( I )
c$$$               D( I ) = D( J )
c$$$               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_t_swap -- specialized swap function for t amplitude sort
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_t_swap( D1, D2, I1, I2)
      Implicit NONE
      Double Precision D1
      Double Precision D2
      Integer I1(4)
      Integer I2(4)
C
C DESCRIPTION
C     Swap amplitude values and indices for sorting.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIBLES
      Double Precision DTmp
      Integer ITmp, i
C
      DTmp = D1
      D1 = D2
      D2 = DTmp
C
      Do i = 1, 4
         ITmp = I1(i)
         I1(i) = I2(i)
         I2(i) = ITmp
      EndDo
C     
      Return
      End
