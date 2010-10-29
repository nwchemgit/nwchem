C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_denom2 -- Apply denominators to a global array of integrals
C
C REVISION
C     $Id$
C
      Subroutine RIMp2_denom2(IJSym, IJTri, G_T, ILo, IHi, JLo, JHi,
     $   NVirI, NVirJ,
     $   OccEigI, OccEigJ, VirEigI, VirEigJ)
      Implicit NONE
      Logical IJSym, IJTri
      Integer G_T, ILo, IHi, JLo, JHI, NVirI, NVirJ
      Double Precision OccEigI(*), OccEigJ(*),
     $   VirEigI(NVirI), VirEigJ(NVirJ)
C
C ARGUMENTS
C     IJSym    Take advantage of IJ permutational symmetry [IN]
C     IJTri    IJ is triangularized in incoming data [IN]
C     G_T      GA handle for matrix to be denominatorified [INOUT]
C     NVirI    Number of virtuals in row blocks [IN]
C     NVirJ    Number of virtuals in column blocks [IN]
C     OccEigI  Occupied eigevalues for I spin case [IN]
C     OccEigJ  Occupied eigevalues for J spin case [IN]
C     VirEigI  Virtual eigevalues for I spin case [IN]
C     VirEigJ  Virtual eigevalues for J spin case [IN]
C
C DESCRIPTION
C     Applies denominators to a global array of integrals.  The GA
C     is organized as multiple blocks of virtual orbitals in each
C     dimension.  Each block corresponds to a different occupied
C     index.  The number of virtuals in each direction may be different
C     to accomodate mixed spin cases.  The matrix may include a
C     subset of occupied orbitals in either direction.  The occupied
C     eigenvalues OccEigI and OccEigJ should correspond to the subset
C     of orbitals represented in the matrix.  In other words, if the
C     upper left corner of the matrix is at (I,J) = (3,5) then the
C     first element of OccEigI should be the eigenvalue of orbital 3
C     and the first element of OccEigJ should be orbital 5.
C
C     This routine takes a data parallel approach to applying the
C     denominators, each processor accessing directly its local
C     patch of G_T using GA routines.  Load balance is determined
C     solely by the distribution of G_T.
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer AILo, AIHi, BJLo, BJHi, MyILo, MyIHi, MyJLo, MyJHi,
     $   ALo, AHi, BLo, BHi, Offset, LD, I, J, A, B, ABot, ATop,
     $   BBot, BTop, H_Int, TrueI, TrueJ
      Double precision DenIJ, Den
      Logical LowOnly
C
C     See if we should do only the lower triangle of the incoming data
C
      LowOnly = IJSym .AND. ILo .eq. JLo .AND. IHi. eq. JHi
     $   .AND. .not. IJTri
C
C     First, find out what patch of the matrix we have.  First index
C     is combined index AI, second is BJ
C
      Call GA_Distribution( G_T, GA_NodeID(), AIlo, AIhi, BJlo, BJhi)
C     
C     GA_Access doesn't like it if we request a bogus
C     patch (as might be returned if this node has no data
C     on it).
C     
      If ( (AIhi-AIlo+1) * (BJHi-BJLo+1) .ne. 0) then
C        
         Call GA_Access( G_T, AIlo, AIhi, BJlo, BJhi, H_Int, LD)
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
C                    Form half of the denominator
C                 
                     If ( IJTri ) then
                        Call Util_Detriangularize(J, TrueI, TrueJ)
                        DenIJ = OccEigI(TrueI) + OccEigJ(TrueJ)
                     Else
                        DenIJ = OccEigI(I) + OccEigJ(J)
                     EndIf
C                 
                     Do A = ABot, ATop
                        Offset = Offset + 1
C                    
                        Den = DenIJ - VirEigI( A ) - VirEigJ( B )
                        Dbl_MB( H_Int + Offset )
     $                     =  Dbl_MB( H_Int + Offset ) / Den
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
         Call GA_Release_Update( G_T, AIlo, AIhi, BJlo, BJhi)
      EndIf
C     
      Return
      End
