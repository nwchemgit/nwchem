C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_denom -- Apply denominators to a global array of integrals
C
C REVISION
C     $Id$
C
      Subroutine RIMP2_Denom(G_Ints, BaseI, BaseJ, NVirI, NVirJ,
     $   OccEigI, OccEigJ, VirEigI, VirEigJ)
      Implicit NONE
      Integer G_Ints, BaseI, BaseJ, NVirI, NVirJ
      Double Precision OccEigI(*), OccEigJ(*),
     $   VirEigI(NVirI), VirEigJ(NVirJ)
C
C ARGUMENTS
C     G_Ints   GA handle for matrix to be denominatorified [INOUT]
C     BaseI    Real value of lowest I in matrix [IN]
C     BaseJ    Real value of lowest J in matrix [IN]
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
C     subset of occupied orbitals in either direction.  The "real"
C     occupied indices of the top left corner of the matrix are
C     (BaseI, BaseJ).  If these are not correct, the wrong occupied
C     eigenvalues may be used!
C
C     This routine takes a data parallel approach to applying the
C     denominators, each processor accessing directly its local
C     patch of G_Ints using GA routines.  Load balance is determined
C     solely by the distribution of G_Ints.
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer   H_Int,Ailo,AIhi,BJlo,BJhi,LD
#ifdef  BAD_GACCESS 
      integer l_local
#endif
C
C     First, find out what patch of the matrix we have.  First index
C     is combined index AI, second is BJ
C
      Call GA_Distribution( G_Ints, GA_NodeID(), AIlo, AIhi, BJlo, BJhi)
C     
C     GA_Access doesn't like it if we request a bogus
C     patch (as might be returned if this node has no data
C     on it).
C     
      If ( (AIhi-AIlo+1) * (BJHi-BJLo+1) .ne. 0) then
C        
#ifdef  BAD_GACCESS
      if(.not.ma_push_get(MT_DBL,(aihi-ailo+1)*(bjhi-bjlo+1),
     $  'scratch buff', l_local, h_int)) call
     $  errquit('rimp2denom: pushget failed',0,0)
      ld=aihi-ailo+1
         write(0,*) ' a ',ailo,aihi,aihi-ailo+1
         write(0,*) ' b ',bjlo,bjhi,bjhi-bjlo+1
         write(0,*) ' ld ',ld,(aihi-ailo+1)*(bjhi-bjlo+1)
      call ga_get(g_ints,ailo,aihi,bjlo,bjhi,dbl_mb(H_Int),ld)
#else
         Call GA_Access( G_Ints, AIlo, AIhi, BJlo, BJhi, H_Int, LD)
#endif
      call rimp2_denom_a(g_ints,
     ,     AIlo, AIhi, BJlo, BJhi, dbl_mb(H_Int), LD,
     ,     nviri,nvirj,basei,basej,
     ,     occeigi,occeigj,vireigi,vireigj)
C        
#ifdef  BAD_GACCESS 
      call ga_put(g_ints,ailo,aihi,bjlo,bjhi,
     $     dbl_mb(h_int),ld)
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('denom: popstack failed',0,0)
#else
         Call GA_Release_Update( G_Ints, AIlo, AIhi, BJlo, BJhi)
#endif
      EndIf
C     
      Return
      End
c
      subroutine rimp2_denom_a(g_ints,
     ,     AIlo, AIhi, BJlo, BJhi, H_Int, LD,
     ,     nviri,nvirj,basei,basej,
     ,     occeigi,occeigj,vireigi,vireigj)
      implicit none
#include "mafdecls.fh"
      integer g_ints,AIlo, AIhi, BJlo, BJhi,LD 
      integer nviri,nvirj
      integer basei,basej
      double precision H_int(*),occeigi(*),occeigj(*),
     .     vireigi(*),vireigj(*)
c
      integer alo,ahi,blo,bhi
      integer ilo,ihi,jlo,jhi
      integer i,j,a,b,offset,bbot,btop,abot,atop
      double precision denij,den
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
      ILo = Int( (AILo-1) / NVirI ) + 1
      IHi = Int( (AIHi-1) / NVirI ) + 1
C     
      JLo = Int( (BJLo-1) / NVirJ ) + 1
      JHi = Int( (BJHi-1) / NVirJ ) + 1
C     
      Offset = 0
      Do J = JLo, JHi
C     
C     First and last J might not cover the whole range
C     
         BBot = 1
         BTop = NVirJ
         If ( J .eq. JLo) BBot = BLo
         If ( J .eq. JHi) BTop = BHi
C     
         Do B = BBot, BTop
C     
C     These two loops walk down a column of our patch
C     
            Do I = ILo, IHi
C     
C     Form half of the denominator
C     
               DenIJ = OccEigI(BaseI+I-1) + OccEigJ(BaseJ+J-1)
C     
C     As with B, A doesn't always cover the whole
C     range
C     
               ABot = 1
               ATop = NVirI
               If ( I .eq. ILo) ABot = ALo
               If ( I .eq. IHi) ATop = AHi
C     
               Do A = ABot, ATop
                  Offset = Offset + 1
C     
                  Den = DenIJ - VirEigI( A ) - VirEigJ( B )
                  H_Int( Offset )
     $                 =  H_int( Offset ) / Den
C     
               EndDo            ! A
            EndDo               ! I
C     
C     Jump to the proper offset for the head of the next
C     column.
C     
            Offset = Offset + LD - (AIhi - AIlo + 1)
C     
         EndDo                  ! B
      EndDo                     ! J
      return
      end
