C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_e_final -- Complete evaluation of energy
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_E_Final(TopSpin, G_PairE, Contrib, E2)
      Implicit NONE
      Integer TopSpin ![in]
      Integer G_PairE(2*TopSpin-1) ![in]
      Double precision Contrib(TopSpin, TopSpin) ![out]
      Double precision E2 ![out]
C
C DESCRIPTION
C
C EXTERNAL ROUTINES
      Double precision Sum_Pairs
      External Sum_Pairs
C
C INCLUDE FILES
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Me, ISpin, JSpin
C
      Me = GA_NodeID()
C
C     Print pair energies if requested
C
      If ( Util_Print('pair energies', Print_Debug) ) then
         If ( Me .eq. 0 ) then
            Call Banner(LuOut, 'Pair Correlation Energies',
     $         '-', .TRUE., .TRUE., .FALSE.)
         EndIf
C
         Do ISpin = 1, TopSpin
            Do JSpin = ISpin, TopSpin
               Call GA_Print( G_PairE(ISpin+JSpin-1) )
            EndDo
         EndDo
      EndIf
C
C     Sum pair energies into total contributions for each spin case.
C     Note that ba spin case is identical to ab.  Up to this point,
C     we only carried ab terms, but now we make it explicit.  There is
C     redundant work here in calling Sum_Pairs for both ab and ba, but
C     it is small so we will not worry.
C     

      Do ISpin = 1, TopSpin
         Do JSpin = 1, TopSpin
            Contrib(ISpin, JSpin) =
     $         Sum_Pairs( G_PairE(ISpin+JSpin-1) ,
c$$$     $         ISpin .eq. JSpin)
     $         .False.)
         EndDo
      EndDo
C
C     Scale to account for spin integration factors
C
C     The contraction is E = 1/4 Sum(ijab) [(ia|jb) - (ibja)] t(ij,ab)
C     Non-zero spin cases are as follows:
C     i j a b
C     a a a a [(IA|JB) - (IB|JA)] t(IJ,AB)
C
C     b b b b [(ia|jb) - (ib|ja)] t(ij,ab)
C
C     a b a b [(IA|jb) -    0   ] t(Ij,Ab)
C     a b b a [   0    - (IB|ja)] t(Ij,aB) --> +(IA|jb) t(Ij,Ab)
C
C     b a b a [(ia|JB) -    0   ] t(iJ,aB)
C     b a a b [   0    - (ib|JA)] t(iJ,Ab) --> +(ia|JB) t(iJ,aB)
C
C     Thus the pure spin case gets a factor of 1/4, while the mixed
C     spins get a factor of 1/2.
C
      Call dscal(TopSpin*TopSpin, 0.5d0, Contrib, 1)
      Call dscal(TopSpin, 0.5d0, Contrib, TopSpin+1)
c$$$      If ( TopSpin .eq. 2) then
c$$$     $   Call dscal(TopSpin, 0.25d0, Contrib, TopSpin+1)
C
      If ( Util_Print('energy spin cases', Print_Debug) )
     $   Call DGEWR(TopSpin, TopSpin, Contrib, TopSpin, 0, 0,
     $   LuOut, 'E(2) spin contributions', 80, 'F20.15')
C
C     Finally, total the spin components
C
      E2 = 0.0d0
      Do JSpin = 1, TopSpin
         Do ISpin = 1, TopSpin
            E2 = E2 + Contrib(ISpin, JSpin)
         EndDo
      EndDo
C
      Return
      End
