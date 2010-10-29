C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     panel_diagonal -- Locate the piece of the diagonal running through
C     the panel.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical function Panel_Diagonal( Ilo, Ihi, Jlo, Jhi, Dlo, Dhi)
      Implicit NONE
      Integer Ilo, Ihi, Jlo, Jhi, Dlo, Dhi
C
C ARGUMENTS
C     The panel is described by (Ilo:Ihi, Jlo:Jhi) [IN], and the
C     diagonal contained within the panel is (Dlo:Dhi) [OUT].
C     If the panel does not contain any elements of the diagonal
C     the function returns .FALSE.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer UR, LL
C
C     Check that the ranges overlap.  These two will be the same sign
C     if the do NOT overlap
C
      UR = Ihi - Jlo
      LL = Ilo - Jhi
C
      If ( UR*LL .gt. 0) then
         Panel_Diagonal = .FALSE.
         Return
      Else
         Panel_Diagonal = .TRUE.
      EndIf
C
C     Now we must locate the overlap of the ranges, since that must
C     be where the diagonal occurs.
C
      Dlo = Max( Ilo, Jlo)
      Dhi = Min( Ihi, Jhi)
C
      Return
      End
