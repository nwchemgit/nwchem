C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     print_correl_map -- Print mapping from correlated orbitals to
C     the full SCF reference, along with occupations and eigenvalues.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Print_Correl_Map(LUN, M, N, Map, LDMap, MFmt,
     $   Occ, LDOcc, OFmt, Evl, LDEvl, EFmt, LinWid)
      Implicit NONE
#include "errquit.fh"
      Integer LUN, M, LDMap, LDOcc, LDEvl, LinWid
      Integer N(M), Map(LDMap, M)
      Double precision Occ(LDOcc, M), Evl(LdEvl, M)
      Character*(*) MFmt, OFmt, EFmt
C
C EXTERNAL ROUTINES
      Integer AtoI
      External AtoI
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C FIXED PARAMETERS
C
C     Allowed format specifiers
C
      Character*(*) FPFmt, IFmt
      Parameter (FPFmt = 'DdEeFfGg', IFmt = 'Ii')
C
C     Row labels for output
C
      Character*(*) CorLab, SCFLab, OccLab, EvlLab
      Parameter( CorLab = 'Cor. Orb. ' )
      Parameter( SCFLab = 'SCF Orb.  ' )
      Parameter( OccLab = 'Occup. Nr.' )
      Parameter( EvlLab = 'Eigenval. ' )
C
      Integer MaxStr
      Parameter (MaxStr = 255)
      Character*(MaxStr) MOutFmt, OOutFmt, EOutFmt
C
      Integer I, J, K, JJ, J2, LLab, NPerL
      Integer MWid, OWid, EWid, MxWid, MPad, OPad, EPad
C
C     *****************
C     * Check formats *
C     *****************
C
      If ( Index(IFmt, MFmt(1:1)) .eq. 0) Call ErrQuit(
     $   'Print_Correl_Map: Bad format spec for map', 0, INPUT_ERR)
      If ( Index(FPFmt, OFmt(1:1)) .eq. 0) Call ErrQuit(
     $   'Print_Correl_Map: Bad format spec for occupation', 0,
     &       INPUT_ERR)
      If ( Index(FPFmt, EFmt(1:1)) .eq. 0) Call ErrQuit(
     $   'Print_Correl_Map: Bad format spec for eigenvalues', 0,
     &       INPUT_ERR)
C
C     Get the width of each format spec, find the largest
C
      Call IndInt( MFmt, I, J)
      MWid = AtoI( MFmt(i:j) )
      Call IndInt( OFmt, I, J)
      OWid = AtoI( OFmt(i:j) )
      Call IndInt( EFmt, I, J)
      EWid = AtoI( EFmt(i:j) )
C
      MxWid = Max( MWid, OWid, EWid)
      MPad = MXWid - MWid
      OPad = MxWid - OWid
      EPad = MxWid - EWid
C
      If ( MWid .lt. 1 .OR. OWid .lt. 1 .OR. EWid .lt. 1) Call ErrQuit(
     $   'Print_Correl_Map: Bad format specification', 0, INPUT_ERR)
C
C     ********************************
C     * Determine space requirements *
C     ********************************
C     Width of the row label
C
      LLab = Max( Len(CorLab), Len(SCFLab), Len(OccLab), Len(EvlLab))
C
C     Determine now much we can fit on a line under the following
C     guidelines:
C        Column 1 blank
C        Column LinWid blank (to avoid extra newlines getting put in)
C        Elements separated by two spaces
C
      NPerL = Int( (LinWid-LLab-2) / (MxWid+2) )
C
C     Make sure we can fit something!!!
C
      If (NPerL .eq. 0) Call ErrQuit(
     $   'Print_Correl_Map: Can''t fit anything on a line', MxWid,
     &       INPUT_ERR)
C
C     *****************************
C     * Create the format strings *
C     *****************************
C     Orbital numbers
C
      If ( MPad+2 .ne. 0) then
         Write (MOutFmt, 9000) LLab, NPerL, MPad+2, MFmt
      Else
         Write (MOutFmt, 9010) LLab, NPerL, MFmt
      EndIf
C
C     Occupation numbers
C
      If ( OPad+2 .ne. 0) then
         Write (OOutFmt, 9000) LLab, NPerL, OPad+2, OFmt
      Else
         Write (OOutFmt, 9010) LLab, NPerL, OFmt
      EndIf
C
C     Eigenvalues
C
      If ( EPad+2 .ne. 0) then
         Write (EOutFmt, 9000) LLab, NPerL, EPad+2, EFmt
      Else
         Write (EOutFmt, 9010) LLab, NPerL, EFmt
      EndIf
C
C     When this prints out it looks like:
C          (1X, A<number>, <number>( <number>X, <format_spec>))
C     and since the F77 standard prohibits 0 repeat factors,
C     we need a second version in case the padding is zero (which
C     it will be for one of the lines).
C
 9000 Format( '(1X, A', I10, ',', I10,'(',I10, 'X,', A, '))')
 9010 Format( '(1X, A', I10, ',', I10,'(', A, '))')
C
      Do K = 1, M
         Do j = 1, N(K), NPerL
            J2 = Min( J+NPerL-1, N(K))
            Write (LUN, '(A)') ' '
            Write (LUN, MOutFmt) CorLab, (jj, jj = J, J2)
            Write (LUN, OOutFmt) OccLab, (Occ(jj, k), jj = J, J2)
            Write (LUN, EOutFmt) EvlLab, (Evl(jj, k), jj = J, J2)
            Write (LUN, MOutFmt) SCFLab, (Map(jj, k), jj = J, J2)
         EndDo
      EndDo
C
      Return
      End
