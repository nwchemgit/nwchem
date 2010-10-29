      Subroutine RIMP2_PStat_Init
C$Id$
      Implicit NONE
#include "errquit.fh"
C
#include "pstat.fh"
#include "rimp2_timers.fh"
C
      Logical Status
      Integer Accum(1)
C
      Status = .TRUE.
      Status = Status .AND. PStat_Allocate('Overall',
     $   PStat_QStat, 0, Accum, PS_Overall)
      Status = Status .AND. PStat_Allocate('3c integral xform.',
     $   PStat_QStat, 0, Accum, PS_xform)
      Status = Status .AND. PStat_Allocate('Energy',
     $   PStat_QStat, 0, Accum, PS_Energy)
      Status = Status .AND. PStat_Allocate('3c integral generation',
     $   PStat_QStat, 0, Accum, PS_int_gen)
      Status = Status .AND. PStat_Allocate('First 2/3 xform.',
     $   PStat_QStat, 0, Accum, PS_step12)
      Status = Status .AND. PStat_Allocate('Final 1/3 xform.',
     $   PStat_QStat, 0, Accum, PS_step3)
      Status = Status .AND. PStat_Allocate('Write 3c integrals',
     $   PStat_QStat, 0, Accum, PS_wr_3ci)
      Status = Status .AND. PStat_Allocate('Read 3c integrals',
     $   PStat_QStat, 0, Accum, PS_rd_3ci)
      Status = Status .AND. PStat_Allocate('Fitting basis xform',
     $   PStat_QStat, 0, Accum, PS_Fit_xf)
      Status = Status .AND. PStat_Allocate('Schwarz init',
     $   PStat_QStat, 0, Accum, PS_Schwarz)
C
      Status = Status .AND. PStat_Allocate('Form 4ci',
     $   PStat_QStat, 0, Accum, PS_Form_4ci)
      Status = Status .AND. PStat_Allocate('Form A & B',
     $   PStat_QStat, 0, Accum, PS_Form_ab)
      Status = Status .AND. PStat_Allocate('Denominators',
     $   PStat_QStat, 0, Accum, PS_Denom)
      Status = Status .AND. PStat_Allocate('Pair Energies',
     $   PStat_QStat, 0, Accum, PS_PrEn)
C
      If (.NOT. Status)
     $   Call ErrQuit('Rimp2_PStat_Init: Unable to allocate timers', 0,
     &       UNKNOWN_ERR)
C
      Return
      End
      Subroutine RIMP2_PStat_Print
      Implicit NONE
#include "pstat.fh"
C
      Call PStat_Print_All
C
      Return
      End
C
      Subroutine RIMP2_PStat_Term
      Implicit NONE
#include "errquit.fh"
C
#include "pstat.fh"
#include "rimp2_timers.fh"
C
      Logical Status
C
      Status = .TRUE.
      Status = Status .AND. PStat_Free(PS_Overall)
      Status = Status .AND. PStat_Free(PS_xform)
      Status = Status .AND. PStat_Free(PS_Energy)
      Status = Status .AND. PStat_Free(PS_int_gen)
      Status = Status .AND. PStat_Free(PS_step12)
      Status = Status .AND. PStat_Free(PS_step3)
      Status = Status .AND. PStat_Free(PS_wr_3ci)
      Status = Status .AND. PStat_Free(PS_rd_3ci)
      Status = Status .AND. PStat_Free(PS_Fit_xf)
      Status = Status .AND. PStat_Free(PS_Schwarz)
C
      Status = Status .AND. PStat_Free(PS_Form_4ci)
      Status = Status .AND. PStat_Free(PS_Form_ab)
      Status = Status .AND. PStat_Free(PS_Denom)
      Status = Status .AND. PStat_Free(PS_PrEn)
C
      If (.NOT. Status)
     $   Call ErrQuit('RIMP2_PStat_Term: Unable to free timers', 0,
     &       UNKNOWN_ERR)
C
      Return
      End
