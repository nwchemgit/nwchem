      Subroutine RIMP2_PStat_Init
C$Id: rimp2_timers.F,v 1.4 1995-03-29 06:05:52 gg502 Exp $
      Implicit NONE
C
#include "pstat.fh"
#include "rimp2_timers.fh"
C
      Logical Status
      Integer NAcc, Accum(1)
      Character*1 Names(1)
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
      Status = Status .AND. PStat_Allocate('2c integral generation',
     $   PStat_QStat, 0, Accum, PS_v_gen)
      Status = Status .AND. PStat_Allocate('V^{-1/2}',
     $   PStat_QStat, 0, Accum, PS_VmHalf)
      Status = Status .AND. PStat_Allocate('Schwarz 4c',
     $   PStat_QStat, 0, Accum, PS_Schwarz_4c)
      Status = Status .AND. PStat_Allocate('Schwarz 2c',
     $   PStat_QStat, 0, Accum, PS_Schwarz_2c)
C
      Status = Status .AND. PStat_Allocate('GA Create',
     $   PStat_QStat, 0, Accum, PS_Create)
      Status = Status .AND. PStat_Allocate('Form 4ci',
     $   PStat_QStat, 0, Accum, PS_Form_4ci)
      Status = Status .AND. PStat_Allocate('Form A & B',
     $   PStat_QStat, 0, Accum, PS_Form_ab)
      Status = Status .AND. PStat_Allocate('Denominators',
     $   PStat_QStat, 0, Accum, PS_Denom)
      Status = Status .AND. PStat_Allocate('Pair Energies',
     $   PStat_QStat, 0, Accum, PS_PrEn)
      Status = Status .AND. PStat_Allocate('Det. Energy Blocking',
     $   PStat_QStat, 0, Accum, PS_E_Blocking)

C
      If (.NOT. Status)
     $   Call ErrQuit('Rimp2_PStat_Init: Unable to allocate timers', 0)
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
      Status = Status .AND. PStat_Free(PS_v_gen)
      Status = Status .AND. PStat_Free(PS_VmHalf)
      Status = Status .AND. PStat_Free(PS_Schwarz_4c)
      Status = Status .AND. PStat_Free(PS_Schwarz_2c)
C
      Status = Status .AND. PStat_Free(PS_Create)
      Status = Status .AND. PStat_Free(PS_Form_4ci)
      Status = Status .AND. PStat_Free(PS_Form_ab)
      Status = Status .AND. PStat_Free(PS_Denom)
      Status = Status .AND. PStat_Free(PS_PrEn)
C
      If (.NOT. Status)
     $   Call ErrQuit('RIMP2_PStat_Term: Unable to free timers', 0)
C
      Return
      End
