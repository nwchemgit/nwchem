      Logical Function PStat_Allocate( Name, Functions, NAcc,
     $   Accums, Handle )
C$Id$
      Implicit NONE
      Character*(*) Name
      Integer Functions, Handle, NAcc, Accums(*)
C
#include "pstatP.fh"
C
      Integer I
#include "bitops.fh"
C
C     Locate next free timer
C
      do handle = 1, pstat_max_items
         if (.not. pstat_active(handle)) goto 10
      enddo
      pstat_allocate = .false.  ! No handles available
      return
 10   continue
c
      PStat_Name( Handle ) = Name
      pstat_active(handle) = .true.
C     
      Do I = 1, PStat_Nr_Stats
         PStat_Functions( Handle, I) =
     &        1.eq.(IAnd( Functions, 2**(I-1) ))
      EndDo
      pstat_times(handle,1) = 0.0d0
      pstat_times(handle,2) = 0.0d0
      pstat_nr_calls(handle) = 0
      pstat_unmatched_starts(handle) = 0
C     
      PStat_Allocate = .TRUE.
c
      End
