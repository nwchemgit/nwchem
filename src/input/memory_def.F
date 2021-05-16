C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     memory_defaults -- Set default memory partition based on platform
C     and other factors
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Memory_defaults(Total, TSource, Heap, HSource,
     $   Stack, SSource, Global, GSource)
      Implicit NONE
#include "errquit.fh"
#include "mem_what.fh"
      Integer Total  ! [inout] Total memory limit (in doubles)
      Integer Heap   ! [inout] Heap memory limit (in doubles)
      Integer Stack  ! [inout] Stack memory limit (in doubles)
      Integer Global ! [inout] GA memory limit (in doubles)
      Character*1 TSource ! [input] Source of memory limit
      Character*1 HSource ! [input] Source of memory limit
      Character*1 SSource ! [input] Source of memory limit
      Character*1 GSource ! [input] Source of memory limit
C
C DESCRIPTION
C     Combines user-specified memory limits with (hopefully intelligent)
C     defaults to produce a complete memory specification.
C
C     As memory limits may now be specified by compiled in defaults,
C     configuration files (i.e. .nwchemrc), environment variables,
C     and/or the input file a priority system is needed to resolve
C     inconsistencies. The priorities are (from highest to lowest):
C     1. User input file
C     2. Environment variables
C     3. Configuration files (i.e. .nwchemrc)
C     4. Defaults
C     5. Not specified (i.e. always override)
C     If inconsistencies are discovered in the memory specification
C     the lowest priority quantit(y/ies) will be modified to resolve
C     them.
C
C     On input, Source = 'U' indicates user-provided input.  All others
C     must be 'N' for "not specified".  On output, notations may be:
C        N)ot specified
C        F)ile .nwchemrc
C        E)nvironment
C        U)ser-specified (input file)
C        D)efault (for total only)
C        I)nferred from default partitioning and/or user input
C        C)alculated from other inputs
C        R)eset of user-specified or default value based on other inputs 
C          (for total only)
C
C EXTERNAL ROUTINES
#include "inp.fh"
      logical  memory_env_get, memory_nwchemrc_get, memory_override
      external memory_env_get, memory_nwchemrc_get, memory_override
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
C
      Integer Total_Dflt   ! Default for total memory (doubles)
C
C     Allow the default to be specified as a compiler option
C     at build time.  If it is not specified, we have sensible
C     platform-specific defaults below.
c
c     RJH: adjusted defaults to provide more local memory in 
c          light of new integrals, caching and replicated data
C
#ifdef DFLT_TOT_MEM
      Parameter (Total_Dflt = DFLT_TOT_MEM )   ! from compile flag
#else
      Parameter ( Total_Dflt = 52428800) ! 400 on 512 MB (per core)
#endif
C
C     Fraction of total memory which goes into each memory category
C     by default.  This is (should be) an empirical estimate of what
C     distribution of memory will serve most applications.
C
      Double Precision Heap_Frc, Stack_Frc, Global_Frc
#ifdef MORE_MA_MEM
      Parameter (Heap_Frc=7d0/18d0,Stack_Frc=7d0/18d0,
     ,     Global_Frc=4d0/18d0)
#else
      Parameter (Heap_Frc = 0.25, Stack_Frc = 0.25, Global_Frc = 0.50 )
#endif
C
C LOCAL VARIABLES
      Logical TGiven, HGiven, SGiven, GGiven
      Logical TOverride, HOverride, SOverride, GOverride
C-----------------------------------------------------------------------
C
C     Check inputs
C
      TGiven = inp_compare(.false., TSource, 'U')
      SGiven = inp_compare(.false., SSource, 'U')
      HGiven = inp_compare(.false., HSource, 'U')
      GGiven = inp_compare(.false., GSource, 'U')
      If ( .NOT. TGiven ) Then
        TGiven = memory_env_get(MEM_TOTAL,Total)
        If (TGiven) TSource = 'E'
      EndIf
      If ( .NOT. TGiven ) Then
        TGiven = memory_nwchemrc_get(MEM_TOTAL,Total)
        If (TGiven) TSource = 'F'
      EndIf
      If ( .NOT. HGiven ) Then
        HGiven = memory_env_get(MEM_HEAP,Heap)
        If (HGiven) HSource = 'E'
      EndIf
      If ( .NOT. HGiven ) Then
        HGiven = memory_nwchemrc_get(MEM_HEAP,Heap)
        If (HGiven) HSource = 'F'
      EndIf
      If ( .NOT. SGiven ) Then
        SGiven = memory_env_get(MEM_STACK,Stack)
        If (SGiven) SSource = 'E'
      EndIf
      If ( .NOT. SGiven ) Then
        SGiven = memory_nwchemrc_get(MEM_STACK,Stack)
        If (SGiven) SSource = 'F'
      EndIf
      If ( .NOT. GGiven ) Then
        GGiven = memory_env_get(MEM_GLOBAL,Global)
        If (GGiven) GSource = 'E'
      EndIf
      If ( .NOT. GGiven ) Then
        GGiven = memory_nwchemrc_get(MEM_GLOBAL,Global)
        If (GGiven) GSource = 'F'
      EndIf
      If ( .NOT. 
     &    (( TGiven .AND..NOT. inp_compare(.false.,TSource, 'N')).OR.
     &     ( .NOT.TGiven .AND. inp_compare(.false.,TSource, 'N'))))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 2,
     &       INPUT_ERR)
      If ( .NOT. 
     &    (( SGiven .AND..NOT. inp_compare(.false.,SSource, 'N')).OR.
     &     ( .NOT.SGiven .AND. inp_compare(.false.,SSource, 'N'))))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 4,
     &       INPUT_ERR)
      If ( .NOT. 
     &    (( HGiven .AND..NOT. inp_compare(.false.,HSource, 'N')).OR.
     &     ( .NOT.HGiven .AND. inp_compare(.false.,HSource, 'N'))))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 6,
     &       INPUT_ERR)
      If ( .NOT. 
     &    (( GGiven .AND..NOT. inp_compare(.false.,GSource, 'N')).OR.
     &     ( .NOT.GGiven .AND. inp_compare(.false.,GSource, 'N'))))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 8,
     &       INPUT_ERR)
C
C     Use the default total limit if the user did not provide one.
C
      If ( .NOT. TGiven) then
         Total = Total_Dflt
         TSource = 'D' ! Default
      EndIf
c
c     Set up defaults, the remaining inputs are considered as
c     modifications of the defaults
c
      If (inp_compare(.false.,GSource, 'N')) Then
        Global = Global_Frc*Total_Dflt
      EndIf
      If (inp_compare(.false.,HSource, 'N')) Then
        Heap = Heap_Frc*Total_Dflt
      EndIf
      If (inp_compare(.false.,SSource, 'N')) Then
        Stack = Stack_Frc*Total_Dflt
      EndIf
c
C
C     Work out which variable(s) may be overridden based on the sources
C     of the others.
C
      TOverride = memory_override(TSource,HSource,SSource,GSource)
      SOverride = memory_override(SSource,TSource,HSource,GSource)
      HOverride = memory_override(HSource,TSource,SSource,GSource)
      GOverride = memory_override(GSource,TSource,SSource,HSource)
c
c     There are 16 potential cases which we work out below
c
      if (TOverride) then
        if (SOverride) then
          if (HOverride) then
            if (GOverride) then
c             T T T T: all settings have equal priority so they must
c                      match
              if (Total .ne. Stack + Heap + Global) then
                Call ErrQuit(
     &          'Memory_Defaults: Inconsistent memory specification',
     &          Total-Stack-Heap-Global, INPUT_ERR)
              endif
            else  ! GOverride
c             T T T F: global is set, all others to be computed.
c                      compute total as sum of the rest.
              TSource = 'C'
              Total = Stack+Heap+Global
            endif ! GOverride
          else  ! HOverride
            if (GOverride) then
c             T T F T: heap is set, all others to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            else  ! GOverride
c             T T F F: heap and global are set all others to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            endif ! GOverride
          endif ! HOverride
        else  ! SOverride
          if (HOverride) then
            if (GOverride) then
c             T F T T: stack is set, all others to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            else  ! GOverride
c             T F T F: stack and global are set, all others to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            endif ! GOverride
          else  ! HOverride
            if (GOverride) then
c             T F F T: stack and heap are set, all others to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            else  ! GOverride
c             T F F F: stack, heap and global are set, total to be computed.
c                      compute total as sum of the rest
              TSource = 'C'
              Total = Stack+Heap+Global
            endif ! GOverride
          endif ! HOverride
        endif ! SOVerride
      else  ! TOverride
        if (SOverride) then
          if (HOverride) then
            if (GOverride) then
c             F T T T: total is set, all others are to be computed.
c                      split total according to standard fractions
              SSource = 'C'
              HSource = 'C'
              GSource = 'C'
              Stack   = Stack_Frc*Total
              Heap    = Heap_Frc*Total
              Global  = Global_Frc*Total
            else  ! GOverride
c             F T T F: total and global are set, all others are to be computed.
c                      split remainder of total according to standard fractions
              SSource = 'C'
              HSource = 'C'
              Stack   = Stack_Frc*(Total-Global)/(Stack_Frc+Heap_Frc)
              Heap    = Heap_Frc*(Total-Global)/(Stack_Frc+Heap_Frc)
            endif ! GOverride
          else  ! HOverride
            if (GOverride) then
c             F T F T: total and heap are set, all others are to be computed.
c                      split remainder of total according to standard fractions
              SSource = 'C'
              GSource = 'C'
              Stack   = Stack_Frc*(Total-Heap)/(Stack_Frc+Global_Frc)
              Global  = Global_Frc*(Total-Heap)/(Stack_Frc+Global_Frc)
            else  ! GOverride
c             F T F F: total, heap and global are set, stack is to be computed.
              SSource = 'C'
              Stack   = Total-Heap-Global
            endif ! GOverride
          endif ! HOverride
        else  ! SOverride
          if (HOverride) then
            if (GOverride) then
c             F F T T: total and stack are set, all others are to be computed.
c                      split remainder of total according to standard fractions
              HSource = 'C'
              GSource = 'C'
              Heap    = Heap_Frc*(Total-Stack)/(Global_Frc+Heap_Frc)
              Global  = Global_Frc*(Total-Stack)/(Global_Frc+Heap_Frc)
            else  ! GOverride
c             F F T F: total, stack and global are set, heap is to be computed.
              HSource = 'C'
              Heap    = Total-Stack-Global
            endif ! GOverride
          else  ! HOverride
            if (GOverride) then
c             F F F T: total, stack and heap are set, global is to be computed.
              GSource = 'C'
              Global  = Total-Stack-Heap
            else  ! GOverride
c             F F F F: total, stack, heap and global are set, the settings
c                      better match
              if (Total .ne. Stack + Heap + Global) then
                Call ErrQuit(
     &          'Memory_Defaults: Inconsistent memory specification',
     &          Total-Stack-Heap-Global, INPUT_ERR)
              endif
            endif ! GOverride
          endif ! HOverride
        endif ! SOVerride
      endif ! TOverride
c
      Heap   = max(0,Heap)
      Stack  = max(0,Stack)
      Global = max(0,Global)
      Total  = max(0,Total)
      if (Total .ne. Stack + Heap + Global) then
        Total   = Stack + Heap + Global
        TSource = 'R'
      endif
*
*     rjh: ensure that heap+stack are at least 6 mb
*
      if (heap+stack .lt. 750000) then
         heap = 750000/2
         stack= 750000/2
         global=total-heap-stack
      endif
*
c$$$      Write (0, *) 'memory_defaults: T H S G ',
c$$$     $     Total, TSource, Heap, HSource, Stack, SSource, Global,
c$$$     $   GSource
C
C     Check that what the user entered (or what we guessed) makes
C     sense.

      If ( Total .lt. Stack + Heap + Global ) Call ErrQuit(
     $   'Memory_Defaults: Inconsistent memory specification',
     $   Total-Stack-Heap-Global, INPUT_ERR)
C
      If ( Total .lt. 0) Call ErrQuit(
     $   'Memory_Defaults: Total memory limit unreasonable', Total,
     &       INPUT_ERR)
      If ( Stack .lt. 0) Call ErrQuit(
     $   'Memory_Defaults: Stack memory limit unreasonable', Stack,
     &       INPUT_ERR)
      If ( Heap .lt. 0) Call ErrQuit(
     $   'Memory_Defaults: Heap memory limit unreasonable', Heap,
     &       INPUT_ERR)
      If ( Global .lt. 0) Call ErrQuit(
     $   'Memory_Defaults: Global memory limit unreasonable', Global,
     &       INPUT_ERR)
C
      Return
      End
c
c-----------------------------------------------------------------------
c
      logical function memory_override(Tcandidate,Tvar1,Tvar2,Tvar3)
      implicit none
c
c     Work out if Tcandidate has the lowest priority in which case
c     it may be overridden, otherwise it cannot be overriden.
c     See subroutine memory_defaults for the definitions of the 
c     priority levels.
c
      character*1 Tcandidate ! [input] The candidate variable
      character*1 Tvar1      ! [input] An other variable
      character*1 Tvar2      ! [input] An other variable
      character*1 Tvar3      ! [input] An other variable
c
c     Local variables
c
      character*1 Highest    ! the highest priority among Tvar*
      character*1 Lowest     ! the lowest priority among Tvar*
c
      Highest = 'N'
      Lowest  = 'U'
      if ((Highest.eq.'N').and.(Tvar1.eq.'F'.or.Tvar1.eq.'E'.or.
     &                          Tvar1.eq.'U')) then
        Highest = Tvar1
      else if ((Highest.eq.'D').and.(Tvar1.eq.'F'.or.Tvar1.eq.'E'.or.
     &                               Tvar1.eq.'U')) then
        Highest = Tvar1
      else if ((Highest.eq.'F').and.(Tvar1.eq.'E'.or.Tvar1.eq.'U')) then
        Highest = Tvar1
      else if ((Highest.eq.'E').and.(Tvar1.eq.'U')) then
        Highest = Tvar1
      endif
      if ((Highest.eq.'N').and.(Tvar2.eq.'F'.or.Tvar2.eq.'E'.or.
     &                          Tvar2.eq.'U')) then
        Highest = Tvar2
      else if ((Highest.eq.'D').and.(Tvar2.eq.'F'.or.Tvar2.eq.'E'.or.
     &                               Tvar2.eq.'U')) then
        Highest = Tvar2
      else if ((Highest.eq.'F').and.(Tvar2.eq.'E'.or.Tvar2.eq.'U')) then
        Highest = Tvar2
      else if ((Highest.eq.'E').and.(Tvar2.eq.'U')) then
        Highest = Tvar2
      endif
      if ((Highest.eq.'N').and.(Tvar3.eq.'F'.or.Tvar3.eq.'E'.or.
     &                          Tvar3.eq.'U')) then
        Highest = Tvar3
      else if ((Highest.eq.'D').and.(Tvar3.eq.'F'.or.Tvar3.eq.'E'.or.
     &                               Tvar3.eq.'U')) then
        Highest = Tvar3
      else if ((Highest.eq.'F').and.(Tvar3.eq.'E'.or.Tvar3.eq.'U')) then
        Highest = Tvar3
      else if ((Highest.eq.'E').and.(Tvar3.eq.'U')) then
        Highest = Tvar3
      endif
c
      if ((Lowest.eq.'U').and.(Tvar1.eq.'N'.or.Tvar1.eq.'D'.or.
     &                         Tvar1.eq.'F'.or.Tvar1.eq.'E')) then
        Lowest = Tvar1
      else if ((Lowest.eq.'E').and.(Tvar1.eq.'N'.or.Tvar1.eq.'D'.or.
     &                              Tvar1.eq.'F')) then
        Lowest = Tvar1
      else if ((Lowest.eq.'F').and.(Tvar1.eq.'N'.or.Tvar1.eq.'D')) then
        Lowest = Tvar1
      else if ((Lowest.eq.'D').and.(Tvar1.eq.'N')) then
        Lowest = Tvar1
      endif
      if ((Lowest.eq.'U').and.(Tvar2.eq.'N'.or.Tvar2.eq.'D'.or.
     &                         Tvar2.eq.'F'.or.Tvar2.eq.'E')) then
        Lowest = Tvar2
      else if ((Lowest.eq.'E').and.(Tvar2.eq.'N'.or.Tvar2.eq.'D'.or.
     &                              Tvar2.eq.'F')) then
        Lowest = Tvar2
      else if ((Lowest.eq.'F').and.(Tvar2.eq.'N'.or.Tvar2.eq.'D')) then
        Lowest = Tvar2
      else if ((Lowest.eq.'D').and.(Tvar2.eq.'N')) then
        Lowest = Tvar2
      endif
      if ((Lowest.eq.'U').and.(Tvar3.eq.'N'.or.Tvar3.eq.'D'.or.
     &                         Tvar3.eq.'F'.or.Tvar3.eq.'E')) then
        Lowest = Tvar3
      else if ((Lowest.eq.'E').and.(Tvar3.eq.'N'.or.Tvar3.eq.'D'.or.
     &                              Tvar3.eq.'F')) then
        Lowest = Tvar3
      else if ((Lowest.eq.'F').and.(Tvar3.eq.'N'.or.Tvar3.eq.'D')) then
        Lowest = Tvar3
      else if ((Lowest.eq.'D').and.(Tvar3.eq.'N')) then
        Lowest = Tvar3
      endif
c
      memory_override = .false.
      if (Tcandidate.eq.'N') then
        memory_override = .true.
      else if (Lowest.eq.Highest) then
        if (Tcandidate.eq.Lowest) then
          memory_override = .false.
        else if ((Lowest.eq.'U').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D'.or.Tcandidate.eq.'F'.or.
     &           Tcandidate.eq.'E')) then
          memory_override = .true.
        else if ((Lowest.eq.'E').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D'.or.Tcandidate.eq.'F')) then
          memory_override = .true.
        else if ((Lowest.eq.'F').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D')) then
          memory_override = .true.
        else if ((Lowest.eq.'D').and.(Tcandidate.eq.'N')) then
          memory_override = .true.
        endif
      else
        if (Tcandidate.eq.Lowest) then
          memory_override = .true.
        else if ((Lowest.eq.'U').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D'.or.Tcandidate.eq.'F'.or.
     &           Tcandidate.eq.'E')) then
          memory_override = .true.
        else if ((Lowest.eq.'E').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D'.or.Tcandidate.eq.'F')) then
          memory_override = .true.
        else if ((Lowest.eq.'F').and.(Tcandidate.eq.'N'.or.
     &           Tcandidate.eq.'D')) then
          memory_override = .true.
        else if ((Lowest.eq.'D').and.(Tcandidate.eq.'N')) then
          memory_override = .true.
        else if ((Lowest.eq.'N').and.(Tcandidate.eq.'D')) then
          memory_override = .true.
        endif
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function memory_env_get(iwhat,iwords)
      implicit none
c
c     Get the number of words of the memory quantity specified by iwhat
c     from environment variables. If an environment variable was found
c     with a valid value the function returns .true. and the
c     corresponding value in iwords.
c     Otherwise the function returns .false. and iwords is unchanged.
c
#include "mem_what.fh"
#include "util.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
c
      integer iwhat  ! [input] id specifying what memory quantity to get
      integer iwords ! [input/output] how many words of memory
c
      logical success
      integer ival
      character*32 strval
      character*32 strenv
      logical debug
c
      debug = .false.
c
      success = .false.
      memory_env_get = success
      if (iwhat.eq.MEM_TOTAL) then
        strenv = "NWCHEM_MEMORY_TOTAL"
      else if (iwhat.eq.MEM_HEAP) then
        strenv = "NWCHEM_MEMORY_HEAP"
      else if (iwhat.eq.MEM_STACK) then
        strenv = "NWCHEM_MEMORY_STACK"
      else if (iwhat.eq.MEM_GLOBAL) then
        strenv = "NWCHEM_MEMORY_GLOBAL"
      else
        if (debug) then
          write(LuOut,
     +    '(i7," :memory_env_get failed: invalid iwhat: ",i8)')
     +    ga_nodeid(),iwhat
        endif
        return
      endif
c
c     Now try and get the value.
c
      call util_getenv(strenv(1:inp_strlen(strenv)),strval)
      if (debug) then
        write(LuOut,*)ga_nodeid(),"memory_env_get: "
     +        //strenv(1:inp_strlen(strenv))//"=<"
     +        //strval(1:inp_strlen(strval))//">"
      endif
      if (inp_strlen(strval).gt.0) then
        read(strval,'(i32)',ERR=900)ival
        if (ival.gt.0) then
          success = .true.
        endif
 900    continue
      endif
c
c     Check the outcome and return results accordingly
c
      memory_env_get = success
      if (success) iwords = ival
      if (debug) then
        write(LuOut,*)ga_nodeid(),"memory_env_get: ",success,
     +                ival,iwords
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function memory_nwchemrc_get(iwhat,iwords)
      implicit none
c
c     Get the number of words of the memory quantity specified by iwhat
c     from nwchemrc files. If an nwchemrc key was found with a valid
c     value the function returns .true. and the corresponding value in
c     iwords.
c     Otherwise the function returns .false. and iwords is unchanged.
c
#include "mem_what.fh"
#include "util.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
c
      integer iwhat  ! [input] id specifying what memory quantity to get
      integer iwords ! [input/output] how many words of memory
c
      logical success
      integer ival
      character*32 strval
      character*32 strkey
      logical debug
c
      debug = .false.
c
      success = .false.
      memory_nwchemrc_get = success
      if (iwhat.eq.MEM_TOTAL) then
        strkey = "memory_total"
      else if (iwhat.eq.MEM_HEAP) then
        strkey = "memory_heap"
      else if (iwhat.eq.MEM_STACK) then
        strkey = "memory_stack"
      else if (iwhat.eq.MEM_GLOBAL) then
        strkey = "memory_global"
      else
        if (debug) then
          write(LuOut,
     +    '(i7," :memory_nwchemrc_get failed: invalid iwhat: ",i8)')
     +    ga_nodeid(),iwhat
        endif
        return
      endif
c
c     Now try and get the value.
c
      if (util_nwchemrc_get(strkey(1:inp_strlen(strkey)),strval)) then
        read(strval,'(i32)',ERR=910)ival
        if (ival.gt.0) then
          success = .true.
        endif
 910    continue
      endif
      if (debug) then
        write(LuOut,*)ga_nodeid(),"memory_nwchemrc_get: "
     +        //strkey(1:inp_strlen(strkey))//"=<"
     +        //strval(1:inp_strlen(strval))//">"
      endif
c
c     Check the outcome and return results accordingly
c
      memory_nwchemrc_get = success
      if (success) iwords = ival
      if (debug) then
        write(LuOut,*)ga_nodeid(),"memory_nwchemrc_get: ",success,
     +                ival,iwords
      endif
c
      return
      end

