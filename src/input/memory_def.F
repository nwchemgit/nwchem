C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     memory_defaults -- Set default memory partition based on platform
C     and other factors
C
C REVISION
C     $Id: memory_def.F,v 1.29 2008-06-20 20:22:35 d3p307 Exp $
C
C SYNOPSIS
      Subroutine Memory_defaults(Total, TSource, Heap, HSource,
     $   Stack, SSource, Global, GSource)
      Implicit NONE
#include "errquit.fh"
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
C     On input, ?Source = 'U' indicates user-provided input.  All others
C     must be 'N' for "not specified".  On output, notations may be:
C        N)ot specified
C        U)ser-specified
C        D)efault (for total only)
C        I)nferred from default partitioning and/or user input
C        C)alculated from other inputs
C        R)eset of user-specified or default value based on other inputs 
C          (for total only)
C
C EXTERNAL ROUTINES
#include "inp.fh"
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
#elif defined(SP1) || defined(LAPI)
C     This is probably a minimum configuration for most IBM
C     platforms, but many will have more.  Assumes dedicated node
C     access.
C
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(LINUX) || defined(WIN32) || defined (MACX)
*
*     parameter ( Total_Dflt = 16777216) ! 128 MB (linux)
*      Parameter ( Total_Dflt = 26214400) ! 200 MB on 256 MB node
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(SUN) || defined(SOLARIS) || defined(HPUX)
*     Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node (sun or solaris)
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(__crayx1)
C
      Parameter ( Total_Dflt = 52428800) ! 512
#elif defined(CRAY_T3E)
C
*too large   Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node
      Parameter ( Total_Dflt = 11010052) ! 83 MB on 128 MB node (t3e)
#elif defined(CRAY_T3D)
C     Note this is a minimum memory configuration for a T3D
C     Code seems to be very large on T3D (~10-15 MB)
C
      Parameter ( Total_Dflt = 5242880 ) ! 40 MB on a 64 MB node (t3d)
#elif defined(KSR)
      Parameter ( Total_Dflt = 2621440 ) ! 20 MB on a 32 MB node  (ksr)
#elif defined(__PARAGON__)
C     The OS Version 1.3 release notes give the following table
C     for memory usage by the OS based on physical memory on the
C     nodes.  The increase is due to larger paging tables.
C     Code itself is about 5 MB with everything configured in.
C     Phys. Mem.   OS Uses  Msg. bufs.   Available memory
C      16 MB        9 MB    1 MB           6 MB      786 432 dbl
C      32 MB       10 MB    1 MB          21 MB    2 763 264 dbl
C      64 MB       11 MB    1 MB          52 MB    6 815 744 dbl
C     128 MB       12 MB    1 MB         115 MB   15 973 289 dbl
C
      Parameter ( Total_Dflt = 2097152 ) ! 16 MB on a 32 MB node (paragon)
#elif defined(DELTA)
C     Note:  RTDB presently needs to be able to malloc to function
C     so this limit leaves it a little room.  It is pretty much
C     an upper limit on memory usage on the Delta, unfortunately.
C
      Parameter (Total_Dflt = 917504 ) ! 7 MB on a 16 MB node (delta)
#elif defined(SGITFP)
*     Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node (sgitfp)
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(DECOSF)
*     Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node (decosf)
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(SGI) && ! defined(SGITFP)
*     Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node (sgi)
      Parameter ( Total_Dflt = 52428800) ! 400 MB on 512 MB 
#elif defined(IBM)  || defined(IBM64)
*     Parameter ( Total_Dflt = 11796480) ! 90 MB on 128 MB node (ibm)
      Parameter ( Total_Dflt = 92428800) ! 400 MB on 512 MB 
#elif defined(FUJITSU_VPP)
      Parameter ( Total_Dflt = 67108864) ! 512 MB (Fujitsu VPP)
#elif defined(BGL)
      Parameter ( Total_Dflt = 52428800) ! 400 on 512 MB (default co mode)
#elif defined(BGP)
      Parameter ( Total_Dflt = 52428800) ! 400 on 512 MB (per core)
#endif
C
C     Fraction of total memory which goes into each memory category
C     by default.  This is (should be) an empirical estimate of what
C     distribution of memory will serve most applications.
C
      Double Precision Heap_Frc, Stack_Frc, Global_Frc
#ifdef ELAN3HACK
      Parameter (Heap_Frc=7d0/18d0,Stack_Frc=7d0/18d0,
     ,     Global_Frc=4d0/18d0)
#else
      Parameter (Heap_Frc = 0.25, Stack_Frc = 0.25, Global_Frc = 0.50 )
#endif
C
C LOCAL VARIABLES
      Integer NSet, Missing, Borrow
      Double precision Renorm
      Logical TGiven, HGiven, SGiven, GGiven
C-----------------------------------------------------------------------
C
C     Check inputs
C
      TGiven = inp_compare(.false., TSource, 'U')
      SGiven = inp_compare(.false., SSource, 'U')
      HGiven = inp_compare(.false., HSource, 'U')
      GGiven = inp_compare(.false., GSource, 'U')
      If ( .NOT. ( TGiven .OR. inp_compare(.false.,TSource, 'N')))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 2,
     &       INPUT_ERR)
      If ( .NOT. ( SGiven .OR. inp_compare(.false.,SSource, 'N')))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 4,
     &       INPUT_ERR)
      If ( .NOT. ( HGiven .OR. inp_compare(.false.,HSource, 'N')))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 6,
     &       INPUT_ERR)
      If ( .NOT. ( GGiven .OR. inp_compare(.false.,GSource, 'N')))
     $   Call ErrQuit('Memory_Defaults: Bad value for argument', 8,
     &       INPUT_ERR)
C
C     Determine how many of (Stack, Heap, Global) the user has
C     specified.
C
      NSet = 0
      If ( SGiven ) NSet = NSet + 1
      If ( HGiven ) NSet = NSet + 1
      If ( GGiven ) NSet = NSet + 1
C
C     Use the default total limit if the user did not provide one.
C     This may be overridden below if NSet = 3.
C
      If ( .NOT. TGiven) then
         Total = Total_Dflt
         TSource = 'D' ! Default
      EndIf
C
C     If S, H, and G are all set by the user, then either T is
C     determined from them, or it is also already set.
C
      If ( NSet .eq. 3 ) then
         If ( .NOT. TGiven .OR. ( TGiven 
     $      .AND. Total .gt. Stack + Heap + Global) ) then
            Total = Stack + Heap + Global
            TSource = 'R'       ! Reset from other inputs
         EndIf
      EndIf
C
C     If two of (Stack, Heap, Global) have been specified, the third
C     can now be found by subtraction.
C
      If ( NSet .eq. 2) then
         Missing = Total - Max( Heap, 0) - Max( Stack, 0)
     $      - Max( Global, 0)
         If ( .NOT. SGiven ) then
            Stack = Missing
            SSource = 'C'
         EndIf
         If ( .NOT. HGiven ) then
            Heap = Missing
            HSource = 'C'
         EndIf
         If ( .NOT. GGiven ) then
            Global = Missing
            GSource = 'C'
         EndIf
      EndIF
C
C     If none or one of (Heap, Stack, Global) are specified, we can
C     use the partitioning formula to figure out the rest.
C
      If ( NSet .le. 1) then
C
C        First pass: Find out how much memory has to be borrowed
C        from the other categories to keep the total constant.
C        Note that Borrow may be negative.  Also, figure out the
C        normalization for the computing the distribution of
C        Borrow among the other categories.
C
         Borrow = 0
         Renorm = 1.0
         If ( HGiven ) then
            Borrow = Heap - (Heap_Frc * Total)
            Renorm = 1.0 - Heap_Frc
         EndIf
         If ( SGiven ) then
            Borrow = Stack - (Stack_Frc * Total)
            Renorm = 1.0 - Stack_Frc
         EndIf
         If ( GGiven ) then
            Borrow = Global - (Global_Frc * Total)
            Renorm = 1.0 - Global_Frc
         EndIf
C
C        Second pass: Set unspecified values based on defaults
C        and what must be borrowed.  Note that if Borrow is zero
C        (if nothing is specified by the user), each category
C        just gets its fraction of the total.
C
         If ( .NOT. HGiven ) then
            Heap   = (Heap_Frc*Total)   - (Borrow*Heap_Frc/Renorm)
            HSource = 'I'
         EndIf
         If ( .NOT. SGiven ) then
            Stack  = (Stack_Frc*Total)  - (Borrow*Stack_Frc/Renorm)
            SSource = 'I'
         EndIf
         If ( .NOT. GGiven ) then
            Global = (Global_Frc*Total) - (Borrow*Global_Frc/Renorm)
            GSource = 'I'
         EndIf
      EndIf
*
*        rjh: ensure that heap+stack are at least 6 mb
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
