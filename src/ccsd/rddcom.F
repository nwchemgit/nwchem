      subroutine rddcom()
      implicit integer (a-z)
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      Logical NonStd
      common/debug/doa(15),dob(15),dog(15),doh(15),dojk(15),
     &             dos(15),dod(15)
C
C     Check for any non-standard setting of the debugging switches
C
      NonStd = .FALSE.
      do i=1,15
         NonStd = NonStd .or. doa(i) .ne. 2
         NonStd = NonStd .or. dob(i) .ne. 2
         NonStd = NonStd .or. dog(i) .ne. 1
         NonStd = NonStd .or. doh(i) .ne. 1
         NonStd = NonStd .or. dojk(i) .ne. 2
         NonStd = NonStd .or. dos(i) .ne. 1
         NonStd = NonStd .or. dod(i) .ne. 1
      enddo
c
C     Print term switches if requested or if anything is non-standard
C
      call ga_sync()
      if ( (NonStd .OR. Util_Print('term debug switches', Print_Debug) )
     $   .AND. ga_nodeid().eq.0)then
         Write (LuOut, '(A)')
         Call Banner(LuOut, 'Term Debugging Switches', '-',
     $      .TRUE., .TRUE., .FALSE.)
         If ( NonStd ) Call Banner(LuOut,
     $      'WARNING: Non-standard calculation requested!', '*',
     $      .TRUE., .TRUE., .TRUE.)
         write(LuOut,57)' doa  ',doa
         write(LuOut,57)' dob  ',dob
         write(LuOut,57)' dog  ',dog
         write(LuOut,57)' doh  ',doh
         write(LuOut,57)' dojk ',dojk
         write(LuOut,57)' dos  ',dos
         write(LuOut,57)' dod  ',dod
   57 format(1x,a6,1x,15i2)
      endif
      call ga_sync()
c
      return
      end

