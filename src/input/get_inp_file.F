      subroutine get_input_filename(input_filename)
C     $Id: get_inp_file.F,v 1.5 2000-05-10 22:11:13 edo Exp $
      implicit none
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
      character*(*) input_filename
c     
c     Process zero ends up with the name of the input file
c     from the command line arguments or nwchem.nw if the first
c     argument begins with a hypen (-).
c
c     If the file cannot be accessed append .nw
c     
      character*255 test
#if defined(CRAY)
      integer ilen,ierror
#endif
#if defined(EXT_INT)
      integer*4 arg, nargs
#else
      integer arg, nargs
#endif
      integer l
      integer iargc
      external iargc
c     
      nargs = iargc()
      do arg = 1, nargs
#if defined(CRAY)
         call pxfgetarg(arg,test,ilen,ierror)
#else 
         call getarg(arg, test)
#endif
         if (arg .eq. 1) input_filename = test
         write(LuOut,1) arg, test(1:inp_strlen(test))
 1       format(' argument',i3,' = ',a)
      enddo
      if (nargs .ge. 1) then
         if (input_filename(1:1) .eq. '-') 
     $        input_filename = 'nwchem.nw'
      else
         input_filename = 'nwchem.nw'
      endif
c
      open(LuIn,file=input_filename,status='old',err=10)
      close(LuIn)
      return
c
 10   l = inp_strlen(input_filename)
      write(LuOut,*) ' Unable to open ',input_filename(1:l), 
     $     ' --- appending .nw'
      input_filename(l+1:l+3) = '.nw'
c     
      end
