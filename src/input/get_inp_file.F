      subroutine get_input_filename(input_filename)
C     $Id: get_inp_file.F,v 1.9 2000-08-02 01:26:00 bjohnson Exp $
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
      integer arg,nargs
#if defined(CRAY)
      integer ilen,ierror
#endif
#if defined(EXT_INT) && !defined(HPUX64)
      integer*4 iargc
#else
      integer iargc
#endif
      integer l
      external iargc
c     
      nargs = iargc()
      do arg = 1, nargs
         call util_getarg(arg,test)
         if (arg .eq. 1) input_filename = test
         write(LuOut,1) arg, test(1:inp_strlen(test))
c
c     hack to get nwchem.nw to be read with mpi
c
         if(inp_strlen(test).eq.0.and.arg.eq.1) 
     .        input_filename(1:1)='-'
 1       format(' argument',i3,' = ',a)
      enddo
 2    continue
      if (nargs .ge. 1) then
         if (input_filename(1:1) .eq. '-') 
     $        input_filename = 'nwchem.nw'
      else
         input_filename = 'nwchem.nw'
      endif
c
#if defined(CYGNUS) || defined(WIN32)
      call fix_windows_path(input_filename)
#endif
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
