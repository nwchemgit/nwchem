      subroutine tce_filename(stub,filename)
c
c $Id: tce_filename.F,v 1.2 2003-07-25 17:56:56 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "stdio.fh"
#include "tce.fh"
      character*(*) stub
      character*(*) filename
c
      if (ioalg.eq.0) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.1) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.2) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.3) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.4) then
        call util_file_name(stub,.true.,.true.,filename)
      endif
      return
      end



      subroutine tce_filenameindexed(index,stub,filename)
c
c $Id: tce_filename.F,v 1.2 2003-07-25 17:56:56 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "stdio.fh"
#include "tce.fh"
#include "inp.fh"
      integer index
      integer length
      character*(*) stub
      character*255 temporary
      character*(*) filename
c
      if ((index.lt.1).or.(index.gt.999))
     1  call errquit('tce_filenameindexed: illegal index',index)
      call tce_filename(stub,temporary)
      length = inp_strlen(temporary)
      if (length+4.gt.len(filename))
     1  call errquit('tce_filenameindexed: name too long',length)
      filename = temporary
      write(filename(length+1:length+4),1000) index
 1000 format('.',i3.3)
      return
      end
