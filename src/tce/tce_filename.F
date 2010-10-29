      subroutine tce_filename(stub,filename)
c
c $Id$
c
      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
#include "util.fh"
c#include "sf.fh"
c#include "stdio.fh"
c#include "tce.fh"
#include "tce_main.fh"
      character*(*) stub
      character*(*) filename
c
      if (stub.eq.'hbar') then
        call util_file_name(stub,.true.,.true.,filename)
      else if (ioalg.eq.0) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.1) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.2) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.3) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.4) then
        call util_file_name(stub,.true.,.true.,filename)
      else if (ioalg.eq.5) then
        call util_file_name(stub,.true.,.false.,filename)
      else if (ioalg.eq.6) then
        call util_file_name(stub,.true.,.true.,filename)
      endif
      return
      end



      subroutine tce_filenameindexed(index,stub,filename)
c
c $Id$
c
      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
c#include "util.fh"
c#include "sf.fh"
c#include "stdio.fh"
c#include "tce.fh"
c#include "tce_main.fh"
#include "errquit.fh"
#include "inp.fh"
      integer index
      integer length
      character*(*) stub
      character*255 temporary
      character*(*) filename
c
      if ((index.lt.1).or.(index.gt.999))
     1  call errquit('tce_filenameindexed: illegal index',index,
     2  UNKNOWN_ERR)
      call tce_filename(stub,temporary)
      length = inp_strlen(temporary)
      if (length+4.gt.len(filename))
     1  call errquit('tce_filenameindexed: name too long',length,
     2  UNKNOWN_ERR)
      filename = temporary
      write(filename(length+1:length+4),1000) index
 1000 format('.',i3.3)
      return
      end



      subroutine tce_filename_4ind(index,stub,filename)
c
c $Id$
c
      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
#include "util.fh"
c#include "sf.fh"
c#include "stdio.fh"
c#include "tce.fh"
c#include "tce_main.fh"
#include "errquit.fh"
#include "inp.fh"
      integer index
      integer length
      character*(*) stub
      character*255 temporary
      character*(*) filename
c
      if ((index.lt.1).or.(index.gt.999))
     1  call errquit('tce_filenameindexed: illegal index',index,
     2  UNKNOWN_ERR)
ccx      call tce_filename(stub,temporary)
      call util_file_name(stub,.false.,.false.,temporary)
      length = inp_strlen(temporary)
      if (length+4.gt.len(filename))
     1  call errquit('tce_filenameindexed: name too long',length,
     2  UNKNOWN_ERR)
      filename = temporary
      write(filename(length+1:length+4),1000) index
 1000 format('.',i3.3)
      return
      end





      subroutine tce_filename_sf(stub,filename)
c
c $Id$
c
      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
#include "util.fh"
c#include "sf.fh"
c#include "stdio.fh"
c#include "tce.fh"
c#include "tce_main.fh"
      character*(*) stub
      character*(*) filename
c
c      else if (ioalg.eq.3) then
        call util_file_name(stub,.false.,.false.,filename)
c      endif
      return
      end

