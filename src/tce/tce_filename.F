      subroutine tce_filename(stub,filename)
c
c $Id: tce_filename.F,v 1.1 2002-11-11 18:25:55 sohirata Exp $
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
