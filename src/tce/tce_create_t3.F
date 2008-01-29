      subroutine tce_create_t3(d_t3,l_t3_offset,k_t3_offset,size_t3,
     &                         irrep,title)
c
c $Id: tce_create_t3.F,v 1.1 2008-01-29 18:26:12 jhammond Exp $
c
c Written by Jeff Hammond, January 2008.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
c#include "rtdb.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
#include "tce_amps.fh"
c#include "tce_diis.fh"
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer irrep
      character*8 title
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
      cpu = - util_cpusec()
      call tce_filename(title,filename)
      call tce_t3_offset_new(l_t3_offset,k_t3_offset,size_t3,irrep)
      call createfile(filename,d_t3,size_t3)
      call gatoeaf(d_t3)
      cpu = cpu + util_cpusec()
      if (nodezero.and.util_print('t3',print_default)) then
        write(LuOut,*)
        write(LuOut,9000) title,size_t3
        write(LuOut,9010) title,filename(1:120)
c       write(LuOut,9090) title,d_t3
c       write(LuOut,9020) cpu
        call util_flush(LuOut)
      endif
c
      return
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
c 9030 format(/,1x,'MBPT(2) correlation energy / hartree = ',f25.15)
c 9040 format(1x,'MBPT(2) total energy / hartree       = ',f25.15)
 9090 format(1x,A,' file handle = ',i10)
      end
