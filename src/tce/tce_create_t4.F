      subroutine tce_create_t4(d_t4,l_t4_offset,k_t4_offset,size_t4,
     &                         irrep,title)
c
c $Id$
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
c#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
c#include "tce_amps.fh"
c#include "tce_diis.fh"
      integer d_t1             ! SF handle for t1 amplitudes
      integer l_t1_offset      ! Offset for t1 file
      integer k_t1_offset      ! Offset for t1 file
      integer size_t1          ! File size in doubles
      integer d_t2             ! SF handle for t2 amplitudes
      integer l_t2_offset      ! Offset for t2 file
      integer k_t2_offset      ! Offset for t2 file
      integer size_t2          ! File size in doubles
      integer d_t3             ! SF handle for t3 amplitudes
      integer l_t3_offset      ! Offset for t3 file
      integer k_t3_offset      ! Offset for t3 file
      integer size_t3          ! File size in doubles
      integer d_t4             ! SF handle for t4 amplitudes
      integer l_t4_offset      ! Offset for t4 file
      integer k_t4_offset      ! Offset for t4 file
      integer size_t4          ! File size in doubles
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer irrep
      character*8 title
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      if (nodezero) write(LuOut,*) '============================='
c
      cpu = - util_cpusec()
      call tce_filename(title,filename)
      call tce_t4_offset_new(l_t4_offset,k_t4_offset,size_t4,irrep)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_t4_offset = ',l_t4_offset
#endif
      call createfile(filename,d_t4,size_t4)
      call gatoeaf(d_t4)
      cpu = cpu + util_cpusec()
      if (nodezero.and.util_print('t4',print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) title,size_t4
        write(LuOut,9010) title,filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) title,d_t4
#endif
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
