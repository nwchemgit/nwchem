      subroutine tce_create_t1(d_t1,l_t1_offset,k_t1_offset,size_t1,
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
      call tce_t1_offset_new(l_t1_offset,k_t1_offset,size_t1,irrep)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_t1_offset = ',l_t1_offset
#endif
c -------
c      if(nodezero) then
c       write(LuOut,*)'T1-number-of-tasks',int_mb(k_t1_offset)
c       call util_flush(LuOut)
c      end if
c -------------
      call createfile(filename,d_t1,size_t1)
      call tce_guess_t1(d_t1,k_t1_offset)
cc      call tce_guess_t1_zero(d_t1,k_t1_offset)
      call reconcilefile(d_t1,size_t1)
      cpu = cpu + util_cpusec()
      if (nodezero.and.util_print('t1',print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) 't1',size_t1
        write(LuOut,9010) 't1',filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) 't1',d_t1
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
