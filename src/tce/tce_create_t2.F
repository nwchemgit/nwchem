      subroutine tce_create_t2(d_t2,l_t2_offset,k_t2_offset,size_t2,
     &                         irrep,title,
     &                         d_t1,k_t1_offset,size_t1,
     &                         d_f1,k_f1_offset,size_1e,
     &                         d_v2,k_v2_offset,size_2e,
     &                         ref,needt1,ioalg,model)
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
      logical needt1
      integer ioalg
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
      integer d_f1             ! SF handle for MO 1e integrals
      integer l_f1_offset      ! Offset for 1e integral file
      integer k_f1_offset      ! Offset for 1e integral file
      integer size_1e          ! File size in doubles
      integer d_v2             ! SF handle for MO 2e integrals
      integer l_v2_offset      ! Offset for 2e integral file
      integer k_v2_offset      ! Offset for 2e integral file
      integer size_2e          ! File size in doubles
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      double precision ref     ! Ground state energy
      double precision corr    ! Correlation energy
      integer irrep
      character*8 title
      character*10 model
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      if (nodezero) write(LuOut,*) '============================='
c
      cpu = - util_cpusec()
      call tce_filename(title,filename)
      call tce_t2_offset_new(l_t2_offset,k_t2_offset,size_t2,irrep)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_t2_offset = ',l_t2_offset
#endif
c ------------
c      if(nodezero) then
c       write(LuOut,*)'T2-number-of-boxes',int_mb(k_t2_offset)
c       call util_flush(LuOut)
c      end if
c -------------
      call createfile(filename,d_t2,size_t2)
c     if(nodezero) then
c       write(LuOut,*)'before tce_guess_t2'
c       call util_flush(LuOut)
c     end if
      call tce_guess_t2(d_v2,k_v2_offset,d_t2,k_t2_offset)
      call reconcilefile(d_t2,size_t2)
c     if(nodezero) then
c       write(LuOut,*)'after tce_guess_t2'
c       call util_flush(LuOut)
c     end if
c starts for t1 amplitudes t1(in the second order of MBPT)
c (only if ga option is used)
ccccc      if(needt1.and.(model.eq.'ccsd')) then
      if(needt1) then ! fix problem with CCSDT?
        if(ioalg.eq.2) then
          call t1mp2(d_f1,d_t1,d_t2,d_v2,k_f1_offset,k_t1_offset,
     &    k_t2_offset,k_v2_offset)
          call reconcilefile(d_t1,size_t1)
          call tce_guess_t1(d_t1,k_t1_offset)
          call reconcilefile(d_t1,size_t1)
        end if
      end if
c -----
      if (util_print('mbpt2',print_debug)) then
        call tce_mbpt2(d_v2,k_v2_offset,d_t2,k_t2_offset,corr)
        if (nodezero) then
          write(LuOut,9030) corr
          write(LuOut,9040) ref + corr
        endif
      endif
      cpu = cpu + util_cpusec()
      if (nodezero.and.util_print(title,print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) title,size_t2
        write(LuOut,9010) title,filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) title,d_t2
#endif
c       write(LuOut,9020) cpu
        call util_flush(LuOut)
      endif
c
      return
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9030 format(/,1x,'MBPT(2) correlation energy / hartree = ',f25.15)
 9040 format(1x,'MBPT(2) total energy / hartree       = ',f25.15)
 9090 format(1x,A,' file handle = ',i10)
      end
