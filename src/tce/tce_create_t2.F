      subroutine tce_create_t2(d_t2,l_t2_offset,k_t2_offset,size_t2,
     &                         irrep,title,
     &                         d_t1,k_t1_offset,size_t1,
     &                         d_f1,k_f1_offset,size_1e,
     &                         d_v2,k_v2_offset,size_2e,
     &                         ref,needt1)
c
c $Id: tce_create_t2.F,v 1.1 2008-01-29 18:26:12 jhammond Exp $
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
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_prop.fh"
#include "tce_ints.fh"
#include "tce_amps.fh"
c#include "tce_diis.fh"
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      double precision ref     ! Ground state energy
      double precision corr    ! Correlation energy
      integer irrep
      character*8 title
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
      cpu = - util_cpusec()
      call tce_filename(title,filename)
      call tce_t2_offset_new(l_t2_offset,k_t2_offset,size_t2,irrep)
c ------------
      if(nodezero) then
       write(LuOut,*)'T2-number-of-boxes',int_mb(k_t2_offset)
       call util_flush(LuOut)
      end if
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
      if(needt1.and.(model.eq.'ccsd')) then
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
        write(LuOut,*)
        write(LuOut,9000) title,size_t2
        write(LuOut,9010) title,filename(1:120)
c       write(LuOut,9090) title,d_t2
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
