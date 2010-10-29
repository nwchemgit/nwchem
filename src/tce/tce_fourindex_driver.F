      subroutine tce_fourindex_driver(rtdb,d_v2,l_v2_offset,k_v2_offset)
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
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_prop.fh"
#include "tce_ints.fh"
c#include "tce_amps.fh"
c#include "tce_diis.fh"
c
c     CI, CC, & MBPT
c
      integer rtdb             ! Run-time database
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer lax_v2_alpha_offset  ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset  ! Offset for 2e integrals if intorb=.true.
      integer lx_o2o_offset,kx_o2o_offset
      integer lx_o2ot_offset,kx_o2ot_offset
      integer lx_v2spin_offset,kx_v2spin_offset
      integer length
c      integer irrep
c      integer irrep_g
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
c     ================================
c     Parallel integral transformation
c     ================================
c
      if(.not.intorb) THEN  !--------------
c --- debug ----
c       if(nodezero) then
c        write(LuOut,*)'old 2e-transformation'
c        call util_flush(LuOut)
c       end if
c --------------
c
c     2-e integrals first half
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      call tce_ao2e(rtdb,d_ao2e)
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()
      if (nodezero.and.util_print('time',print_default))
     1  write(LuOut,9020) cpu, wall
      if (nodezero) call util_flush(LuOut)
c
c     2-e integrals second half
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      call tce_filename('v2',filename)
      call tce_mo2e_offset(l_v2_offset,k_v2_offset,size_2e)
      call createfile(filename,d_v2,size_2e)
      call tce_mo2e(rtdb,d_ao2e,d_v2,k_v2_offset)
      call reconcilefile(d_v2,size_2e)
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()
      if (nodezero.and.util_print('mo2e',print_default)) then
        write(LuOut,*)
        write(LuOut,9000) '2-e',size_2e
        write(LuOut,9010) '2-e',filename(1:20)
c       write(LuOut,9090) '2-e',d_v2
        write(LuOut,9020) cpu, wall
        call util_flush(LuOut)
      endif
c
      end if ! not intorb ------------------
c
c
c
c
c
c
      if(intorb) THEN !----------------------
c
       if(nodezero) then
        write(6,*)'4-electron integrals stored in orbital form'
       end if
       cpu = - util_cpusec()
       wall = - util_wallsec()
       call tce_filename('v2',filename)
       length = 0 
       call tce_mo2e_offset_intorb(lax_v2_alpha_offset,
     &      kax_v2_alpha_offset,size_2e)
c --- debug ---
       call ga_sync()
c -------------
       l_v2_alpha_offset=lax_v2_alpha_offset 
       k_v2_alpha_offset=kax_v2_alpha_offset 
       call createfile(filename,d_v2,size_2e)
       call ga_zero(d_v2)
      if (nodezero) then
        write(LuOut,*)
        write(LuOut,9000) 'v2   ',size_2e
        call util_flush(LuOut)
      endif
       d_v2orb = d_v2
       cpu = - util_cpusec()
       wall = - util_wallsec()
       if(fast2e.eq.1) then
        if(nodezero) 
     &   write(LuOut,9123)fast2e 
        call tce_mo2e_trans(rtdb,d_v2,kax_v2_alpha_offset,size_2e)
       else if(fast2e.eq.2) then
        if(nodezero) then
         write(LuOut,9123)fast2e
         write(LuOut,9124)imaxsize
        end if
        call atpart_nalength
        call tce_mo2e_zones_4a_disk_2s_new(rtdb,d_v2,
     1                                kax_v2_alpha_offset,
     1                                size_2e)
       else if(fast2e.eq.3) then
        if(nodezero) then
         write(LuOut,9123)fast2e
         write(LuOut,9124)imaxsize
        end if
        call atpart_nalength
        call tce_mo2e_zones_4a_disk_ga(rtdb,d_v2,
     1                                kax_v2_alpha_offset,
     1                                size_2e)
       else if(fast2e.eq.4) then
        if(nodezero) then
         write(LuOut,9123)fast2e
         write(LuOut,9124)imaxsize
        end if
        call atpart_nalength
        call  tce_mo2e_zones_4a_disk_ga_chop(rtdb,d_v2,
     1                                kax_v2_alpha_offset,
     1                                size_2e)
       else
        write(LuOut,9123)fast2e 
        call errquit('tce_property: wrong 2emet',fast2e,CALC_ERR)
       end if      


       cpu = cpu + util_cpusec()
       wall = wall + util_wallsec()
      if (nodezero) then
c        write(LuOut,*)
c        write(LuOut,9000) '2-e size',size_2e
c        call util_flush(LuOut)
c        write(LuOut,9010) '2-e',filename(1:20)
c       write(LuOut,9090) '2-e',d_v2
        write(LuOut,9020) cpu, wall
        call util_flush(LuOut)
      endif
      call tce_mo2e_offset_size(length)
c      if(nodezero) then
c      write(LuOut,*)'length of v2spin_offset:',length
c      call util_flush(LuOut)
c      endif
      call e2_offset_2_offset(lx_o2o_offset,kx_o2o_offset,
     &      lx_o2ot_offset,kx_o2ot_offset,
     &      lx_v2spin_offset,kx_v2spin_offset,length) 
c
       l_o2o_offset    = lx_o2o_offset
       k_o2o_offset    = kx_o2o_offset
       l_o2ot_offset   = lx_o2ot_offset
       k_o2ot_offset   = kx_o2ot_offset
       l_v2spin_offset = lx_v2spin_offset
       k_v2spin_offset = kx_v2spin_offset
c
       k_v2_offset = k_v2spin_offset
c
c
      end if ! intorb ----------------------- 
c
c     =========
c     Terminate
c     =========
c
      return
c
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9090 format(1x,A,' file handle = ',i10)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9480 format(1x,'Cpu & wall time / sec for ',A,2f15.1)
 9110 format(1x,'MBPT(0) energy / hartree             = ',f25.15)
 9030 format(/,1x,'MBPT(2) correlation energy / hartree = ',f25.15)
 9040 format(1x,'MBPT(2) total energy / hartree       = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9123 format(1x,'4-index algorithm nr.',i4,1x,'is used') 
 9124 format(1x,'imaxsize = ',2x,i6)
 9125 format(1x,'imaxsize ichop = ',2x,2i6)
 9230 format(1x,f17.13,f18.13,f11.5,2f8.1)
 9420 format(1x,i4,f18.13,2f8.1)
      end
