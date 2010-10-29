      subroutine tce_fock_create(rtdb,d_f1,l_f1_offset,k_f1_offset)
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
      logical recompf          ! True if recompute Fock
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer length
c      integer irrep
c      integer irrep_g
      double precision ref     ! Ground state energy
c      integer dummy            ! Dummy argument for DIIS
c      integer l_shell,l_subshell,l_sh_size
      character*255 filename
      character*4 irrepname
      logical dft_energy
      external dft_energy
      logical scf
      external scf
      integer i
      integer j
      double precision ddotfile
      external ddotfile
c      integer sym_abelian_axis ! Move this and the following to "sym.fh"
c      external sym_abelian_axis
c
      nodezero=(ga_nodeid().eq.0)
c
c     1-e integrals
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      if (.not.rtdb_get(rtdb,'tce:recompf',mt_log,1,recompf))
     1  call errquit('tce_property: failed reading from rtdb',0,
     2  RTDB_ERR)
      if (recompf) then
        call tce_ao1e_fock2e(rtdb,g_ao1e)
      endif
c
c->d3p975
c      do i = 1, ipol
c        if (.not.ga_destroy(g_movecs(i)))
c     1    call errquit('tce_property: GA problem',0,GA_ERR)
c      enddo
c<-d3p975
c
      call tce_filename('f1',filename)
      call tce_mo1e_offset(l_f1_offset,k_f1_offset,size_1e)
      call createfile(filename,d_f1,size_1e)
      if (recompf) then
        call tce_mo1e(g_ao1e,d_f1,k_f1_offset)
      else
        call tce_mo1e_epsilon(d_f1)
      endif
      call gatoeaf(d_f1)
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()
      if (nodezero.and.util_print('mo1e',print_default)) then
        write(LuOut,*)
        if (recompf) then
          write(LuOut,9120) 'Fock matrix recomputed'
        else
          write(LuOut,9120) 'Epsilons used as Fock matrix'
        endif
        write(LuOut,9000) '1-e',size_1e
        write(LuOut,9010) '1-e',filename(1:20)
c       write(LuOut,9090) '1-e',d_f1
        write(LuOut,9020) cpu, wall
        call util_flush(LuOut)
      endif
      return
c
c     ======
c     Format
c     ======
c
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9090 format(1x,A,' file handle = ',i10)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9480 format(1x,'Cpu & wall time / sec for ',A,2f15.1)
 9080 format(1x,A,' total energy / hartree       = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9210 format(/,1x,'Iteration ',i3,' using ',i4,' trial vectors')
 9230 format(1x,f17.13,f18.13,f11.5,2f8.1)
 9240 format(1x,
     1'--------------------------------------------------------------'
     2,/,1x,'Iterations converged')
 9410 format(
     1  1x,'--------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f18.13,2f8.1)
      end
