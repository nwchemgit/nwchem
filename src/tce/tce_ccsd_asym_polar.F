      double precision function tce_ccsd_asym_polar(
     &  irrep_a1,d_a1,k_a1_offset,irrep_b1,
     &  d_t1,k_t1_offset,d_t2,k_t2_offset,
     &  d_y1,k_y1_offset,d_y2,k_y2_offset,
     &  d_tr1,k_tr1_offset,d_tr2,k_tr2_offset,
     &  d_yr1,k_yr1_offset,d_yr2,k_yr2_offset)
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
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
c#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
c#include "tce_amps.fh"
c#include "tce_diis.fh"
c
c     CI, CC, & MBPT
c
      integer iter,maxiter
      integer irrep_g
      integer d_a1,k_a1_offset,irrep_a1,irrep_b1
      integer d_a0                ! Dipole polarizability component handle
      integer l_a0_offset         ! Offset for a0 file
      integer k_a0_offset         ! Offset for a0 file
      integer size_a0             ! File size in doubles
      integer d_yr1             ! SF handle for yr1 amplitudes
      integer k_yr1_offset      ! Offset for yr1 file
      integer size_yr1          ! File size in doubles
      integer d_yr2             ! SF handle for yr2 amplitudes
      integer k_yr2_offset      ! Offset for yr2 file
      integer size_yr2          ! File size in doubles
      integer d_yrr1            ! SF handle for yrr1 amplitudes
      integer d_yrr2            ! SF handle for yrr2 amplitudes
      integer d_tr1             ! SF handle for tr1 amplitudes
      integer k_tr1_offset      ! Offset for tr1 file
      integer size_tr1          ! File size in doubles
      integer d_tr2             ! SF handle for tr2 amplitudes
      integer k_tr2_offset      ! Offset for tr2 file
      integer size_tr2          ! File size in doubles
      integer d_rr1             ! SF handle for rr1 amplitudes
      integer d_rr2             ! SF handle for rr2 amplitudes
      integer d_o1
      integer k_o1_offset
      integer size_o1
      integer d_f1             ! SF handle for MO 1e integrals
      integer k_f1_offset      ! Offset for 1e integral file
      integer size_1e          ! File size in doubles
      integer d_v2             ! SF handle for MO 2e integrals
      integer k_v2_offset      ! Offset for 2e integral file
      integer size_2e          ! File size in doubles
      integer d_t1             ! SF handle for t1 amplitudes
      integer k_t1_offset      ! Offset for t1 file
      integer size_t1          ! File size in doubles
      integer d_t2             ! SF handle for t2 amplitudes
      integer k_t2_offset      ! Offset for t2 file
      integer size_t2          ! File size in doubles
      integer d_y1             ! Lambda 1 file
      integer d_y2             ! Lambda 2 file
      integer k_y1_offset      ! Offset for Lambda 1 file
      integer size_y1          ! File size in doubles
      integer k_y2_offset      ! Offset for Lambda 2 file
      integer size_y2          ! File size in doubles
      logical nodezero
      double precision cpu
      double precision wall
      double precision omega
      double precision alpha0,alpha1,alpha2!,alpha3,alpha4,alpha5,alpha6
      external ddotfile
      character*255 filename
      character*20 title
      character*4 irrepname
      nodezero=(ga_nodeid().eq.0)
      cpu=-util_cpusec()
      wall=-util_wallsec()
      alpha1=0.0d0
      alpha2=0.0d0
c      alpha3=0.0d0
c      alpha4=0.0d0
c      alpha5=0.0d0
c      alpha6=0.0d0
      irrep_g=0
c
      call tce_e_offset(l_a0_offset,k_a0_offset,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_a0_offset = ',l_a0_offset
#endif
      call tce_filename('a0',filename)
      call createfile(filename,d_a0,size_a0)
c
      irrep_d=irrep_a1
      irrep_tr=irrep_b1
      irrep_y=irrep_g
c
      call alpha_1(d_a1,d_a0,d_t1,d_t2,d_tr1,d_tr2,
     2     d_y1,d_y2,k_a1_offset,k_a0_offset,
     4     k_t1_offset,k_t2_offset,k_tr1_offset,
     5     k_tr2_offset,k_y1_offset,k_y2_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha1,1,0)
      call tce_zero(d_a0,size_a0)
c
      irrep_o=irrep_a1
      irrep_y=irrep_b1
c
      call alpha_2_yr(d_a0,d_a1,d_t1,d_t2,d_yr1,d_yr2,
     2     k_a0_offset,k_a1_offset,k_t1_offset,
     3     k_t2_offset,k_yr1_offset,k_yr2_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha2,1,0)
      call tce_zero(d_a0,size_a0)
c
      alpha0=(-0.5d0)*(alpha1+alpha2)
      call deletefile(d_a0)
c
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
#if defined(DEBUG_PRINT)
      if (nodezero) write(LuOut,9020) cpu, wall
#endif
      if (.not.ma_pop_stack(l_a0_offset))
     1     call errquit("tce_ccsd_asym_polar: MA problem",
     2                  l_a0_offset,MA_ERR)
      tce_ccsd_asym_polar=alpha0
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
      end
