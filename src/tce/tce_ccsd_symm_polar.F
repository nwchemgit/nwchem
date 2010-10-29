      double precision function tce_ccsd_symm_polar(omega,
     &  irrep_a1,d_a1,k_a1_offset,
     &  irrep_b1,d_b1,k_b1_offset,
     &  d_f1,k_f1_offset,d_v2,k_v2_offset,
     &  d_t1,k_t1_offset,d_t2,k_t2_offset,
     &  d_y1,k_y1_offset,d_y2,k_y2_offset,
     &  d_tr1ap,d_tr1am,k_tr1a_offset,
     &  d_tr2ap,d_tr2am,k_tr2a_offset,
     &  d_tr1bp,d_tr1bm,k_tr1b_offset,
     &  d_tr2bp,d_tr2bm,k_tr2b_offset)
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
      integer iter,maxiter
      integer d_a1,k_a1_offset,d_b1,k_b1_offset
      integer d_a0                ! Dipole polarizability component handle
      integer l_a0_offset         ! Offset for a0 file
      integer k_a0_offset         ! Offset for a0 file
      integer size_a0             ! File size in doubles
      integer d_tr1ap,d_tr1am,k_tr1a_offset
      integer d_tr2ap,d_tr2am,k_tr2a_offset
      integer d_tr1bp,d_tr1bm,k_tr1b_offset
      integer d_tr2bp,d_tr2bm,k_tr2b_offset
      integer d_f1             ! SF handle for MO 1e integrals
      integer k_f1_offset      ! Offset for 1e integral file
      integer d_v2             ! SF handle for MO 2e integrals
      integer k_v2_offset      ! Offset for 2e integral file
      integer d_t1             ! SF handle for t1 amplitudes
      integer k_t1_offset      ! Offset for t1 file
      integer d_t2             ! SF handle for t2 amplitudes
      integer k_t2_offset      ! Offset for t2 file
      integer d_y1             ! Lambda 1 file
      integer k_y1_offset      ! Offset for Lambda 1 file
      integer d_y2             ! Lambda 2 file
      integer k_y2_offset      ! Offset for Lambda 2 file
      integer irrep_a1,irrep_b1
      logical nodezero
      double precision cpu
      double precision wall
      double precision omega
      double precision alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6
      external ddotfile
      character*255 filename
      character*20 title
      character*4 irrepname
      nodezero=(ga_nodeid().eq.0)
      cpu=-util_cpusec()
      wall=-util_wallsec()
! #if defined(DEBUG_PRINT)
!         if(nodezero) then
!         write(LuOut,*) '================================='
!         write(LuOut,*) 'irrep=',irrep_a1,irrep_b1
!         write(LuOut,*) 'd_d1=',d_a1,d_b1
!         write(LuOut,*) 'k_d1=',k_a1_offset,k_b1_offset
!         write(LuOut,*) 'f1=',d_f1,k_f1_offset
!         write(LuOut,*) 'v2=',d_v2,k_v2_offset
!         write(LuOut,*) 't1=',d_t1,k_t1_offset
!         write(LuOut,*) 't2=',d_t2,k_t2_offset
!         write(LuOut,*) 'y1=',d_y1,k_y1_offset
!         write(LuOut,*) 'y2=',d_y2,k_y2_offset
!         write(LuOut,*) 'tr1a=',d_tr1ap, d_tr1am,k_tr1a_offset
!         write(LuOut,*) 'tr2a=',d_tr2ap,d_tr2am,k_tr2a_offset
!         write(LuOut,*) 'tr1b=',d_tr1bp,d_tr1bm,k_tr1b_offset
!         write(LuOut,*) 'tr2b=',d_tr2bp,d_tr2bm,k_tr2b_offset
!         write(LuOut,*) '================================='
!         endif
! #endif
      alpha1=0.0d0
      alpha2=0.0d0
      alpha3=0.0d0
      alpha4=0.0d0
      alpha5=0.0d0
      alpha6=0.0d0
      call tce_e_offset(l_a0_offset,k_a0_offset,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_a0_offset = ',l_a0_offset
#endif
      call tce_filename('a0',filename)
      call createfile(filename,d_a0,size_a0)
c
c     <0|(1+Y)[A,T(B,+)]|0>
c
      irrep_d=irrep_a1
      irrep_tr=irrep_b1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha1 A'
#endif
      call alpha_1(d_a1,d_a0,d_t1,d_t2,d_tr1bp,d_tr2bp,d_y1,d_y2,
     3     k_a1_offset,k_a0_offset,k_t1_offset,k_t2_offset,
     5     k_tr1b_offset,k_tr2b_offset,k_y1_offset,k_y2_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha1,1,0)
      call tce_zero(d_a0,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha1 = ',alpha1
#endif
c
c     <0|(1+Y)[A,T(B,-)]|0>
c
      irrep_d=irrep_a1
      irrep_tr=irrep_b1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha1 B'
#endif
      call alpha_1(d_a1,d_a0,d_t1,d_t2,d_tr1bm,d_tr2bm,d_y1,d_y2,
     3     k_a1_offset,k_a0_offset,k_t1_offset,k_t2_offset,
     5     k_tr1b_offset,k_tr2b_offset,k_y1_offset,k_y2_offset)
c
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha2,1,0)
      call tce_zero(d_a0,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha2 = ',alpha2
#endif
c
c     <0|(1+Y)[B,T(A,+)]|0>
c
      irrep_d=irrep_b1
      irrep_tr=irrep_a1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha1 C'
#endif
      call alpha_1(d_b1,d_a0,d_t1,d_t2,d_tr1ap,d_tr2ap,d_y1,d_y2,
     3     k_b1_offset,k_a0_offset,k_t1_offset,k_t2_offset,
     5     k_tr1a_offset,k_tr2a_offset,k_y1_offset,k_y2_offset)
c
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha3,1,0)
      call tce_zero(d_a0,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha3 = ',alpha3
#endif
c
c     <0|(1+Y)[B,T(A,-)]|0>
c
      irrep_d=irrep_b1
      irrep_tr=irrep_a1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha1 D'
#endif
      call alpha_1(d_b1,d_a0,d_t1,d_t2,d_tr1am,d_tr2am,
     2     d_y1,d_y2,k_b1_offset,k_a0_offset,k_t1_offset,k_t2_offset,
     5     k_tr1a_offset,k_tr2a_offset,k_y1_offset,k_y2_offset)
c
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha4,1,0)
      call tce_zero(d_a0,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha4 = ',alpha4
#endif
c
c     <0|(1+Y)[[H,T(A,+)],T(B,-)]|0>
c
      irrep_tra=irrep_a1
      irrep_trb=irrep_b1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha2 A'
#endif
      call alpha_2(d_f1,d_a0,d_t1,d_t2,d_tr1ap,d_tr2ap,d_tr1bm,d_tr2bm,
     3     d_v2,d_y1,d_y2,k_f1_offset,k_a0_offset,
     4     k_t1_offset,k_t2_offset,k_tr1a_offset,k_tr2a_offset,
     6     k_tr1b_offset,k_tr2b_offset,k_v2_offset,
     7     k_y1_offset,k_y2_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha5,1,0)
      call tce_zero(d_a0,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha5 = ',alpha5
#endif
c
c     <0|(1+Y)[[H,T(A,-)],T(B,+)]|0>
c
      irrep_tra=irrep_a1
      irrep_trb=irrep_b1
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha2 B'
#endif
      call alpha_2(d_f1,d_a0,d_t1,d_t2,d_tr1am,d_tr2am,d_tr1bp,d_tr2bp,
     3     d_v2,d_y1,d_y2,k_f1_offset,k_a0_offset,
     4     k_t1_offset,k_t2_offset,k_tr1a_offset,k_tr2a_offset,
     6     k_tr1b_offset,k_tr2b_offset,k_v2_offset,
     7     k_y1_offset,k_y2_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha6,1,0)
      call deletefile(d_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'alpha6 = ',alpha6
#endif
c
      alpha0=-0.5d0*(alpha1+alpha2+alpha3+alpha4+alpha5+alpha6)
c
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu, wall
      if (.not.ma_pop_stack(l_a0_offset))
     1     call errquit("tce_ccsd_symm_polar: MA problem",
     2                  l_a0_offset,MA_ERR)
      tce_ccsd_symm_polar=alpha0
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
