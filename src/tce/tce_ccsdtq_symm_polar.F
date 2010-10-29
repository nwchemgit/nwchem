      double precision function tce_ccsdtq_symm_polar(omega,
     &  irrep_a1,d_a1,k_a1_offset,irrep_b1,d_b1,k_b1_offset,
     &  d_f1,k_f1_offset,d_v2,k_v2_offset,
     &  d_t1,k_t1_offset,d_t2,k_t2_offset,
     &  d_t3,k_t3_offset,d_t4,k_t4_offset,
     &  d_y1,k_y1_offset,d_y2,k_y2_offset,
     &  d_y3,k_y3_offset,d_y4,k_y4_offset,
     &  d_tr1ap,d_tr1am,k_tr1a_offset,d_tr2ap,d_tr2am,k_tr2a_offset,
     &  d_tr3ap,d_tr3am,k_tr3a_offset,d_tr4ap,d_tr4am,k_tr4a_offset,
     &  d_tr1bp,d_tr1bm,k_tr1b_offset,d_tr2bp,d_tr2bm,k_tr2b_offset,
     &  d_tr3bp,d_tr3bm,k_tr3b_offset,d_tr4bp,d_tr4bm,k_tr4b_offset)
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
      integer d_tr3ap,d_tr3am,k_tr3a_offset
      integer d_tr4ap,d_tr4am,k_tr4a_offset
      integer d_tr1bp,d_tr1bm,k_tr1b_offset
      integer d_tr2bp,d_tr2bm,k_tr2b_offset
      integer d_tr3bp,d_tr3bm,k_tr3b_offset
      integer d_tr4bp,d_tr4bm,k_tr4b_offset
      integer d_rr1,d_rr2,d_rr3,d_rr4
      integer d_o1,k_o1_offset
      integer d_f1,k_f1_offset
      integer d_v2,k_v2_offset
      integer d_t1,k_t1_offset,size_t1
      integer d_t2,k_t2_offset,size_t2
      integer d_t3,k_t3_offset,size_t3
      integer d_t4,k_t4_offset,size_t4
      integer d_y1,k_y1_offset
      integer d_y2,k_y2_offset
      integer d_y3,k_y3_offset
      integer d_y4,k_y4_offset
      integer irrep_a1,irrep_b1
      logical nodezero,debugprint
      double precision cpu
      double precision wall
      double precision omega
      double precision alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6
      external ddotfile
      character*255 filename
      character*20 title
      character*4 irrepname
      alpha0=0.0d0
#if defined(CCSDTQ)
      alpha1=0.0d0
      alpha2=0.0d0
      alpha3=0.0d0
      alpha4=0.0d0
      alpha5=0.0d0
      alpha6=0.0d0
      nodezero=(ga_nodeid().eq.0)
      cpu=-util_cpusec()
      wall=-util_wallsec()
      call tce_e_offset(l_a0_offset,k_a0_offset,size_a0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_a0_offset = ',l_a0_offset
#endif
      call tce_filename('a0',filename)
      call createfile(filename,d_a0,size_a0)
      irrep_c=irrep_a
      irrep_oa=irrep_a
      irrep_ob=irrep_b
      irrep_tra=irrep_a
      irrep_trb=irrep_b
c      if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
      call ccsdtq_lr_alpha(d_f1,d_a0,d_a1,d_b1,
     1     d_t1,d_t2,d_t3,d_t4,d_tr1am,d_tr2am,d_tr3am,d_tr4am,
     4     d_tr1bp,d_tr2bp,d_tr3bp,d_tr4bp,d_v2,d_y1,d_y2,d_y3,d_y4,
     6     k_f1_offset,k_a0_offset,k_a1_offset,
     7     k_b1_offset,k_t1_offset,k_t2_offset,
     8     k_t3_offset,k_t4_offset,k_tr1a_offset,
     9     k_tr2a_offset,k_tr3a_offset,k_tr4a_offset,k_tr1b_offset,
     2     k_tr2b_offset,k_tr3b_offset,k_tr4b_offset,k_v2_offset,
     3     k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
      call reconcilefile(d_a0,size_a0)
      call get_block(d_a0,alpha0,1,0)
      call tce_zero(d_a0,size_a0)
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu, wall
      if (.not.ma_pop_stack(l_a0_offset))
     1     call errquit("tce_ccsdtq_symm_polar: MA problem",
     2                  l_a0_offset,MA_ERR)
#endif
      tce_ccsdtq_symm_polar=alpha0
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
