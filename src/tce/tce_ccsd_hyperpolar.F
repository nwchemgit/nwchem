      subroutine tce_ccsd_hyperpolar(beta,respaxis,
     &  irrep_d1,d_d1,k_d1_offset,
     &  d_f1,k_f1_offset,d_v2,k_v2_offset,
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
      integer i,j,k,l,ii,jj,kk,ll,iii,jjj,kkk,lll
      integer iter,maxiter
      integer irrep_g,irrep_d1(3)
      integer axisA,axisB,axisC
c      integer d_b0              ! Dipole polarizability component handle
c      integer l_b0_offset       ! Offset for b0 file
c      integer k_b0_offset       ! Offset for b0 file
c      integer size_b0           ! File size in doubles
      integer d_d1(3)           ! MO dipole integral file
      integer k_d1_offset(3)    ! Offset for d1 integral file
      integer d_yr1(12)             ! SF handle for yr1 amplitudes
      integer k_yr1_offset(3)      ! Offset for yr1 file
      integer d_yr2(12)             ! SF handle for yr2 amplitudes
      integer k_yr2_offset(3)      ! Offset for yr2 file
      integer d_tr1(12)             ! SF handle for tr1 amplitudes
      integer k_tr1_offset(3)      ! Offset for tr1 file
      integer d_tr2(12)             ! SF handle for tr2 amplitudes
      integer k_tr2_offset(3)      ! Offset for tr2 file
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
      logical nodezero,respaxis(3)
      double precision cpu,wall
!       double precision omega
      double precision ptol
      double precision beta(3,3,3)
      double precision beta0,beta1,beta2,beta3,beta4!,beta5,beta6
      external ddotfile
      character*255 filename
!       character*20 title
      character*4 irrepname
      character*1 axisname(3)
      data axisname/'X','Y','Z'/
      nodezero=(ga_nodeid().eq.0)
      cpu=-util_cpusec()
      wall=-util_wallsec()
      irrep_g=0
      ptol=1e-10
c
      do axisA = 1, 3
        do axisB = 1, 3
          do axisC = 1, 3
            beta(axisA,axisB,axisC)=0.0d0
          enddo
        enddo
      enddo
c
      do axisA = 1, 3
        irrep_a=irrep_d1(axisA)
        do axisB = 1, 3
          irrep_b=irrep_d1(axisB)
          do axisC = 1, 3
            irrep_c=irrep_d1(axisC)
            if ( ( respaxis(axisA).and.respaxis(axisB) )
     &                            .and.respaxis(axisC) ) then
c
#if defined(DEBUG_PRINT)
                write(LuOut,*) "beta_",axisname(axisA),
     &                                 axisname(axisB),axisname(axisC)
#endif
            if (.false.) then
              call tce_ccsd_hyperpolar_static(beta,axisA,axisB,axisC,
     1             d_d1(axisA),k_d1_offset(axisA),
     2             d_d1(axisB),k_d1_offset(axisB),d_f1,k_f1_offset,
     3             d_v2,k_v2_offset,d_t1,k_t1_offset,d_t2,k_t2_offset,
     4             d_y1,k_y1_offset,d_y2,k_y2_offset,
     5             d_tr1(axisA),k_tr1_offset(axisA),
     6             d_tr2(axisA),k_tr2_offset(axisA),
     7             d_tr1(axisB),k_tr1_offset(axisB),
     8             d_tr2(axisB),k_tr2_offset(axisB),
     9             d_tr1(axisC),k_tr1_offset(axisC),
     1             d_tr2(axisC),k_tr2_offset(axisC),
     2             d_yr1(axisA),k_yr1_offset(axisA),
     3             d_yr2(axisA),k_yr2_offset(axisA))
            else
              call tce_ccsd_hyperpolar_shg(beta,axisA,axisB,axisC,
     1             d_d1,k_d1_offset,
     2             d_f1,k_f1_offset,d_v2,k_v2_offset,
     3             d_t1,k_t1_offset,d_t2,k_t2_offset,
     4             d_y1,k_y1_offset,d_y2,k_y2_offset,
     5             d_tr1,k_tr1_offset,d_tr2,k_tr2_offset,
     6             d_yr1,k_yr1_offset,d_yr2,k_yr2_offset)
            endif
c
            endif ! respaxis(axisA,axisB,axisC)
          enddo ! axisB loop
        enddo ! axisB loop
      enddo ! axisA loop
c
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu, wall
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
c
      subroutine tce_ccsd_hyperpolar_static(beta,axisA,axisB,axisC,
     &  d_a1,k_a1_offset,d_b1,k_b1_offset,d_f1,k_f1_offset,
     &  d_v2,k_v2_offset,d_t1,k_t1_offset,d_t2,k_t2_offset,
     &  d_y1,k_y1_offset,d_y2,k_y2_offset,
     &  d_tr1a,k_tr1a_offset,d_tr2a,k_tr2a_offset,
     &  d_tr1b,k_tr1b_offset,d_tr2b,k_tr2b_offset,
     &  d_tr1c,k_tr1c_offset,d_tr2c,k_tr2c_offset,
     &  d_yr1a,k_yr1a_offset,d_yr2a,k_yr2a_offset)
c
c $Id$
c
c Written by Jeff Hammond, February 2008.
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
      integer axisA,axisB,axisC
      integer d_b0              ! Hyperpolarizability component handle
      integer l_b0_offset       ! Offset for b0 file
      integer k_b0_offset       ! Offset for b0 file
      integer size_b0           ! File size in doubles
      integer d_a1              ! SF handle for dipole integrals
      integer k_a1_offset       ! Offset for dipole integral file
      integer d_b1              ! SF handle for dipole integrals
      integer k_b1_offset       ! Offset for dipole integral file
      integer d_yr1a            ! SF handle for yr1 amplitudes
      integer k_yr1a_offset     ! Offset for yr1 file
      integer d_yr2a            ! SF handle for yr2 amplitudes
      integer k_yr2a_offset     ! Offset for yr2 file
      integer d_tr1a,d_tr1b,d_tr1c
      integer k_tr1a_offset,k_tr1b_offset,k_tr1c_offset
      integer d_tr2a,d_tr2b,d_tr2c
      integer k_tr2a_offset,k_tr2b_offset,k_tr2c_offset
      integer d_f1              ! SF handle for MO 1e integrals
      integer k_f1_offset       ! Offset for 1e integral file
      integer d_v2              ! SF handle for MO 2e integrals
      integer k_v2_offset       ! Offset for 2e integral file
      integer d_t1              ! SF handle for t1 amplitudes
      integer k_t1_offset       ! Offset for t1 file
      integer d_t2              ! SF handle for t2 amplitudes
      integer k_t2_offset       ! Offset for t2 file
      integer d_y1              ! Lambda 1 file
      integer k_y1_offset       ! Offset for Lambda 1 file
      integer d_y2              ! Lambda 2 file
      integer k_y2_offset       ! Offset for Lambda 2 file
      logical nodezero
      double precision cpu
      double precision wall
      double precision omega
      double precision beta(3,3,3)
      double precision beta0,beta1,beta2,beta3,beta4
      external ddotfile
      character*255 filename
      character*20 title
      character*4 irrepname
      nodezero=(ga_nodeid().eq.0)
c
      call tce_e_offset(l_b0_offset,k_b0_offset,size_b0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_b0_offset = ',l_b0_offset
#endif
      call tce_filename('b0',filename)
      call createfile(filename,d_b0,size_b0)
c
      beta1=0.0d0
      beta2=0.0d0
      beta3=0.0d0
      beta4=0.0d0
c
      irrep_g=0
c
      irrep_y=irrep_g
      irrep_tra=irrep_a
      irrep_trb=irrep_b
      irrep_trc=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_1"
#endif
      call beta_1(d_b0,d_t1,d_tr1a,d_tr2a,
     1     d_tr1b,d_tr2b,d_tr1c,d_tr2c,
     2     d_v2,d_y1,d_y2,k_b0_offset,k_t1_offset,
     3     k_tr1a_offset,k_tr2a_offset,
     4     k_tr1b_offset,k_tr2b_offset,
     5     k_tr1c_offset,k_tr2c_offset,
     6     k_v2_offset,k_y1_offset,k_y2_offset)
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta1,1,0)
      call tce_zero(d_b0,size_b0)
c
c
      irrep_y=irrep_g
      irrep_d=irrep_a
      irrep_tra=irrep_b
      irrep_trb=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_2"
#endif
      call beta_2(d_a1,d_b0,d_tr1b,d_tr2b, d_tr1c,d_tr2c,d_y1,d_y2,
     1     k_a1_offset,k_b0_offset,k_tr1b_offset,k_tr2b_offset,
     2     k_tr1c_offset,k_tr2c_offset,k_y1_offset,k_y2_offset)
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta2,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_a
      irrep_tra=irrep_b
      irrep_trb=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_3"
#endif
      call beta_3(d_f1,d_b0,d_t1,d_t2,d_tr1b,d_tr2b,d_tr1c,d_tr2c,d_v2,
     1     d_yr1a,d_yr2a,k_f1_offset,k_b0_offset,k_t1_offset,
     2     k_t2_offset,k_tr1b_offset,k_tr2b_offset,k_tr1c_offset,
     3     k_tr2c_offset,k_v2_offset,k_yr1a_offset,k_yr2a_offset)
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta3,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_a
      irrep_d=irrep_b
      irrep_tr=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4"
#endif
      call beta_4(d_b1,d_b0,d_t1,d_t2,d_tr1c,d_tr2c,d_yr1a,d_yr2a,
     1     k_b1_offset,k_b0_offset,k_t1_offset,k_t2_offset,
     2     k_tr1c_offset,k_tr2c_offset,k_yr1a_offset,k_yr2a_offset)
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta4,1,0)
c
      call deletefile(d_b0)
      if (.not.ma_pop_stack(l_b0_offset))
     1     call errquit("tce_ccsd_hyperpolar: MA problem",
     2                  l_b0_offset,MA_ERR)
c
c     This term contributes to only one component
c
      beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta1
c
c     These terms both contribute to three components
c
      beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta2
      beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta2
      beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta2
c
      beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta3
      beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta3
      beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta3
c
c     This term contributes to six components
c
      beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta4
      beta(axisA,axisC,axisB)=beta(axisA,axisC,axisB)-beta4
      beta(axisB,axisA,axisC)=beta(axisB,axisA,axisC)-beta4
      beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta4
      beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta4
      beta(axisC,axisB,axisA)=beta(axisC,axisB,axisA)-beta4
c
      return
      end


      subroutine tce_ccsd_hyperpolar_shg(beta,axisA,axisB,axisC,
     &  d_d1,k_d1_offset,
     &  d_f1,k_f1_offset,d_v2,k_v2_offset,
     &  d_t1,k_t1_offset,d_t2,k_t2_offset,
     &  d_y1,k_y1_offset,d_y2,k_y2_offset,
     &  d_tr1,k_tr1_offset,d_tr2,k_tr2_offset,
     &  d_yr1,k_yr1_offset,d_yr2,k_yr2_offset)
c
c $Id$
c
c Written by Jeff Hammond, February 2008.
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
      integer axisA,axisB,axisC
      integer d_b0              ! Hyperpolarizability component handle
      integer l_b0_offset       ! Offset for b0 file
      integer k_b0_offset       ! Offset for b0 file
      integer size_b0           ! File size in doubles
      integer d_d1(3)           ! SF handle for dipole integrals
      integer k_d1_offset(3)    ! Offset for dipole integral file
      integer d_yr1(12)         ! SF handle for yr1 amplitudes
      integer k_yr1_offset(3)   ! Offset for yr1 file
      integer d_yr2(12)         ! SF handle for yr2 amplitudes
      integer k_yr2_offset(3)   ! Offset for yr2 file
      integer d_tr1(12)         ! SF handle for tr1 amplitudes
      integer k_tr1_offset(3)   ! Offset for tr1 file
      integer d_tr2(12)         ! SF handle for tr2 amplitudes
      integer k_tr2_offset(3)   ! Offset for tr2 file
      integer d_f1              ! SF handle for MO 1e integrals
      integer k_f1_offset       ! Offset for 1e integral file
      integer d_v2              ! SF handle for MO 2e integrals
      integer k_v2_offset       ! Offset for 2e integral file
      integer d_t1              ! SF handle for t1 amplitudes
      integer k_t1_offset       ! Offset for t1 file
      integer d_t2              ! SF handle for t2 amplitudes
      integer k_t2_offset       ! Offset for t2 file
      integer d_y1              ! Lambda 1 file
      integer k_y1_offset       ! Offset for Lambda 1 file
      integer d_y2              ! Lambda 2 file
      integer k_y2_offset       ! Offset for Lambda 2 file
      logical nodezero
      double precision cpu
      double precision wall
      double precision omega
      double precision beta1,beta2,beta3,beta4,beta5,beta6
      double precision beta7,beta8,beta9,beta10,beta11,beta12
      double precision beta(3,3,3)
      external ddotfile
      character*255 filename
      character*20 title
      character*4 irrepname
      character*1 axisname(3)
      data axisname/'X','Y','Z'/
      nodezero=(ga_nodeid().eq.0)
c
      call tce_e_offset(l_b0_offset,k_b0_offset,size_b0)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_b0_offset = ',l_b0_offset
#endif
      call tce_filename('b0',filename)
      call createfile(filename,d_b0,size_b0)
c
      beta1=0.0d0
      beta2=0.0d0
      beta3=0.0d0
      beta4=0.0d0
c
      irrep_g=0
c
      irrep_y=irrep_a
      irrep_d=irrep_b
      irrep_tr=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_a(-2w)*D_b*T_c(+w)"
#endif
      call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1     d_tr1(axisC+0),d_tr2(axisC+0),d_yr1(axisA+9),d_yr2(axisA+9),
     2     k_d1_offset(axisB),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisC),k_tr2_offset(axisC),
     4     k_yr1_offset(axisA),k_yr2_offset(axisA))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta1,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_a(+2w)*D_b*T_c(-w)"
#endif
      call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1     d_tr1(axisC+3),d_tr2(axisC+3),d_yr1(axisA+6),d_yr2(axisA+6),
     2     k_d1_offset(axisB),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisC),k_tr2_offset(axisC),
     4     k_yr1_offset(axisA),k_yr2_offset(axisA))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta2,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_a
      irrep_d=irrep_c
      irrep_tr=irrep_b
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_a(-2w)*D_c*T_b(+w)"
#endif
      call beta_4(d_d1(axisC),d_b0,d_t1,d_t2,
     1     d_tr1(axisB+0),d_tr2(axisB+0),d_yr1(axisA+9),d_yr2(axisA+9),
     2     k_d1_offset(axisC),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta3,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_a(+2w)*D_c*T_b(-w)"
#endif
      call beta_4(d_d1(axisC),d_b0,d_t1,d_t2,
     1     d_tr1(axisB+3),d_tr2(axisB+3),d_yr1(axisA+6),d_yr2(axisA+6),
     2     k_d1_offset(axisC),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta4,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_b
      irrep_d=irrep_a
      irrep_tr=irrep_c
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_b(-w)*D_a*T_c(+w)"
#endif
      call beta_4(d_d1(axisA),d_b0,d_t1,d_t2,
     1     d_tr1(axisC+0),d_tr2(axisC+0),d_yr1(axisB+3),d_yr2(axisB+3),
     2     k_d1_offset(axisA),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisC),k_tr2_offset(axisC),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta5,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_b(+w)*D_a*T_c(-w)"
#endif
      call beta_4(d_d1(axisA),d_b0,d_t1,d_t2,
     1     d_tr1(axisC+3),d_tr2(axisC+3),d_yr1(axisB+0),d_yr2(axisB+0),
     2     k_d1_offset(axisA),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisC),k_tr2_offset(axisC),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta6,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_b
      irrep_d=irrep_c
      irrep_tr=irrep_a
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_b(-w)*D_c*T_a(+2w)"
#endif
      call beta_4(d_d1(axisC),d_b0,d_t1,d_t2,
     1     d_tr1(axisA+6),d_tr2(axisA+6),d_yr1(axisB+3),d_yr2(axisB+3),
     2     k_d1_offset(axisC),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta7,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_b(+w)*D_c*T_a(-2w)"
#endif
      call beta_4(d_d1(axisC),d_b0,d_t1,d_t2,
     1     d_tr1(axisA+9),d_tr2(axisA+9),d_yr1(axisB+0),d_yr2(axisB+0),
     2     k_d1_offset(axisC),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisB),k_yr2_offset(axisB))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta8,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_c
      irrep_d=irrep_a
      irrep_tr=irrep_b
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_c(-w)*D_a*T_b(+w)"
#endif
      call beta_4(d_d1(axisA),d_b0,d_t1,d_t2,
     1     d_tr1(axisB+0),d_tr2(axisB+0),d_yr1(axisC+3),d_yr2(axisC+3),
     2     k_d1_offset(axisA),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisB),k_tr2_offset(axisB),
     4     k_yr1_offset(axisC),k_yr2_offset(axisC))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta9,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_c(+w)*D_a*T_b(-w)"
#endif
      call beta_4(d_d1(axisA),d_b0,d_t1,d_t2,
     1     d_tr1(axisB+3),d_tr2(axisB+3),d_yr1(axisC+0),d_yr2(axisC+0),
     2     k_d1_offset(axisA),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisB),k_tr2_offset(axisB),
     4     k_yr1_offset(axisC),k_yr2_offset(axisC))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta10,1,0)
      call tce_zero(d_b0,size_b0)
c
      irrep_y=irrep_c
      irrep_d=irrep_b
      irrep_tr=irrep_a
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_c(-w)*D_b*T_a(+w)"
#endif
      call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1     d_tr1(axisA+0),d_tr2(axisA+0),d_yr1(axisC+3),d_yr2(axisC+3),
     2     k_d1_offset(axisb),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisC),k_yr2_offset(axisC))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta11,1,0)
      call tce_zero(d_b0,size_b0)
#if defined(DEBUG_PRINT)
      write(LuOut,*) "beta_4 Y_c(+w)*D_b*T_a(-w)"
#endif
      call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1     d_tr1(axisA+3),d_tr2(axisA+3),d_yr1(axisC+0),d_yr2(axisC+0),
     2     k_d1_offset(axisB),k_b0_offset,k_t1_offset,k_t2_offset,
     3     k_tr1_offset(axisA),k_tr2_offset(axisA),
     4     k_yr1_offset(axisC),k_yr2_offset(axisC))
      call reconcilefile(d_b0,size_b0)
      call get_block(d_b0,beta12,1,0)
      call tce_zero(d_b0,size_b0)
c
      call deletefile(d_b0)
      if (.not.ma_pop_stack(l_b0_offset))
     1     call errquit("tce_ccsd_hyperpolar: MA problem",
     2                  l_b0_offset,MA_ERR)
c
! #if defined(DEBUG_PRINT)
      if (abs( min( min(beta1,beta3),min( min(beta5,beta7),
     &         min(beta9,beta11) ) ) ).gt.(1e-10)) then
        write(LuOut,*) "beta_",axisname(axisA),
     &                         axisname(axisB),axisname(axisC)
        write(LuOut,*) "beta1 = ",beta1
        write(LuOut,*) "beta2 = ",beta2
        write(LuOut,*) "beta3 = ",beta3
        write(LuOut,*) "beta4 = ",beta4
        write(LuOut,*) "beta1 = ",beta5
        write(LuOut,*) "beta2 = ",beta6
        write(LuOut,*) "beta3 = ",beta7
        write(LuOut,*) "beta4 = ",beta8
        write(LuOut,*) "beta1 = ",beta9
        write(LuOut,*) "beta2 = ",beta10
        write(LuOut,*) "beta3 = ",beta11
        write(LuOut,*) "beta4 = ",beta12
      endif
! #endif
! c
! c     This term contributes to only one component
! c
!       beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta1
! c
! c     These terms both contribute to three components
! c
!       beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta2
!       beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta2
!       beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta2
! c
!       beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta3
!       beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta3
!       beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta3
! c
! c     This term contributes to six components
! c
!       beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta4
!       beta(axisA,axisC,axisB)=beta(axisA,axisC,axisB)-beta4
!       beta(axisB,axisA,axisC)=beta(axisB,axisA,axisC)-beta4
!       beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta4
!       beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta4
!       beta(axisC,axisB,axisA)=beta(axisC,axisB,axisA)-beta4
c
      return
      end
