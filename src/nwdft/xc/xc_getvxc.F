cc AJL/Begin/FDE
C> \ingroup nwdft
C> @{
C>
C> \file xc_getvxc.F
C> Calculate exchange-correlation energy
C> 
C> \brief Calculate the exchange-correlation energy and Fock matrix
C> contributions
C>
C> This driver routine solves for the XC energy and potential (Vxc) via 
C> numerical quadrature methods. The results are obtained either by
C> direct numerical integration or by means of a LSQ fit of the Vxc to
C> a set of Gaussian functions. This fitted function can be used to
C> evaluate Vxc via a summation of a series of 3-center overlap
C> integrals (3OIs). The algorithms are formulated in terms of matrix
C> products. See subsequent subroutines for further explanation.
C>
      Subroutine xc_getvxc(rtdb, Exc, nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, rho_n, fde_option,
     &                   g_dens_fde, cd_coef)
c
      implicit none
#include "errquit.fh"
c      
      integer nExc      !< [Input] The number of energy terms
                        !< - nExc=1: Exc(1) = exchange + correlation 
                        !< - nExc=2: Exc(1) = exchange,
                        !<           Exc(2) = correlation 
      integer iVxc_opt  !< [Input] If 1 then do density fitting for 
                        !< exchange
      integer g_xcinv   !< [Work] GA for the inversion of the fitting
                        !< matrix
      integer g_dens(2) !< [Input] The density matrices, if ipol=1
                        !< g_dens(1)=\f$D^\alpha+D^\beta\f$, else
                        !< g_dens(1)=\f$D^\alpha\f$ and 
                        !< g_dens(2)=\f$D^\beta\f$.
      integer g_vxc(4)  !< [Output] DFT Fock matrix contributions, if
                        !< ipol=1 g_vxc(1)=\f$F^\alpha+F^\beta\f$, else
                        !< g_vxc(1)=\f$F^\alpha\f$ and
                        !< g_vxc(2)=\f$F^\beta\f$.
      logical IOLGC     !< [Input] .TRUE. do not use disk for exchange
                        !< fitting, .FALSE. store data on disk
      integer rtdb      !< [Input] The RTDB handle
      integer fde_option!< [Input] Are we doing an FDE contribution?
                        !< Should be zero for standard, non-zero for
                        !< other values: 
                        !< 1: ESP
                        !< 2: QM/XC_{Nadd}
                        !< 3: FD/XC_{Nadd}
                        !< 4: QM+FD/XC_{Nadd}
                        !< 5: QM/Ts_{Nadd}
                        !< 6: FD/XC_{Nadd}
                        !< 7: QM_FD/XC_{Nadd}
c
#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "bas.fh"
#include "global.fh"
c#include "tcgmsg.fh"
#include "cdft.fh"
c#include "oep.fh"
c#include "dftpara.fh"
#include "util.fh"
c#include "sym.fh"
#include "stdio.fh"
c#include "case.fh"
c
      double precision rho_n    !< [Output] The number of electrons
                                !< obtained by integrating the density
      logical havehfxc
c     
      logical xc_gotxc
      external xc_gotxc
c
c--> XC Energy
c
      double precision Exc(2) !< [Output] The energy terms
                              !< - nExc=1: Exc(1) = exchange +
                              !<   correlation 
                              !< - nExc=2: Exc(1) = exchange,
                              !<           Exc(2) = correlation 

      double precision cd_coef(nbf_cd,*) !< CD fitting coefficients plus
                                         !< Additional space for XC
                                         !< fitting coefficients in ADFT
                                         
c
c This driver routine solves for the XC energy and potential (Vxc) via 
c numerical quadrature methods. The results are obtained either by direct 
c numerical integration or by means of a LSQ fit of the Vxc to a set of 
c Gaussian functions. This fitted function can be used to evaluate Vxc
c via a summation of a series of 3-center overlap integrals (3OIs). The
c algorithms are formulated in terms of matrix products. See subsequent
c subroutines for further explanation.
c 
c              XC Energy and Potential Index Key, Vxc(pq,i)
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      
c              2        1   |   Vxc^up   Vxc^dw   
c              1        2   |   Vxc      
c              2        2   |   Vxc^up   Vxc^dw   
c
c           nTcols = ipol 
c
      integer me,nTrows,nTcols
      integer lTmat,iTmat,g_truevxc(2)
      logical oprint_intermediate_xc, oprint_time 
c
      double precision time1_xc,time2_xc
C AJL/BEGIN
      integer g_dens_fde(2)
C AJL/END
c
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
c      Exc(1)=0.d0
c      Exc(2)=0.d0
      iTmat=0
c
      me=ga_nodeid()
      havehfxc=abs(xfac(1)).gt.1d-8
c
      if (oprint_intermediate_xc)then
c        write(luout,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
c    &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
c        write(luout,*)'g_dens(1),g_vxc(1),IOLGC:',
c    &               g_dens(1),g_vxc(1),IOLGC
         if (me.eq.0) write(luout,*)'Fock XC matrix entering xc_getvxc:'
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
c        call ga_print(g_dens(1))
c        if(ipol.eq.2)call ga_print(g_dens(2))
         call util_flush(luout)
      endif
c
c     Get the DFT exchange-correlation contribution
      if(oprint_time) 
     &      time1_xc=util_cpusec()   ! start xc build time
      if (xc_gotxc()) then  

         ! XCFIT will be true for ADFT as well
         if(xcfit) then       
            nTrows = nbf_xc
            nTcols = ipol

            ! Allocated memory twice as large to hold both
            ! CD and XC coefs for an eventual ADFT calculation
            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols*2,'Tmat',lTmat,
     &           iTmat))call errquit('xc_getv: cannot allocate Tmat',0,
     &       MA_ERR)

            ! ADFT needs the CD coefs to generate the fitted density
            ! on the grid, and will use the rest of the memory to
            ! store the XC coefs
            if (ADFT) then
              call dcopy(nTrows*nTcols,cd_coef,1,dbl_mb(iTmat),1)
              call dfill(nTrows*nTcols,0.D0,
     $                   dbl_mb(iTmat+nTrows*nTcols),1)
            else
              call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
            endif
         endif
         
         if(havehfxc.or.(.not.cdfit)) then
               if(.not.ga_duplicate(g_vxc(1),g_truevxc(1),'g vxc 1'))
     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
               call ga_zero(g_truevxc(1))
               if(ipol.eq.2) then
                  if(.not.ga_duplicate(g_vxc(2),g_truevxc(2),'gv21'))
     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
                  call ga_zero(g_truevxc(2))
               endif
         else
               g_truevxc(1)=g_vxc(1)
               g_truevxc(2)=g_vxc(2)
         endif
c
         if (fde_option.eq.0) then
           call grid_quadv0(rtdb, g_dens, g_truevxc,
     &                    nexc,rho_n,  Exc, dbl_mb(itmat))
         else
           call grid_quadv0_fde(rtdb, g_dens, g_truevxc,
     &          nexc,rho_n,  Exc, dbl_mb(itmat),fde_option,
     &          g_dens_fde)
         endif

         if(havehfxc.or.(.not.cdfit)) then
             call ga_dadd(1d0,g_vxc(1),1d0,g_truevxc(1),g_vxc(1))
             if (.not. ga_destroy(g_truevxc(1))) call errquit(
     &           ' xc_getv: ga_destroy failed ',0, GA_ERR)
             if(ipol.eq.2) then
                 call ga_dadd(1d0,g_vxc(2),1d0,g_truevxc(2),g_vxc(2))
                 if (.not. ga_destroy(g_truevxc(2))) call errquit(
     &               ' xc_getv: ga_destroy failed ',0, GA_ERR)
             endif
         endif
         if(oprint_time) then 
           time2_xc=util_cpusec()   ! end xc build time
c
c        print fock_xc build time
           if (me.eq.0) then
             write(luout,"(4x,'Fock_xc Build Time:',F13.1,'s')") 
     &                 time2_xc-time1_xc
           endif
         end if
c
c        In case we are performing an xc fit calculation
c        including the ADFT flavor
         if(xcfit) then
c     
c     symmetrize the "T" vector
c    
            if (oskel)then
               if (adft) then
                 call sym_vec_symmetrize(
     .                geom,xc_bas_han,Dbl_MB(iTmat+nbf_xc*ipol))
                 if (ipol.gt.1)then
                    call sym_vec_symmetrize(geom, xc_bas_han, 
     &                      Dbl_MB(iTmat+nbf_xc*ipol+nbf_xc))
                 endif
               else
                 call sym_vec_symmetrize(
     .                geom,xc_bas_han,Dbl_MB(iTmat))
                 if (ipol.gt.1)then
                    call sym_vec_symmetrize(geom, xc_bas_han, 
     &                      Dbl_MB(iTmat+nbf_xc))
                 endif
               endif
            endif

            ! ADFT had the XC coefs at the end, we need to copy them
            ! to the beginning for xc_fitv (or change xc_fitv to operate
            ! on the last nTrows*nTcols)
            if (ADFT) then
              call dcopy(nTrows*nTcols,dbl_mb(itmat+nTrows*nTcols),1,
     &                   dbl_mb(iTmat),1)
            endif


            call xc_fitv(rtdb,Dbl_MB(iTmat), nTrows, nTcols,
     &           g_vxc, g_xcinv, IOLGC)  

            if (ADFT) then
              ! Save XC coeficients for a potential gradient evaluation
              call dcopy(ipol*nbf_cd,dbl_mb(iTmat),1,
     $                   cd_coef(1,ipol+1),1)

              ! In order to build the right Fock matrices in dft_fitvc
              ! we need cd_coef_up + cd_coef_dw + xc_coef_up and
              !         cd_coef_up + cd_coef_dw + xc_coef_dw
              if (ipol.gt.1) then
                call daxpy(nbf_cd,1.0d0,cd_coef(1,2),1,cd_coef(1,1),1)
                call dcopy(nbf_cd,cd_coef(1,1),1,cd_coef(1,2),1)
              endif
              call daxpy(ipol*nbf_cd,1.0d0,dbl_mb(iTmat),1,cd_coef,1)
            endif

            if (.not.ma_pop_stack(lTmat))
     &             call errquit('xc_getv: cannot pop stack',0, MA_ERR)
c     
         ! Standard cd fitting needs the coefficients for the full
         ! density in each spin channel to get the right Fock
         ! matrices
         endif
      endif
      if ((.not.adft).and.CDFIT .and. (ipol.gt.1)) then
         call daxpy(nbf_cd,1.0d0,cd_coef(1,2),1,cd_coef(1,1),1)
         call dcopy(nbf_cd,cd_coef(1,1),1,cd_coef(1,2),1)
      endif
c
      if (oprint_intermediate_xc)then
         if (me.eq.0) write(luout,*)'Fock XC matrix leaving xc_getvxc:'
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
         call util_flush(luout)
      endif

      return
      end
C>
C> @}
cc AJL/End
