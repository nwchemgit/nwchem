*
* $Id: v_bwexc.F 22964 2012-10-08 16:26:17Z bylaska $
*
      

*    ************************************
*    *                                  *
*    *          v_bwexc_all             *       
*    *                                  *
*    ************************************
      subroutine v_bwexc_all(gga,n2ft3d,ispin,dn,xcp,xce)
      implicit none
      integer gga
      integer n2ft3d     
      integer  ispin     
      real*8  dn(n2ft3d,2)
      real*8  xcp(n2ft3d,2),xce(n2ft3d)
      
      
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwxc.fh"

      integer tmp1(2)
      
*     **** lda's ****
      if ((.not.nwxc_is_on().and.gga.eq.0).or.
     +    (nwxc_is_on().and.nwxc_is_lda())) then
         if (.not.MA_push_get(mt_dbl,(ispin*n2ft3d),'tmp1',
     >                        tmp1(2),tmp1(1)))
     >      call errquit('v_bwexc_all: out of stack memory',0,MA_ERR)

         call vxc(n2ft3d,ispin,dn,xcp,xce,dbl_mb(tmp1(1)))

        if (.not.MA_pop_stack(tmp1(2)))
     >     call errquit('v_bwexc_all: error popping stack',0,MA_ERR)



*     **** gga's ****
      else if ((.not.nwxc_is_on().and.(gga.ge.10).and.(gga.lt.100)).or.
     +         (nwxc_is_on().and.nwxc_is_gga())) then
         call v_bwexc(gga,n2ft3d,ispin,dn,1.0d0,1.0d0,xcp,xce)
      
*     **** meta-gga's ****
      else if (nwxc_is_mgga()) then
         call v_mexc(gga,n2ft3d,ispin,dn,1.0d0,1.0d0,xcp,xce)

*     **** hybrid gga's ****
      else if (gga.eq.110) then
         call v_bwexc(10,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.111) then
         call v_bwexc(11,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.112) then
         call v_bwexc(12,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.114) then
         call v_bwexc(14,n2ft3d,ispin,dn,1.00d0,1.0d0,xcp,xce)
      else if (gga.eq.115) then
         call v_bwexc(15,n2ft3d,ispin,dn,0.80d0,1.0d0,xcp,xce)
      end if
        
      return
      end



*    ************************************
*    *                                  *
*    *          v_bwexc_all_tmp1        *
*    *                                  *
*    ************************************
      subroutine v_bwexc_all_tmp1(gga,n2ft3d,ispin,dn,xcp,xce,tmp1)
      implicit none
      integer gga
      integer n2ft3d
      integer  ispin
      real*8  dn(n2ft3d,2)
      real*8  xcp(n2ft3d,2),xce(n2ft3d)
      real*8  tmp1(n2ft3d)


#include "mafdecls.fh"
#include "errquit.fh"
#include "nwxc.fh"

*     **** lda's ****
      if ((.not.nwxc_is_on().and.gga.eq.0).or.
     +    (nwxc_is_on().and.nwxc_is_lda())) then

         call vxc(n2ft3d,ispin,dn,xcp,xce,tmp1)

*     **** gga's ****
      else if ((.not.nwxc_is_on().and.(gga.ge.10).and.(gga.lt.100)).or.
     +         (nwxc_is_on().and.nwxc_is_gga())) then

         call v_bwexc(gga,n2ft3d,ispin,dn,1.0d0,1.0d0,xcp,xce)
      
*     **** hybrid gga's ****
      else if (gga.eq.110) then
         call v_bwexc(10,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.111) then
         call v_bwexc(11,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.112) then
         call v_bwexc(12,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.114) then
         call v_bwexc(14,n2ft3d,ispin,dn,1.00d0,1.0d0,xcp,xce)
      else if (gga.eq.115) then
         call v_bwexc(15,n2ft3d,ispin,dn,0.80d0,1.0d0,xcp,xce)
      end if
        
      return
      end

*    ************************************
*    *                                  *
*    *          v_bwexc_print           *
*    *                                  *
*    ************************************
      subroutine v_bwexc_print(luout,gga)
      implicit none
      integer luout,gga

#include "errquit.fh"
#include "nwxc.fh"

      logical  control_has_disp
      external control_has_disp
      character*80 control_options_disp
      external     control_options_disp

      if (nwxc_is_on()) then
         call nwxc_print_nwpw()
         return
      endif
      IF (gga.eq.-1) THEN
         write(luout,1131) 'off'
      ELSE IF (gga.eq.0) THEN
         write(luout,1131) 'LDA (Vosko et al) parameterization'
      ELSE IF (gga.eq.10) THEN
         write(luout,1131)
     >   'PBE96 (White and Bird) parameterization'
      ELSE IF (gga.eq.11) THEN
         write(luout,1131)
     >   'BLYP (White and Bird) parameterization'
      ELSE IF (gga.eq.12) THEN
         write(luout,1131)
     >   'revPBE (White and Bird) parameterization'
      ELSE IF (gga.eq.13) THEN
         write(luout,1131)
     >   'PBEsol (White and Bird) parameterization'

      ELSE IF (gga.eq.110) THEN
         write(luout,1131)
     >   'PBE0 (White and Bird) parameterization'
      ELSE IF (gga.eq.111) THEN
         write(luout,1131)
     >   'BLYP0 (White and Bird) parameterization'
      ELSE IF (gga.eq.112) THEN
         write(luout,1131)
     >   'revPBE0 (White and Bird) parameterization'
      ELSE IF (gga.eq.113) THEN
         write(luout,1131)
     >   'BNL (White and Bird) parameterization'
      ELSE IF (gga.eq.114) THEN
         write(luout,1131)
     >   'HSE (White and Bird) parameterization'
      ELSE IF (gga.eq.115) THEN
         write(luout,1131)
     >   'B3LYP (White and Bird) parameterization'

      ELSE IF (gga.eq.200) THEN
         write(luout,1131) 'Hartree-Fock'
      ELSE
         write(luout,1131) 'unknown parameterization'
         call errquit('bad exchange_correlation',0, INPUT_ERR)
      END IF
      if (control_has_disp()) then
        if (index(control_options_disp(),'-old').ne.0) 
     >    write(luout,1132) 'Grimme2'
        if (index(control_options_disp(),'-zero').ne.0) 
     >    write(luout,1132) 'Grimme3'
        if (index(control_options_disp(),'-bj').ne.0) 
     >    write(luout,1132) 'Grimme4'
      end if
      

      return
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' dispersion correction= ',A)
      end
