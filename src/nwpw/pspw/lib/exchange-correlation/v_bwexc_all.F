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

      integer tmp1(2)
      
*     **** lda's ****
      if (gga.eq.0) then
         if (.not.MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)))
     >      call errquit('v_bwexc_all: out of stack memory',0,MA_ERR)

         call vxc(n2ft3d,ispin,dn,xcp,xce,dbl_mb(tmp1(1)))

        if (.not.MA_pop_stack(tmp1(2)))
     >     call errquit('v_bwexc_all: error popping stack',0,MA_ERR)



*     **** gga's ****
      else if ((gga.ge.10).and.(gga.lt.100)) then
         call v_bwexc(gga,n2ft3d,ispin,dn,1.0d0,1.0d0,xcp,xce)
      
*     **** hybrid gga's ****
      else if (gga.eq.110) then
         call v_bwexc(10,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
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

*     **** lda's ****
      if (gga.eq.0) then
         call vxc(n2ft3d,ispin,dn,xcp,xce,tmp1)

*     **** gga's ****
      else if ((gga.ge.10).and.(gga.lt.100)) then
         call v_bwexc(gga,n2ft3d,ispin,dn,1.0d0,1.0d0,xcp,xce)
      
*     **** hybrid gga's ****
      else if (gga.eq.110) then
         call v_bwexc(10,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.112) then
         call v_bwexc(12,n2ft3d,ispin,dn,0.75d0,1.0d0,xcp,xce)
      else if (gga.eq.114) then
         call v_bwexc(14,n2ft3d,ispin,dn,1.00d0,1.0d0,xcp,xce)
      else if (gga.eq.115) then
         call v_bwexc(15,n2ft3d,ispin,dn,0.80d0,1.0d0,xcp,xce)
      end if
        
      return
      end

