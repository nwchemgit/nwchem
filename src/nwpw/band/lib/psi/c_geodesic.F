* $Id: c_geodesic.F,v 1.15 2006-12-13 21:53:19 d3p708 Exp $
*


*     ***********************************
*     *			                *
*     *		c_geodesic_init	        *
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value
      integer npack1,nemax,nbrill
      

*     **** external functions **** 
      logical  Pneb_w_allocate
      integer  cpsi_neq,cpsi_nbrillq
      external Pneb_w_allocate
      external cpsi_neq,cpsi_nbrillq


      call Cram_max_npack(npack1)
      nemax  = cpsi_neq(1)+cpsi_neq(2)
      nbrill = cpsi_nbrillq()
      
      value = MA_alloc_get(mt_dcpl,npack1*nemax*nbrill,
     >                     'U',U(2),U(1))
      value = value.and.Pneb_w_allocate(0,0,Vt)
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax*nbrill,
     >                     'S',S(2),S(1))
      if (.not. value) 
     > call errquit('c_geodesic_init:out of heap memory',0, MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *		c_geodesic_finalize	*
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*
      subroutine c_geodesic_finalize()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

      value =           MA_free_heap(S(2))
      value = value.and.Pneb_w_free(Vt)
      value = value.and.MA_free_heap(U(2))
      if (.not. value) 
     > call errquit(
     >  'c_geodesic_finalize:error freeing of heap memory',0, MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		c_geodesic_start        *
*     *					*
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_start(A,max_sigma,dE)
      implicit none
      complex*16 A(*)
      real*8     max_sigma,dE

#include "errquit.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      integer V(2),i
      integer asize,nbrill,nex(2),neall,npack1,indx
      real*8 tmp
      logical value
*     **** external functions **** 
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      real*8   c_electron_eorbit
      integer Pneb_nbrillq,Pneb_ispinq
      external Pneb_w_push_get,Pneb_w_pop_stack
      external c_electron_eorbit
      external Pneb_nbrillq,Pneb_ispinq


      call nwpw_timing_start(10)   
      call Cram_max_npack(npack1)

*     **** allocate tmp space ****
      if (.not.Pneb_w_push_get(0,0,V)) 
     > call errquit('c_geodesic_start: out of stack memory',0,MA_ERR)

*     **** HomeGrown SVD ****
      nbrill=Pneb_nbrillq()
      call Pneb_neq(nex)
      neall=nex(1)+(nex(2)*(Pneb_ispinq()-1))

      call dcopy(nbrill*neall,0.0d0,0,dbl_mb(S(1)),1)
        
      call Pneb_SVD(0,0,npack1, 
     >              A, 
     >              dcpl_mb(U(1)),
     >              dbl_mb(S(1)),
     >              dcpl_mb(V(1)) )


      call Pneb_w_copy_dagger(0,0,dcpl_mb(V(1)),dcpl_mb(Vt(1)))

*     **** find max_sigma ****
      max_sigma = -1.0d200
      indx=S(1)
      do i=1,neall*nbrill
        tmp=dabs(dbl_mb(indx))
        if (tmp.gt.max_sigma) max_sigma=tmp
        indx=indx+1
      end do

*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*c_electron_eorbit(A)
       
*     **** deallocate tmp space ****     
      if (.not.Pneb_w_pop_stack(V)) 
     > call errquit("c_geo_start:error popping stack memory",0,MA_ERR)
     
      call nwpw_timing_end(10)
      return
      end



*     *******************************
*     *		            	    *
*     *		c_geodesic_get	    *
*     *			            *
*     *******************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical    value
      integer    npack1,asize
      integer    tmp1(2),tmp2(2),tmp3(2)


*     **** external functions ****
      logical Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack
 
      value=MA_set_auto_verify(.true.)
      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)

*     **** allocate tmp space ****
      value = Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      call Pneb_SCVtrans1(t,0,0,dbl_mb(S(1)),dcpl_mb(Vt(1)),
     > dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp3(1)))
      call Pneb_www_multiply(0,0,dcpl_mb(Vt(1)),dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp2(1)))
      call BGrsm_f_zero(Ynew)
      call Pneb_fwf_multiply(0,0,Yold,dcpl_mb(tmp2(1)),npack1,Ynew)
      call Pneb_fwf_multiplyAdd(0,0,dcpl_mb(U(1)),dcpl_mb(tmp3(1)),
     > npack1, Ynew )

*     **** Orthonormality Check ****

      call Pneb_orthocheckmake(.true.,0,0,npack1,Ynew)     

*     **** deallocate tmp space ****     
      value = Pneb_w_pop_stack(tmp3)
      value = value.and.
     >        Pneb_w_pop_stack(tmp2)
      value = value.and.
     >        Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end

*     ***********************************
*     *					*
*     *		c_geodesic_transport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic_transport(t,Yold,Ynew)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      integer    npack1
      logical    value
      integer    tmp1(2),tmp2(2),tmp3(2)

*     **** external functions ****
      logical Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)

*     **** allocate tmp space ****
      value = Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      call Pneb_SCVtrans2(t,0,0,dbl_mb(S(1)),dcpl_mb(Vt(1)),
     > dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp3(1)))
      call Pneb_www_multiply(0,0,dcpl_mb(Vt(1)),dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp2(1)))
      call Pneb_fwf_multiply(0,0,Yold,dcpl_mb(tmp2(1)),npack1,
     >  Ynew)
      call Pneb_fwf_multiplyAdd(0,0,dcpl_mb(U(1)),dcpl_mb(tmp3(1)),
     > npack1,Ynew)

*     **** deallocate tmp space ****     
      value = Pneb_w_pop_stack(tmp3)
      value = value.and.
     >        Pneb_w_pop_stack(tmp2)
      value = value.and.
     >        Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *					*
*     *		c_geodesic_Gtransport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic_Gtransport(t,Yold,tG)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 tG(*)

#include "mafdecls.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      integer npack1
      integer tmp1(2),tmp2(2),tmp3(2)
      logical value
*     **** external functions ****
      logical Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack

      
      call nwpw_timing_start(10)

      call Cram_max_npack(npack1)

*     **** allocate tmp space ****

      value = Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

*     ***** calculate U^t * G ****

      call Pneb_ffw_multiply(0,0,
     >                  dcpl_mb(U(1)),
     >                  tG,
     >                  npack1,
     >                  dcpl_mb(tmp1(1)))


      call Pneb_w_copy(0,0,
     >           dcpl_mb(tmp1(1)),
     >           dcpl_mb(tmp3(1)))


      call Pneb_SCVtrans3(t,0,0,dbl_mb(S(1)),
     > dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp3(1)))



      call Pneb_www_multiply(0,0,dcpl_mb(Vt(1)),dcpl_mb(tmp1(1)),
     > dcpl_mb(tmp2(1)))


      call Pneb_fwf_multiply(0,0,Yold,dcpl_mb(tmp2(1)),
     >  npack1,tG)


      call Pneb_fwf_multiplyAdd(0,0,dcpl_mb(U(1)),
     > dcpl_mb(tmp3(1)),
     > npack1,tG)


*     **** deallocate tmp space ****     
      value = Pneb_w_pop_stack(tmp3)
      value = value.and.
     >        Pneb_w_pop_stack(tmp2)
      value = value.and.
     >        Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('error popping stack memory',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end
