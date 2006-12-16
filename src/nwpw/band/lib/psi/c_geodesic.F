* $Id: c_geodesic.F,v 1.17 2006-12-16 01:54:00 bylaska Exp $
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
#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      integer V(2),i
      integer nbrill,nex(2),neall,npack1,indx
      real*8 tmp
      
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


*     **** calculate Vt ****
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
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      logical    value
      integer    npack1,nemax
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)


*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,Pneb_nbrillq
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,Pneb_nbrillq
 
      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax = cpsi_ne(1)+cpsi_ne(2)
      nemax = nemax*Pneb_nbrillq()

*     **** allocate tmp space ****
      value = Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_get: out of stack memory',0,MA_ERR)

      call Pneb_SCVtrans1(0,0,t,
     >                    dbl_mb(S(1)),
     >                    dcpl_mb(Vt(1)),
     >                    dcpl_mb(tmp1(1)),
     >                    dcpl_mb(tmp3(1)),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))

      call Pneb_www_Multiply2(0,0,
     >                        one,
     >                        dcpl_mb(Vt(1)),
     >                        dcpl_mb(tmp1(1)),
     >                        zero,
     >                        dcpl_mb(tmp2(1)))

      call Pneb_fwf_Multiply(0,0,
     >                       one,
     >                       Yold,npack1,
     >                       dcpl_mb(tmp2(1)),
     >                       zero,
     >                       Ynew)

      call Pneb_fwf_Multiply(0,0,
     >                       one,
     >                       dcpl_mb(U(1)),npack1,
     >                       dcpl_mb(tmp3(1)),
     >                       one,
     >                       Ynew)

*     **** Orthonormality Check ****
      call Pneb_orthocheckmake(.true.,0,0,npack1,Ynew)     

*     **** deallocate tmp space ****     
      value =           MA_pop_stack(tmpS(2))
      value = value.and.MA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
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
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))
      integer    npack1,nemax
      logical    value
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,Pneb_nbrillq
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,Pneb_nbrillq

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax = cpsi_ne(1)+cpsi_ne(2)
      nemax = nemax*Pneb_nbrillq()

*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)

      call Pneb_SCVtrans2(0,0,t,
     >                   dbl_mb(S(1)),
     >                   dcpl_mb(Vt(1)),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))
      call Pneb_www_Multiply2(0,0,
     >                       one,
     >                       dcpl_mb(Vt(1)),
     >                       dcpl_mb(tmp1(1)),
     >                       zero,
     >                       dcpl_mb(tmp2(1)))

      call Pneb_fwf_Multiply(0,0,
     >                       mone,
     >                       Yold,npack1,
     >                       dcpl_mb(tmp2(1)),
     >                       zero,
     >                       Ynew)

      call Pneb_fwf_Multiply(0,0,
     >                       one,
     >                       dcpl_mb(U(1)),npack1,
     >                       dcpl_mb(tmp3(1)),
     >                       one,
     >                       Ynew)

*     **** deallocate tmp space ****     
      value =           MA_pop_stack(tmpS(2))
      value = value.and.MA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('c_geodesic_transport:error popping stack',0,MA_ERR)

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
      real*8     t
      complex*16 Yold(*)
      complex*16 tG(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))
      integer    npack1,nemax
      logical    value
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,Pneb_nbrillq
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,Pneb_nbrillq

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax = cpsi_ne(1)+cpsi_ne(2)
      nemax = nemax*Pneb_nbrillq()


*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.MA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value)
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)


      call Pneb_ffw_Multiply(0,0,
     >                       dcpl_mb(U(1)),
     >                       tG,npack1,
     >                       dcpl_mb(tmp2(1)))



      call Pneb_SCVtrans3(0,0,t,
     >                   dbl_mb(S(1)),
     >                   dcpl_mb(tmp2(1)),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))

      call Pneb_www_Multiply2(0,0,
     >                       one,
     >                       dcpl_mb(Vt(1)),
     >                       dcpl_mb(tmp1(1)),
     >                       zero,
     >                       dcpl_mb(tmp2(1)))

      call Pneb_fwf_Multiply(0,0,
     >                       mone,
     >                       Yold,npack1,
     >                       dcpl_mb(tmp2(1)),
     >                       zero,
     >                       tG)

      call Pneb_fwf_Multiply(0,0,
     >                       mone,
     >                       dcpl_mb(U(1)),npack1,
     >                       dcpl_mb(tmp3(1)),
     >                       one,
     >                       tG)


*     **** deallocate tmp space ****
      value =           MA_pop_stack(tmpS(2))
      value = value.and.MA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value)
     > call errquit('c_geodesic_transport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end
