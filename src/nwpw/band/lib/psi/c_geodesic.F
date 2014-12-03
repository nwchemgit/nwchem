* $Id$
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      integer npack1,neall,nemax,nbrillq

*     **** external functions **** 
      integer  cpsi_ne,cpsi_neq,cpsi_nbrillq,Pneb_w_size
      integer  cpsi_data_alloc
      external cpsi_ne,cpsi_neq,cpsi_nbrillq,Pneb_w_size
      external cpsi_data_alloc

      call Cram_max_npack(npack1)
      neall  = cpsi_neq(1)+cpsi_neq(2)
      nemax  = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()
      
      U_tag  = cpsi_data_alloc(nbrillq,neall,2*npack1)
      Vt_tag = cpsi_data_alloc(nbrillq,1,2*Pneb_w_size(0,1))
      S_tag  = cpsi_data_alloc(nbrillq,nemax,1)
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"


      call cpsi_data_dealloc(U_tag)
      call cpsi_data_dealloc(Vt_tag)
      call cpsi_data_dealloc(S_tag)
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

      subroutine c_geodesic_start(A_tag,max_sigma,dE)
      implicit none
      integer A_tag
      real*8  max_sigma,dE

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      integer nb,i,nbrillq,neall,npack1
      integer ashift,ushift,sshift,vshift,vtshift,V_tag
      real*8 tmp
      
*     **** external functions **** 
      integer  cpsi_nbrillq,cpsi_neq,Pneb_w_size
      integer  cpsi_data_push_stack,cpsi_data_get_chnk
      real*8   c_electron_eorbit
      external cpsi_nbrillq,cpsi_neq,Pneb_w_size
      external cpsi_data_push_stack,cpsi_data_get_chnk
      external c_electron_eorbit


      call nwpw_timing_start(10)   

      call Cram_max_npack(npack1)
      nbrillq = cpsi_nbrillq()
      neall   = cpsi_neq(1)+cpsi_neq(2)

*     **** allocate tmp space ****
      V_tag = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,1))

      max_sigma = -1.0d200
      do nb=1,nbrillq 
         ashift  = cpsi_data_get_chnk(A_tag,nb)
         ushift  = cpsi_data_get_chnk(U_tag,nb)
         vshift  = cpsi_data_get_chnk(V_tag,nb)
         sshift  = cpsi_data_get_chnk(S_tag,nb)
         vtshift = cpsi_data_get_chnk(Vt_tag,nb)

*        **** HomeGrown SVD ****
         call Pneb_SVD(0,nb,npack1, 
     >              dbl_mb(ashift),
     >              dbl_mb(ushift),
     >              dbl_mb(sshift),
     >              dbl_mb(vshift) )

*        **** calculate Vt ****
         call Pneb_w_copy_dagger(0,nb,dbl_mb(vshift),dbl_mb(vtshift))

*        **** find max_sigma ****
         do i=1,neall
            tmp=dabs(dbl_mb(sshift))
            if (tmp.gt.max_sigma) max_sigma=tmp
            sshift=sshift+1
         end do

      end do
      call K1dB_MaxAll(max_sigma)

*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*c_electron_eorbit(A_tag)
       
*     **** deallocate tmp space ****     
      call cpsi_data_pop_stack(V_tag)
     
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

      subroutine c_geodesic_get(t,Yold_tag,Ynew_tag)
      implicit none
      real*8     t
      integer    Yold_tag
      integer    Ynew_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      logical value
      integer nb,npack1,nemax,nbrillq
      integer tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
      integer yoldshift,ynewshift,ushift,sshift,vtshift


*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
 
      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_get: out of stack memory',0,MA_ERR)

      do nb=1,nbrillq
         yoldshift = cpsi_data_get_chnk(Yold_tag,nb)
         ynewshift = cpsi_data_get_chnk(Ynew_tag,nb)
         ushift    = cpsi_data_get_chnk(U_tag,nb)
         sshift    = cpsi_data_get_chnk(S_tag,nb)
         vtshift   = cpsi_data_get_chnk(Vt_tag,nb)
         call Pneb_SCVtrans1(0,nb,t,
     >                    dbl_mb(sshift),
     >                    dbl_mb(vtshift),
     >                    dcpl_mb(tmp1(1)),
     >                    dcpl_mb(tmp3(1)),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))

         call Pneb_www_Multiply2(0,nb,
     >                    one,
     >                    dbl_mb(vtshift),
     >                    dcpl_mb(tmp1(1)),
     >                    zero,
     >                    dcpl_mb(tmp2(1)))

         call Pneb_fwf_Multiply(0,nb,
     >                    one,
     >                    dbl_mb(yoldshift),npack1,
     >                    dcpl_mb(tmp2(1)),
     >                    zero,
     >                    dbl_mb(ynewshift))

         call Pneb_fwf_Multiply(0,nb,
     >                    one,
     >                    dbl_mb(ushift),npack1,
     >                    dcpl_mb(tmp3(1)),
     >                    one,
     >                    dbl_mb(ynewshift))

*        **** Orthonormality Check ****
         call Pneb_orthoCheckMake(.true.,0,nb,npack1,dbl_mb(ynewshift))
      end do

*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
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

      subroutine c_geodesic_transport(t,Yold_tag,Ynew_tag)
      implicit none
      real*8   t
      integer  Yold_tag
      integer  Ynew_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))

      logical    value
      integer    nb,npack1,nemax,nbrillq
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
      integer    yoldshift,ynewshift,ushift,sshift,vtshift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)

      do nb=1,nbrillq
         yoldshift = cpsi_data_get_chnk(Yold_tag,nb)
         ynewshift = cpsi_data_get_chnk(Ynew_tag,nb)
         ushift    = cpsi_data_get_chnk(   U_tag,nb)
         sshift    = cpsi_data_get_chnk(   S_tag,nb)
         vtshift   = cpsi_data_get_chnk(  Vt_tag,nb)

         call Pneb_SCVtrans2(0,nb,t,
     >                   dbl_mb(sshift),
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))
         call Pneb_www_Multiply2(0,nb,
     >                   one,
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   zero,
     >                   dcpl_mb(tmp2(1)))

         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(yoldshift),npack1,
     >                   dcpl_mb(tmp2(1)),
     >                   zero,
     >                   dbl_mb(ynewshift))

         call Pneb_fwf_Multiply(0,nb,
     >                   one,
     >                   dbl_mb(ushift),npack1,
     >                   dcpl_mb(tmp3(1)),
     >                   one,
     >                   dbl_mb(ynewshift))
      end do
*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
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

      subroutine c_geodesic_Gtransport(t,Yold_tag,tG_tag)
      implicit none
      real*8   t
      integer  Yold_tag
      integer  tG_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))

      logical  value
      integer  nb,npack1,nemax,nbrillq
      integer  tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
      integer  yoldshift,tGshift,ushift,sshift,vtshift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()


*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value)
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)


      do nb=1,nbrillq
         yoldshift = cpsi_data_get_chnk(Yold_tag,nb)
         tGshift   = cpsi_data_get_chnk(  tG_tag,nb)
         ushift    = cpsi_data_get_chnk(   U_tag,nb)
         sshift    = cpsi_data_get_chnk(   S_tag,nb)
         vtshift   = cpsi_data_get_chnk(  Vt_tag,nb)
         
         call Pneb_ffw_Multiply(0,nb,
     >                   dbl_mb(ushift),
     >                   dbl_mb(tGshift),npack1,
     >                   dcpl_mb(tmp2(1)))
         call Pneb_SCVtrans3(0,nb,t,
     >                   dbl_mb(sshift),
     >                   dcpl_mb(tmp2(1)),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))
         call Pneb_www_Multiply2(0,nb,
     >                   one,
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   zero,
     >                   dcpl_mb(tmp2(1)))
         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(yoldshift),npack1,
     >                   dcpl_mb(tmp2(1)),
     >                   one,
     >                   dbl_mb(tGshift))
         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(ushift),npack1,
     >                   dcpl_mb(tmp3(1)),
     >                   one,
     >                   dbl_mb(tGshift))
      end do

*     **** deallocate tmp space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value)
     > call errquit('c_geodesic_transport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     *******************************************
*     *						*
*     *		c_geodesic_transport_junk	*
*     *						*
*     *******************************************
*
*   Temporary code until BGrsm_list fixed
*     Uses - geodesic common block
*

      subroutine c_geodesic_transport_junk(t,Yold_tag,Ynew)
      implicit none
      real*8   t
      integer Yold_tag
      complex*16  Ynew(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))

      logical    value
      integer    nb,npack1,nemax,nbrillq
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
c      integer    yoldshift,ynewshift,ushift,sshift,vtshift
      integer    yoldshift,ushift,sshift,vtshift
      integer    nbshift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()
      nbshift = nemax*npack1

*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)

      do nb=1,nbrillq
         yoldshift = cpsi_data_get_chnk(Yold_tag,nb)
c         ynewshift = cpsi_data_get_chnk(Ynew_tag,nb)
         ushift    = cpsi_data_get_chnk(   U_tag,nb)
         sshift    = cpsi_data_get_chnk(   S_tag,nb)
         vtshift   = cpsi_data_get_chnk(  Vt_tag,nb)

         call Pneb_SCVtrans2(0,nb,t,
     >                   dbl_mb(sshift),
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))
         call Pneb_www_Multiply2(0,nb,
     >                   one,
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   zero,
     >                   dcpl_mb(tmp2(1)))

         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(yoldshift),npack1,
     >                   dcpl_mb(tmp2(1)),
     >                   zero,
     >                   Ynew(1+(nb-1)*nbshift))

         call Pneb_fwf_Multiply(0,nb,
     >                   one,
     >                   dbl_mb(ushift),npack1,
     >                   dcpl_mb(tmp3(1)),
     >                   one,
     >                   Ynew(1+(nb-1)*nbshift))
      end do
*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('c_geodesic_transport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     *******************************************
*     *						*
*     *		c_geodesic_Gtransport_junk	*
*     *						*
*     *******************************************
*
*   Temporary code until BGrsm_list fixed

*     Uses - geodesic common block
*

      subroutine c_geodesic_Gtransport_junk(t,Yold_tag,tG)
      implicit none
      real*8   t
      integer   Yold_tag
      complex*16 tG(*)
c      integer  tG_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))

      logical  value
      integer  nb,npack1,nemax,nbrillq
      integer  tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
c      integer  yoldshift,tGshift,ushift,sshift,vtshift
      integer  yoldshift,ushift,sshift,vtshift
      integer  nbshift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external Pneb_w_push_get,Pneb_w_pop_stack
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()
      nbshift = nemax*npack1


*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value)
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)


      do nb=1,nbrillq
         yoldshift = cpsi_data_get_chnk(Yold_tag,nb)
c         tGshift   = cpsi_data_get_chnk(  tG_tag,nb)
         ushift    = cpsi_data_get_chnk(   U_tag,nb)
         sshift    = cpsi_data_get_chnk(   S_tag,nb)
         vtshift   = cpsi_data_get_chnk(  Vt_tag,nb)
         
         call Pneb_ffw_Multiply(0,nb,
     >                   dbl_mb(ushift),
     >                   tG(1+(nb-1)*nbshift),npack1,
     >                   dcpl_mb(tmp2(1)))
         call Pneb_SCVtrans3(0,nb,t,
     >                   dbl_mb(sshift),
     >                   dcpl_mb(tmp2(1)),
     >                   dcpl_mb(tmp1(1)),
     >                   dcpl_mb(tmp3(1)),
     >                   dbl_mb(tmpC(1)),
     >                   dbl_mb(tmpS(1)))
         call Pneb_www_Multiply2(0,nb,
     >                   one,
     >                   dbl_mb(vtshift),
     >                   dcpl_mb(tmp1(1)),
     >                   zero,
     >                   dcpl_mb(tmp2(1)))
         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(yoldshift),npack1,
     >                   dcpl_mb(tmp2(1)),
     >                   one,
     >                   tG(1+(nb-1)*nbshift))
         call Pneb_fwf_Multiply(0,nb,
     >                   mone,
     >                   dbl_mb(ushift),npack1,
     >                   dcpl_mb(tmp3(1)),
     >                   one,
     >                   tG(1+(nb-1)*nbshift))
      end do

*     **** deallocate tmp space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value)
     > call errquit('c_geodesic_transport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end

