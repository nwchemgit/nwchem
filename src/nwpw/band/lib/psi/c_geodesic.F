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
      external cpsi_ne,cpsi_neq,cpsi_nbrillq,Pneb_w_size
      integer  cpsi_data_alloc
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
      logical value
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


        
*     ***********************************
*     *                                 *
*     *        c_geodesic_start00       *
*     *                                 *
*     ***********************************
*           
*     Uses - c_geodesic common block
*        
      subroutine c_geodesic_start00(A_tag,max_sigma,dE)
      implicit none
      integer A_tag
      real*8  max_sigma,dE

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value
      integer nb,i,nbrillq,neall,npack1
      integer ashift,ushift,sshift,vshift,vtshift,V_tag
      real*8 tmp

*     **** external functions **** 
      integer  cpsi_nbrillq,cpsi_neq,Pneb_w_size
      external cpsi_nbrillq,cpsi_neq,Pneb_w_size
      integer  cpsi_data_push_stack,cpsi_data_get_chnk
      external cpsi_data_push_stack,cpsi_data_get_chnk
      real*8   c_electron_eorbit00
      external c_electron_eorbit00


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
      dE = 2.0d0*c_electron_eorbit00(A_tag)

*     **** deallocate tmp space ****     
      call cpsi_data_pop_stack(V_tag)

      call nwpw_timing_end(10)
      return
      end



*     ***********************************
*     *                                 *
*     *         c_geodesic_start0       *
*     *                                 *
*     ***********************************
*
*     Uses - c_geodesic common block
*
      subroutine c_geodesic_start0(A_tag,max_sigma,dE_tag)
      implicit none
      integer A_tag
      real*8  max_sigma
      integer dE_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value
      integer nb,i,nbrillq,neall,npack1
      integer ashift,ushift,sshift,vshift,vtshift,V_tag,dE_shift
      real*8 tmp

*     **** external functions **** 
      integer  cpsi_nbrillq,cpsi_neq,Pneb_w_size
      external cpsi_nbrillq,cpsi_neq,Pneb_w_size
      integer  cpsi_data_push_stack,cpsi_data_get_chnk
      external cpsi_data_push_stack,cpsi_data_get_chnk
      !real*8   c_electron_eorbit
      !external c_electron_eorbit


      call nwpw_timing_start(10)

      call Cram_max_npack(npack1)
      nbrillq = cpsi_nbrillq()
      neall   = cpsi_neq(1)+cpsi_neq(2)

      !write(*,*) "FERA0, nb=",npack1,nbrillq,neall,Pneb_w_size(0,1)
      !value = MA_set_auto_verify(.true.)
      !value = MA_verify_allocator_stuff()
      !call MA_summarize_allocated_blocks()
      !write(*,*) "FERA1, nb=",npack1,nbrillq,neall,Pneb_w_size(0,1)


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
      call c_electron_eorbit0_tag(A_tag,dE_tag)
      do nb=1,nbrillq
        dE_shift         = cpsi_data_get_chnk(dE_tag,nb)
        dbl_mb(dE_shift) = 2.0d0*dbl_mb(dE_shift)
      end do


      !dE = 2.0d0*c_electron_eorbit(A_tag)

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



*     *******************************
*     *                             *
*     *         c_geodesic_get0     *
*     *                             *
*     *******************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_get0(nb,t,Yold_tag,Ynew_tag)
      implicit none
      integer    nb
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
      integer npack1,nemax,nbrillq
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
     > call errquit('c_geodesic_get0: out of stack memory',0,MA_ERR)


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
*     *		c_geodesic_transport0	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic_transport0(nb,t,Yold_tag,Ynew_tag)
      implicit none
      integer nb
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
      integer    npack1,nemax,nbrillq
      integer    tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2)
      integer    yoldshift,ynewshift,ushift,sshift,vtshift

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk
      external cpsi_ne,cpsi_nbrillq,cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      nemax   = cpsi_ne(1)+cpsi_ne(2)
      !nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      value =           Pneb_w_push_get(0,1,tmp1)
      value = value.and.Pneb_w_push_get(0,1,tmp2)
      value = value.and.Pneb_w_push_get(0,1,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     > call errquit('c_geodesic_transport:out of stack memory',0,MA_ERR)

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
*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not. value) 
     >call errquit('c_geodesic_transport0:error popping stack',0,MA_ERR)

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



*     *******************************************
*     *                                         *
*     *         c_geodesic_transport_junk0      *
*     *                                         *
*     *******************************************
*
*   Temporary code until BGrsm_list fixed
*     Uses - geodesic common block
*

      subroutine c_geodesic_transport_junk0(t,Yold_tag,Ynew)
      implicit none
      real*8   t(*)
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

         call Pneb_SCVtrans2(0,nb,t(nb),
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
      if (.not.value)
     >call errquit('c_geodesic_transport0:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end




*     *******************************************
*     *                                         *
*     *         c_geodesic_Gtransport_junk0     *
*     *                                         *
*     *******************************************
*
*   Temporary code until BGrsm_list fixed

*     Uses - geodesic common block
*
      subroutine c_geodesic_Gtransport_junk0(t,Yold_tag,tG)
      implicit none
      real*8   t(*)
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
      if (.not.value)
     >call errquit('c_geodesic_transport0:out of stack memory',0,MA_ERR)

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
         call Pneb_SCVtrans3(0,nb,t(nb),
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
      if (.not.value)
     >call errquit('c_geodesic_transport0:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *                                 *
*     *        c_geodesic_get_U_tag     *
*     *                                 *
*     ***********************************
      integer function c_geodesic_get_U_tag()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

      c_geodesic_get_U_tag = U_tag
      return
      end


*     ***********************************
*     *                                 *
*     *       c_geodesic_ABC_init       *
*     *                                 *
*     ***********************************
      subroutine c_geodesic_ABC_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / c_geodesic_ABC_block / ABC,Ermdr

*     **** external functions ****
      logical  Pneb_w_allocate_block
      external Pneb_w_allocate_block

      if (.not.Pneb_w_allocate_block(0,0,3,ABC))
     >   call errquit('c_geodesic_ABC_init:alloc ABC',0,MA_ERR)
      return
      end

*     ***********************************
*     *                                 *
*     *       c_geodesic_ABC_finalize   *
*     *                                 *
*     ***********************************
      subroutine c_geodesic_ABC_finalize()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / c_geodesic_ABC_block / ABC,Ermdr

*     **** external functions ****
      logical  Pneb_w_free
      external Pneb_w_free

      if (.not.Pneb_w_free(ABC))
     > call errquit('c_geodesic_ABC_end:deallocate ABC',0,MA_ERR)
      return
      end



*     ***********************************
*     *                                 *
*     *       c_geodesic_ABC_start      *
*     *                                 *
*     ***********************************
      subroutine c_geodesic_ABC_start(Ermdr0,Y_tag,Z_tag,HY_tag,HZ_tag)
      implicit none
      real*8  Ermdr0
      integer Y_tag,Z_tag,HY_tag,HZ_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / c_geodesic_ABC_block / ABC,Ermdr

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      integer nb,ms,ishift,A,Bdag,C,tmp1(2),nbrillq,npack1
      real*8 etmp,weight

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack
      integer  Pneb_w_size,cpsi_data_get_allptr,cpsi_nbrillq
      external Pneb_w_size,cpsi_data_get_allptr,cpsi_nbrillq
      real*8   c_electron_eorbit,brillioun_weight
      external c_electron_eorbit,brillioun_weight

      Ermdr = Ermdr0
      call Cram_max_npack(npack1)
      nbrillq = cpsi_nbrillq()

      ishift = Pneb_w_size(0,0)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      if (.not.Pneb_w_push_get(0,0,tmp1))
     >   call errquit('c_geodesic_ABC_start:stack',0,MA_ERR)
      
*     **** Copy Hpsi_k to HY then ****
*     **** generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
      call c_electron_get_Hpsi_k(HY_tag)
      do nb=1,nbrillq
         call c_electron_run0(nb,Z_tag)
      end do
      call c_electron_get_Hpsi_k(HZ_tag)


*     **** generate A=<Y|HY>, Bdag=<U|HY>, C=<U|HU> matrices ****
      call Pneb_ffw_hermit_Multiply(0,0,
     >                           dbl_mb(cpsi_data_get_allptr(Y_tag)),
     >                           dbl_mb(cpsi_data_get_allptr(HY_tag)),
     >                           npack1,dcpl_mb(A))
      call Pneb_ffw_Multiply(0,0,dbl_mb(cpsi_data_get_allptr(Z_tag)),
     >                           dbl_mb(cpsi_data_get_allptr(HY_tag)),
     >                           npack1,dcpl_mb(Bdag))
      call Pneb_ffw_hermit_Multiply(0,0,
     >                           dbl_mb(cpsi_data_get_allptr(Z_tag)),
     >                           dbl_mb(cpsi_data_get_allptr(HZ_tag)),
     >                           npack1,dcpl_mb(C))

*     **** compute V'*A*V ****
      call Pneb_www_Multiply3(0,0,one,
     >                        dcpl_mb(A),
     >                        dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                        dcpl_mb(tmp1(1)),
     >                        zero,
     >                        dcpl_mb(A))

*     **** compute B'*V ****
      call Pneb_www_Multiply3(0,0,one,
     >                       dcpl_mb(Bdag),
     >                       dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                       zero,
     >                       dcpl_mb(tmp1(1)))
      call Parallel_shared_vector_copy(.true.,2*ishift,
     >                                 dcpl_mb(tmp1(1)),dcpl_mb(Bdag))

      if (.not.Pneb_w_pop_stack(tmp1))
     >   call errquit('c_geodesic_ABC_start:pop',1,MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *       c_geodesic_ABC_energy     *
*     *                                 *
*     ***********************************
      real*8 function c_geodesic_ABC_energy(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / c_geodesic_ABC_block / ABC,Ermdr

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      logical value
      integer ishift,A,Bdag,C,nemax,nbrillq
      integer tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2),tmp4(2)
      real*8  energy,etmp

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ispin,cpsi_ne,cpsi_nbrillq
      external cpsi_ispin,cpsi_ne,cpsi_nbrillq
      integer  Pneb_w_size,cpsi_data_get_allptr
      external Pneb_w_size,cpsi_data_get_allptr

      ishift = Pneb_w_size(0,0)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      nemax   = cpsi_ne(1) + cpsi_ne(2)
      nbrillq = cpsi_nbrillq()

*     **** push stack space ****
      value =           Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp2)
      value = value.and.Pneb_w_push_get(0,0,tmp3)
      value = value.and.Pneb_w_push_get(0,0,tmp4)
      value = value.and.
     >        BA_push_get(mt_dbl,nbrillq*nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nbrillq*nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) 
     >   call errquit('c_geodesic_ABC_energy:stack',0,MA_ERR)

      call Pneb_SCVtrans1(0,0,t,
     >                    dbl_mb(cpsi_data_get_allptr(S_tag)),
     >                    dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                    dcpl_mb(tmp1(1)),
     >                    dcpl_mb(tmp3(1)),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))

*     **** V*cos*B*sin*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(Bdag),
     >                        dcpl_mb(tmp3(1)),
     >                        zero,
     >                        dcpl_mb(tmp2(1)))
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(tmp2(1)),
     >                        zero,
     >                        dcpl_mb(tmp4(1)))
*     **** V*cos*A*cos*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(A),
     >                        zero,
     >                        dcpl_mb(tmp2(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp1(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))
*     **** V*sin*B'*cos*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(tmp3(1)),
     >                        dcpl_mb(Bdag),
     >                        zero,
     >                        dcpl_mb(tmp2(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp1(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))

*     **** V*sin*C*sin*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(tmp3(1)),
     >                        dcpl_mb(C),
     >                        zero,
     >                        dcpl_mb(tmp2(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp3(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))

      call Pneb_w_trace_weight(0,0,dcpl_mb(tmp4(1)),energy)
      if (cpsi_ispin().eq.1) energy = energy+energy

*     **** pop stack space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Pneb_w_pop_stack(tmp4)
      value = value.and.Pneb_w_pop_stack(tmp3)
      value = value.and.Pneb_w_pop_stack(tmp2)
      value = value.and.Pneb_w_pop_stack(tmp1)
      if (.not.value) 
     >   call errquit('c_geodesic_ABC_energy:pop',1,MA_ERR)

      c_geodesic_ABC_energy = energy + Ermdr
      return
      end

*     ***********************************
*     *                                 *
*     *      c_geodesic_ABC_denergy     *
*     *                                 *
*     ***********************************
      real*8 function c_geodesic_ABC_denergy(t)
      implicit none
      real*8 t
#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic_common.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / c_geodesic_ABC_block / ABC,Ermdr

*     **** local variables ****
      complex*16 zero,one,mone
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
      parameter (mone=(-1.0d0,0.0d0))

      logical value
      integer nemax,nbrillq,ishift,A,Bdag,C
      integer CV(2),SV(2),SSV(2),CSV(2),tmpC(2),tmpS(2)
      integer tmp1(2),tmp4(2)
      real*8  denergy,denergyB,denergyA,denergyC,denergyBdag

*     **** external functions ****
      logical  Pneb_w_push_get,Pneb_w_pop_stack
      external Pneb_w_push_get,Pneb_w_pop_stack
      integer  cpsi_ispin,cpsi_ne,cpsi_nbrillq
      external cpsi_ispin,cpsi_ne,cpsi_nbrillq
      integer  Pneb_w_size,cpsi_data_get_allptr
      external Pneb_w_size,cpsi_data_get_allptr

      ishift = Pneb_w_size(0,0)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      nemax   = cpsi_ne(1) + cpsi_ne(2)
      nbrillq = cpsi_nbrillq()

*     **** push stack space ****
      value =           Pneb_w_push_get(0,0,CV)
      value = value.and.Pneb_w_push_get(0,0,SV)
      value = value.and.Pneb_w_push_get(0,0,SSV)
      value = value.and.Pneb_w_push_get(0,0,CSV)
      value = value.and.Pneb_w_push_get(0,0,tmp1)
      value = value.and.Pneb_w_push_get(0,0,tmp4)
      value = value.and.
     >        BA_push_get(mt_dbl,nbrillq*nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nbrillq*nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) 
     >   call errquit('c_geodesic_ABC_denergy:push',0,MA_ERR)

      call Pneb_SCVtrans1(0,0,t,
     >                    dbl_mb(cpsi_data_get_allptr(S_tag)),
     >                    dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                    dcpl_mb(CV(1)),
     >                    dcpl_mb(SV(1)),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))
      call Pneb_SCVtrans2(0,0,t,
     >                    dbl_mb(cpsi_data_get_allptr(S_tag)),
     >                    dbl_mb(cpsi_data_get_allptr(Vt_tag)),
     >                    dcpl_mb(SSV(1)),
     >                    dcpl_mb(CSV(1)),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))

*     **** -V*S*sin*B*sin*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                       dcpl_mb(Bdag),
     >                       dcpl_mb(SV(1)),
     >                       zero,
     >                       dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(SSV(1)),
     >                        dcpl_mb(tmp1(1)),
     >                        zero,
     >                        dcpl_mb(tmp4(1)))

      call Pneb_w_trace_weight(0,0,dcpl_mb(tmp4(1)),denergyA)

*     ****  V*cos*B*cos*S*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(Bdag),
     >                        dcpl_mb(CSV(1)),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(CV(1)),
     >                        dcpl_mb(tmp1(1)),
     >                        zero,
     >                        dcpl_mb(tmp4(1)))

      call Pneb_w_trace_weight(0,0,dcpl_mb(tmp4(1)),denergyB)
      denergy = denergy + denergyB


*     **** -V*S*sin*A*cos*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(SSV(1)),
     >                        dcpl_mb(A),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,mone,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(CV(1)), 
     >                        zero,
     >                        dcpl_mb(tmp4(1)))
*     **** -V*cos*A*sin*S*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(CV(1)),
     >                        dcpl_mb(A),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,mone,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(SSV(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))

*     **** V*S*cos* B' *cos*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(CSV(1)),
     >                        dcpl_mb(Bdag),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(CV(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))
*     **** -V*sin*B'*sin*S*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(SV(1)),
     >                        dcpl_mb(Bdag),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,mone,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(SSV(1)), 
     >                        one,
     >                        dcpl_mb(tmp4(1)))

*     ****  V*S*cos*C*sin*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(CSV(1)),
     >                        dcpl_mb(C),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(SV(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))
*     ****  V*sin*C*cos*S*V' ****
      call Pneb_www_Multiply2(0,0,one,
     >                        dcpl_mb(SV(1)),
     >                        dcpl_mb(C),
     >                        zero,
     >                        dcpl_mb(tmp1(1)))
      call Pneb_www_Multiply1(0,0,one,
     >                        dcpl_mb(tmp1(1)),
     >                        dcpl_mb(CSV(1)),
     >                        one,
     >                        dcpl_mb(tmp4(1)))

      call Pneb_w_trace_weight(0,0,dcpl_mb(tmp4(1)),denergyC)
      denergy = denergy + denergyC
      if (cpsi_ispin().eq.1) denergy = denergy+denergy


*     **** pop stack space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.BA_pop_stack(tmp4(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.Pneb_w_pop_stack(CSV)
      value = value.and.Pneb_w_pop_stack(SSV)
      value = value.and.Pneb_w_pop_stack(SV)
      value = value.and.Pneb_w_pop_stack(CV)
      if (.not.value) 
     >   call errquit('c_geodesic_ABC_denergy:pop',1,MA_ERR)

      c_geodesic_ABC_denergy = denergy
      return
      end

