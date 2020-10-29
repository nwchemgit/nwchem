*
* $Id$
*

*     ***********************************
*     *					*
*     *	          Cram_Init		*	
*     *					*
*     ***********************************
  
      subroutine Cram_Init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "cram_common.fh"


*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value
      integer nfft3d
      integer i,j,k
      integer k1,k2,k3
      integer q,p,indx
      integer nb
      integer nx,ny,nz
      integer nxh,nyh,nzh      
      integer pack(2),nidb
      real*8  kv(3),ecut

*     **** external functions ****
      logical  cloak_masker
      integer  brillioun_nbrillq
      real*8   lattice_ecut,lattice_wcut,brillioun_k
      external cloak_masker
      external brillioun_nbrillq
      external lattice_ecut,lattice_wcut,brillioun_k
      
      maxsize = brillioun_nbrillq()+1

      value = BA_alloc_get(mt_int,maxsize,
     >                     'nidb_list',
     >                      nidb_list(2),
     >                      nidb_list(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxsize,
     >                     'nidb2_list',
     >                      nidb2_list(2),
     >                      nidb2_list(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxsize,
     >                     'nwaveall_list',
     >                      nwaveall_list(2),
     >                      nwaveall_list(1))
      value = value.and.
     >        BA_alloc_get(mt_int,2*maxsize,
     >                     'pack_list',
     >                      pack_list(2),
     >                      pack_list(1))
      if (.not. value) 
     > call errquit('Cram_init: out of heap memory',0,MA_ERR)


      call cloak_init()
      call Parallel3d_taskid_i(taskid)
      call C3dB_nfft3d(1,nfft3d)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      do nb=0,maxsize-1
         value = BA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2),pack(1))
          if (.not. value) 
     >     call errquit('Cram_init: out of heap memory',0, MA_ERR)

         int_mb(pack_list(1)+ 2*nb)   = pack(1)
         int_mb(pack_list(1)+ 2*nb+1) = pack(2)

         if (nb.eq.0) then
            ecut = lattice_ecut()
            kv(1) = 0.0d0
            kv(2) = 0.0d0
            kv(3) = 0.0d0
         else
            ecut = lattice_wcut()
            kv(1) =  brillioun_k(1,nb)
            kv(2) =  brillioun_k(2,nb)
            kv(3) =  brillioun_k(3,nb)
            !kv(1) =  0.0d0
            !kv(2) =  0.0d0
            !kv(3) =  0.0d0
         end if
         call cloak_set(kv,ecut)

         nidb = 0

*        **** k=(k1,k2,k3) **** 
         do k=(-nzh+1),(nzh-1)
         do j=(-nyh+1),(nyh-1)
         do i=(-nxh+1),(nxh-1)
            k1=i
            k2=j
            k3=k
            if (k1.lt.0) k1 = k1 + nx
            if (k2.lt.0) k2 = k2 + ny
            if (k3.lt.0) k3 = k3 + nz

            !call C3dB_ktoqp(1,k3+1,q,p)
            call C3dB_ijktoindexp(1,k1+1,k2+1,k3+1,indx,p)
              
            if (p.eq.taskid) then
               !indx = (q-1)*(nx)*ny + k2*(nx) + k1 + 1
               if (.not.cloak_masker(indx)) then
                  nidb = nidb + 1
                  int_mb(pack(1)+nidb-1) = indx
               end if
            end if
         end do
         end do
         end do

         int_mb(nidb_list(1)+nb) = nidb
      end do

      call c_Balance_Init(maxsize,
     >                  int_mb(nidb_list(1)),
     >                  int_mb(nidb2_list(1)))

      do nb=0,maxsize-1
        int_mb(nwaveall_list(1)+nb) = int_mb(nidb_list(1)+nb)
        call C3dB_ISumAll(int_mb(nwaveall_list(1)+nb))
      end do

*     **** set npack_max - note that npack_max       **** 
*     **** can be associated with different k-points ****
*     **** on different processors                   ****
      npack_max = 0
      do nb=1,maxsize-1
        if (int_mb(nidb2_list(1)+nb).gt.npack_max) then 
          npack_max = int_mb(nidb2_list(1)+nb)
        end if
      end do

      return 
      end


*     ***********************************
*     *					*
*     *	          Cram_end		*	
*     *					*
*     ***********************************

      subroutine Cram_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      logical value
      integer dum2,nb

      value = .true.
      do nb=0,maxsize-1
         dum2 = int_mb(pack_list(1)+2*nb+1)
         value = value.and.BA_free_heap(dum2)
      end do
      value = value.and.BA_free_heap(nidb_list(2))
      value = value.and.BA_free_heap(nidb2_list(2))
      value = value.and.BA_free_heap(nwaveall_list(2))
      value = value.and.BA_free_heap(pack_list(2))

      if (.not. value) 
     >  call errquit('Cram_end:error freeing heap memory',0, MA_ERR)

      call c_Balance_End()
      call cloak_end()
      return
      end



*     ***********************************
*     *					*
*     *	          Cram_c_pack		*	
*     *					*
*     ***********************************

      subroutine Cram_c_pack(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "bafdecls.fh"  
#include "errquit.fh"
#include "cram_common.fh"

      
*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)
      integer nidb,pack

      call nwpw_timing_start(9)
      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call C3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(2*nfft3d,A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)

      do i=1,(nidb)
        A(i) = dcpl_mb(tmp1(1) + int_mb(pack+i-1) - 1)
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not.value) call errquit('Cram_c_pack:error popping stack',0,
     &       MA_ERR)

      call c_Balance_c_balance(nb,A)

      call nwpw_timing_end(9)

      return
      end



*     ***********************************
*     *                                 *
*     *           Cram_c_pack_start     *
*     *                                 *
*     ***********************************

      subroutine Cram_c_pack_start(nb,A,tmp1,request,reqcnt,msgtype)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 tmp1(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i,nfft3d
      integer nidb,pack

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call C3dB_nfft3d(1,nfft3d)

      call dcopy(2*nfft3d,A,1,tmp1,1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)

      do i=1,(nidb)
        A(i) = tmp1(int_mb(pack+i-1))
      end do

      if (control_balance())
     >  call c_Balance_c_balance_start(nb,A,request,reqcnt,msgtype)

      call nwpw_timing_end(9)
      return
      end


*     ***********************************
*     *                                 *
*     *           Cram_c_pack_end       *
*     *                                 *
*     ***********************************

      subroutine Cram_c_pack_end(nb,tmp1,request,reqcnt)
      implicit none
      integer    nb
      complex*16 tmp1(*)
      integer    request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "cram_common.fh"

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      if (control_balance())
     > call c_Balance_c_balance_end(nb,tmp1,request,reqcnt)

      call nwpw_timing_end(9)
      return
      end



*     ***********************************
*     *					*
*     *	          Cram_r_pack	        *	
*     *					*
*     ***********************************

      subroutine Cram_r_pack(nb,A)
      implicit none
      integer nb
      real*8  A(*)

#include "bafdecls.fh"
#include "cram_common.fh"
#include "errquit.fh"

      
*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)
      integer nidb,pack

      call nwpw_timing_start(9)
      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call C3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) 
     >  call errquit('Cram_r_pack:out of stack memory',0, MA_ERR)

      call dcopy(nfft3d,A,1,dbl_mb(tmp1(1)),1)
      call dcopy(nfft3d,0.0d0,0,A,1)
      do i=1,(nidb)
        A(i) = dbl_mb(tmp1(1) + int_mb(pack+i-1) - 1)
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not.value) 
     >  call errquit('Cram_r_pack:error popping stack',0, MA_ERR)

      call c_Balance_t_balance(nb,A)

      call nwpw_timing_end(9)
      return
      end


*     ***********************************
*     *					*
*     *	          Cram_i_pack	        *	
*     *					*
*     ***********************************

      subroutine Cram_i_pack(nb,A)
      implicit none
#include "errquit.fh"
      integer nb
      integer A(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      
*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)
      integer nidb,pack

      call nwpw_timing_start(9)
      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call C3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_int,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) 
     >  call errquit('Cram_i_pack:out of stack memory',0, MA_ERR)

      call icopy(nfft3d,A,1,int_mb(tmp1(1)),1)
      call icopy(nfft3d,0,0,A,1)
      do i=1,(nidb)
        A(i) = int_mb(tmp1(1) + int_mb(pack+i-1) - 1)
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not.value) 
     >  call errquit('Cram_i_pack:error popping stack',1, MA_ERR)

      call c_Balance_i_balance(nb,A)

      call nwpw_timing_end(9)
      return
      end


*     ***********************************
*     *					*
*     *	          Cram_c_unpack		*	
*     *					*
*     ***********************************

      subroutine Cram_c_unpack(nb,A)
      implicit none
#include "errquit.fh"
      integer    nb
      complex*16 A(*)

#include "bafdecls.fh"
#include "cram_common.fh"


*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)
      integer nidb,pack

      call nwpw_timing_start(9)
      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call c_Balance_c_unbalance(nb,A)

      call C3dB_nfft3d(1,nfft3d)
      
      value = BA_push_get(mt_dcpl,(nidb),
     >                    'tmp1',tmp1(2),tmp1(1))
      if (.not.value) 
     >  call errquit('Cram_c_unpack:out of stack memory',0, MA_ERR)

      call dcopy(2*(nidb),A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)
      do i=1,(nidb)
        A(int_mb(pack+i-1)) = dcpl_mb(tmp1(1)+i-1)
      end do
      value = BA_pop_stack(tmp1(2))
      if (.not.value) 
     >  call errquit('Cram_c_unpack:error popping stack',0, MA_ERR)

      call nwpw_timing_end(9)
      return
      end


*     ***********************************
*     *                                 *
*     *      Cram_c_unpack_start        *
*     *                                 *
*     ***********************************

      subroutine Cram_c_unpack_start(nb,A,tmp1,request,reqcnt,msgtype)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 tmp1(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "cram_common.fh"

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)
      call Cram_c_Copy(nb,A,tmp1)
      if (control_balance())
     > call c_Balance_c_unbalance_start(nb,tmp1,request,reqcnt,msgtype)

      call nwpw_timing_end(9)

      return
      end


*     ***********************************
*     *                                 *
*     *           Cram_c_unpack_end     *
*     *                                 *
*     ***********************************

      subroutine Cram_c_unpack_end(nb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer    nb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i,nfft3d
      integer nidb,pack

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)
      nidb = int_mb(nidb_list(1)+nb)
      pack = int_mb(pack_list(1)+2*nb)

      call C3dB_nfft3d(1,nfft3d)

      call nwpw_timing_start(9)
      if (control_balance())
     >  call c_Balance_c_unbalance_end(nb,tmp1,request,reqcnt)

      call dcopy(2*(nidb),tmp1,1,tmp2,1)
      call dcopy(2*nfft3d,0.0d0,0,tmp1,1)
      do i=1,(nidb)
        tmp1(int_mb(pack+i-1)) = tmp2(i)
      end do

      call nwpw_timing_end(9)

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_npack		*	
*     *					*
*     ***********************************

      subroutine Cram_npack(nb,npack)
      implicit none
      integer nb
      integer npack

#include "bafdecls.fh"
#include "cram_common.fh"


      npack = int_mb(nidb2_list(1)+nb)
      return
      end


*     ***********************************
*     *                                 *
*     *        Cram_max_npack           *
*     *                                 *
*     ***********************************

      subroutine Cram_max_npack(npack)
      implicit none
      integer npack

#include "bafdecls.fh"
#include "cram_common.fh"

      npack = npack_max
      return
      end

*     ***********************************
*     *					*
*     *	         Cram_nwave		*	
*     *					*
*     ***********************************

      integer function Cram_nwave(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "cram_common.fh"

      Cram_nwave =  int_mb(nidb2_list(1)+nb)
      return
      end

*     ***********************************
*     *					*
*     *	         Cram_nwave_brdcst      *	
*     *					*
*     ***********************************
      integer function Cram_nwave_brdcst(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "cram_common.fh"

      integer nbq,taskid_k,p
      integer iw

      if (nb.eq.0) then
         iw = int_mb(nidb2_list(1))
      else
         iw = 0
         call Parallel3d_taskid_k(taskid_k)
         call K1dB_ktoqp(nb,nbq,p)
         if (p.eq.taskid_k) iw = int_mb(nidb2_list(1)+nbq)
         call K1dB_ISumAll(iw)
      end if

      Cram_nwave_brdcst = iw
      return
      end


*     ***********************************
*     *					*
*     *	         Cram_nwave_all		*	
*     *					*
*     ***********************************
      integer function Cram_nwave_all(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "cram_common.fh"

      Cram_nwave_all = int_mb(nwaveall_list(1)+nb)
      return
      end

*     ***********************************
*     *					*
*     *	         Cram_nwave_all_brdcst  *	
*     *					*
*     ***********************************
      integer function Cram_nwave_all_brdcst(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "cram_common.fh"

      integer nbq,taskid_k,p
      integer iw

      if (nb.eq.0) then
         iw = int_mb(nwaveall_list(1))
      else
         iw = 0
         call Parallel3d_taskid_k(taskid_k)
         call K1dB_ktoqp(nb,nbq,p)
         if (p.eq.taskid_k) iw = int_mb(nwaveall_list(1)+nbq)
         call K1dB_ISumAll(iw)
      end if

      Cram_nwave_all_brdcst = iw
      return
      end


*     ***********************************
*     *					*
*     *	         Cram_zero		*	
*     *					*
*     ***********************************

      subroutine Cram_zero(nb,zero,pzero) 
      implicit none
      integer  nb
      integer zero,pzero

      integer qzero

*     *********************************************************
*     **** warning this routine assumes a specific packing ****
*     *********************************************************
*     index = (qzero-1)*(nx/2+1)*ny + (j-1)*(nx/2+1) + i 
      zero = 1
      !call C3dB_ktoqp(1,1,qzero,pzero)
      call C3dB_ijktoindexp(1,1,1,1,zero,pzero)

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_cc_nzdot		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_nzdot(nb,ne,A,B,sum)
      implicit none
      integer    nb
      integer    ne
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer n,np
      integer shift

*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc


      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)

      do n=1,ne
        shift = (n-1)*npack_max

        sum(n) = tzdotc(int_mb(nidb2_list(1)+nb),
     >                         A(1+shift),1,
     >                         B(1),      1)
      end do

      if (np.gt.1) then
         call C3dB_Vector_SumAll(2*ne,sum)
      end if

      call nwpw_timing_end(2)

      return
      end 

*     ***********************************
*     *					*
*     *	         Cram_cc_nzdot		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_nzdot2com(nb,ne,ne1,A,B,sum)
      implicit none
      integer    nb
      integer    ne,ne1
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer n,np
      integer shift,shift2

*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc


      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)

      do n=1,ne
        shift = (n-1)*npack_max
        shift2 = ne1*npack_max
        sum(n) = tzdotc(int_mb(nidb2_list(1)+nb),
     >                         A(1+shift),1,
     >                         B(1),      1)
        sum(n) = sum(n) + tzdotc(int_mb(nidb2_list(1)+nb),
     >                         A(1+shift+shift2),1,
     >                         B(1+shift2),1)			    	
      end do

      if (np.gt.1) then
         call C3dB_Vector_SumAll(2*ne,sum)
      end if

      call nwpw_timing_end(2)

      return
      end 


*     ***********************************
*     *                                 *
*     *          Cram_cc_inzdot          *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_inzdot(nb,ne,A,B,sum)
      implicit none
      integer    nb
      integer    ne
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer n,np
      integer shift

*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc


      call nwpw_timing_start(2)

      do n=1,ne
        shift = (n-1)*npack_max

        sum(n) = tzdotc(int_mb(nidb2_list(1)+nb),
     >                         A(1+shift),1,
     >                         B(1),      1)
      end do


      call nwpw_timing_end(2)

      return
      end

*     ***********************************
*     *                                 *
*     *          Cram_cc_inzdotAdd      *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_inzdotAdd(nb,ne,A,B,sum)
      implicit none
      integer    nb
      integer    ne
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer n,np
      integer shift

*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc


      call nwpw_timing_start(2)

      do n=1,ne
        shift = (n-1)*npack_max

        sum(n) = sum(n)+tzdotc(int_mb(nidb2_list(1)+nb),
     >                         A(1+shift),1,
     >                         B(1),      1)
      end do


      call nwpw_timing_end(2)

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_cc_izdot		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_izdot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum

#include "bafdecls.fh"
#include "cram_common.fh"


*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc

      call nwpw_timing_start(2)

      sum = tzdotc(int_mb(nidb2_list(1)+nb),A,1,B,1)

      call nwpw_timing_end(2)

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_cc_izdotAdd      *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_izdotAdd(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer n,np
      integer shift

*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc


      call nwpw_timing_start(2)

      sum = sum+tzdotc(int_mb(nidb2_list(1)+nb),A,1,B,1)

      call nwpw_timing_end(2)

      return
      end



*     *******************************
*     *				    *
*     *	         Cram_cc_zdot	    *	
*     *				    *
*     *******************************

      subroutine Cram_cc_zdot(nb,A,B,tsum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 tsum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np

*     **** external functions ****
      complex*16   tzdotc
      external     tzdotc

      integer i
      call nwpw_timing_start(2)

      call Parallel3d_np_i(np)

      tsum = tzdotc(int_mb(nidb2_list(1)+nb),A,1,B,1)

      if (np.gt.1) call C3dB_Vector_SumAll(2,tsum)

      call nwpw_timing_end(2)

      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function tzdotc(n,x,incx,y,incy)
      implicit none
      integer n,incx,incy
      complex*16 x(*),y(*)
      integer i,ix,iy
      complex*16 tsum
      tsum = dcmplx(0.0d0,0.0d0)
      if ((incx.eq.1).and.(incy.eq.1)) then
       do i=1,n
         tsum = tsum + dconjg(x(i))*y(i)
       end do
      else
       ix=1
       iy=1
       if (incx.lt.0) ix=1+(-n+1)*incx
       if (incy.lt.0) iy=1+(-n+1)*incy
       do i=1,n
         tsum = tsum + dconjg(x(ix))*y(iy)
         ix=ix+incx
         iy=iy+incy
       end do
      end if
      tzdotc = tsum
      return
      end


*     ***********************************
*     *					*
*     *	         Cram_cc_idot		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_idot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     sum

#include "bafdecls.fh"
#include "cram_common.fh"


*     **** local variables ****
c      integer np

*     **** external functions ****
      real*8   ddot
      external ddot


      call nwpw_timing_start(2)


c      call Parallel3d_np_i(np)

      sum = ddot(2*int_mb(nidb2_list(1)+nb),A,1,B,1)

c     if (np.gt.1) call C3dB_SumAll(sum)

      call nwpw_timing_end(2)

      return
      end 




*     ***********************************
*     *					*
*     *	         Cram_cc_dot		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_dot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     sum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   ddot
      external ddot


      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)

      sum = ddot(2*int_mb(nidb2_list(1)+nb),A,1,B,1)

      if (np.gt.1) call C3dB_SumAll(sum)

      call nwpw_timing_end(2)
      return
      end 




*     ***********************************
*     *                                 *
*     *         Cram_ccm_idgemm         *
*     *                                 *
*     ***********************************

      subroutine Cram_ccm_idgemm(nb,n,m,A,B,alpha,beta,hml)
      implicit none
      integer    nb,n,m
      complex*16 A(*)
      complex*16 B(*)
      real*8     alpha,beta
      real*8     hml(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np,npack1


      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)
      npack1 = int_mb(nidb2_list(1)+nb)

      call DGEMM('T','N',n,m,2*npack1,
     >           alpha,
     >           A, 2*npack_max,
     >           B, 2*npack_max,
     >           beta,
     >           hml,n)


      call nwpw_timing_end(2)
      return
      end 



*     ***********************************
*     *                                 *
*     *          Cram_ccm_izgemm         *
*     *                                 *
*     ***********************************

      subroutine Cram_ccm_izgemm(nb,n,m,A,B,alpha,beta,hml)
      implicit none
      integer    nb,n,m
      complex*16 A(*)
      complex*16 B(*)
      complex*16 alpha,beta
      complex*16 hml(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer npack1


      call nwpw_timing_start(2)
      npack1 = int_mb(nidb2_list(1)+nb)

      call ZGEMM('C','N',n,m,npack1,
     >           alpha,
     >           A, npack_max,
     >           B, npack_max,
     >           beta,
     >           hml,n)


      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_ccm_sym_izgemm    *
*     *                                 *
*     ***********************************

      subroutine Cram_ccm_sym_izgemm(nb,n,A,B,alpha,beta,hml)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      complex*16 alpha,beta
      complex*16 hml(n,n)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer j,k
      integer npack1

      call nwpw_timing_start(2)
      npack1 = int_mb(nidb2_list(1)+nb)

      do k=1,n
        call ZGEMM('C','N',k,1,npack1,
     >           alpha,
     >           A, npack_max,
     >           B(1+(k-1)*npack_max), npack_max,
     >           beta,
     >           hml(1,k),k)
      end do

      do k=1,n
      do j=k+1,n
        hml(j,k) = dconjg(hml(k,j))
      end do
      end do
      do k=1,n
        hml(k,k) = dcmplx(dble(hml(k,k)),0.0d0)
      end do

c      do k=1,n
c      do j=1,n
c         write(*,*) "sym_izgemm=",j,k,hml(j,k)
c      end do
c      end do
      

      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *                                 *
*     *          Cram_cmc_zgemm        *
*     *                                 *
*     ***********************************

      subroutine Cram_cmc_zgemm(nb,n,m,A,B,alpha,beta,hml)
      implicit none
      integer    nb,n,m
      complex*16 A(*)
      complex*16 B(*)
      complex*16 alpha,beta
      complex*16 hml(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np,npack1


      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)
      npack1 = int_mb(nidb2_list(1)+nb)

      call ZGEMM('N','N',npack1,n,m,
     >           alpha,
     >           A, npack_max,
     >           hml,n,
     >           beta,
     >           B, npack_max)


      call nwpw_timing_end(2)
      return
      end




*     ***********************************
*     *					*
*     *	         Cram_rr_idot		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_idot(nb,A,B,sum)
      implicit none
      integer    nb
      real*8  A(*)
      real*8  B(*)
      real*8  sum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** external functions ****
      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

      sum = ddot(int_mb(nidb2_list(1)+nb),A(1),1,B(1),1)

      call nwpw_timing_end(2)
      return
      end 





*     ***********************************
*     *					*
*     *	         Cram_rr_dot		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_dot(nb,A,B,sum)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np)

      sum = ddot(int_mb(nidb2_list(1)+nb),A(1),1,B(1),1)

      if (np.gt.1) call C3dB_SumAll(sum)

      call nwpw_timing_end(2)
      return
      end 


*     ***********************************
*     *					*
*     *	         Cram_r_dsum		*	
*     *					*
*     ***********************************

      subroutine Cram_r_dsum(nb,A,sum)
      implicit none
      integer nb
      real*8  A(*)
      real*8  sum

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   dsum
      external dsum

      call nwpw_timing_start(2)
      call Parallel3d_np_i(np)

      sum = dsum(int_mb(nidb2_list(1)+nb),A(1),1)

      if (np.gt.1) call C3dB_SumAll(sum)

      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *					*
*     *	         Cram_c_Copy		*	
*     *					*
*     ***********************************

      subroutine Cram_c_Copy(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      call dcopy(2*(int_mb(nidb2_list(1)+nb)),A,1,B,1)

      return
      end

*     ***********************************
*     *					*
*     *	         Cram_r_Copy		*	
*     *					*
*     ***********************************

      subroutine Cram_r_Copy(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      call dcopy((int_mb(nidb2_list(1)+nb)),A,1,B,1)

      return
      end

*     ***********************************
*     *					*
*     *	         Cram_rc_Copy		*	
*     *					*
*     ***********************************

      subroutine Cram_rc_Copy(nb,A,B)
      implicit none
      integer nb
      real*8     A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
         B(i) = dcmplx(A(i),0.0d0)
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *          Cram_cr_Copy           *       
*     *                                 *
*     ***********************************

      subroutine Cram_cr_Copy(nb,A,B)
      implicit none
      integer nb
      complex*16 A(*)
      real*8     B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
         B(i) = dble(A(i))
      end do
      return
      end


*     ***********************************
*     *					*
*     *	         Cram_c_Zero		*	
*     *					*
*     ***********************************

      subroutine Cram_c_Zero(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      call dcopy(2*(int_mb(nidb2_list(1)+nb)),0.0d0,0,A,1)

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_cc_Sum		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_Sum(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) + B(i)
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_cc_Sum2           *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        B(i) = B(i) + A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_rr_Sum		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_Sum(nb,A,B,C)
      implicit none
      integer  nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) + B(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Cram_rr_Sum2           *
*     *                                 *
*     ***********************************

      subroutine Cram_rr_Sum2(nb,A,B)
      implicit none
      integer  nb
      real*8 A(*)
      real*8 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        B(i) = B(i) + A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_cc_Sub		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_Sub(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) - B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Cram_rr_Sub		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_Sub(nb,A,B,C)
      implicit none
      integer    nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) - B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Cram_rr_Sqrt		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_Sqrt(nb,A,C)
      implicit none
      integer  nb
      real*8   A(*)
      real*8   C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = dsqrt(A(i))
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_rr_Sqrt1          *
*     *                                 *
*     ***********************************

      subroutine Cram_rr_Sqrt1(nb,A)
      implicit none
      integer  nb
      real*8   A(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        A(i) = dsqrt(A(i))
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_cr_Sqr		*	
*     *					*
*     ***********************************

      subroutine Cram_cr_Sqr(nb,A,C)
      implicit none
      integer    nb
      complex*16 A(*)
      real*8     C(*)

#include "bafdecls.fh"
#include "cram_common.fh"


*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = dble(A(i))**2 + dimag(A(i))**2
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_c_SMul		*	
*     *					*
*     ***********************************

      subroutine Cram_c_SMul(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Cram_c_SMul1           *
*     *                                 *
*     ***********************************

      subroutine Cram_c_SMul1(nb,alpha,A)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        A(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_c_ZMul		*	
*     *					*
*     ***********************************

      subroutine Cram_c_ZMul(nb,alpha,A,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_r_SMul		*	
*     *					*
*     ***********************************

      subroutine Cram_r_SMul(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)
      real*8 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i,nn

      nn = int_mb(nidb2_list(1)+nb)

      do i=1,nn
        C(i) = alpha*A(i)
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_r_SMul1            *
*     *                                 *
*     ***********************************

      subroutine Cram_r_SMul1(nb,alpha,A)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        A(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_cc_daxpy		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_daxpy(nb,alpha,A,B)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

      call daxpy(2*(int_mb(nidb2_list(1)+nb)),alpha,A,1,B,1)

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_rr_daxpy		*	
*     *					*
*     ***********************************

      subroutine Cram_rr_daxpy(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)
      real*8 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = C(i) + alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_cc_zaxpy		*	
*     *					*
*     ***********************************

      subroutine Cram_cc_zaxpy(nb,alpha,A,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = C(i) + alpha*A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_rc_Mul            *
*     *					*
*     ***********************************

      subroutine Cram_rc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,( int_mb(nidb2_list(1)+nb))
        C(i) = A(i)*B(i)
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_rc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Cram_rc_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,( int_mb(nidb2_list(1)+nb))
        B(i) = B(i)*A(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Cram_irc_Mul           *
*     *					*
*     ***********************************

      subroutine Cram_irc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) * dcmplx(-dimag(B(i)),dble(B(i)))
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_icc_Mul           *
*     *					*
*     ***********************************

      subroutine Cram_icc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i) * dcmplx(-dimag(B(i)),dble(B(i)))
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Cram_cc_Mul            *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i)*B(i)
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Cram_cc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_Mul2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        B(i) = B(i)*A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Cram_ccr_conjgMul	*	
*     *					*
*     ***********************************
      subroutine Cram_ccr_conjgMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,( int_mb(nidb2_list(1)+nb))
        C(i) = dble(dconjg(A(i))*B(i))
      end do

      return
      end


*     ***********************************
*     *                			*
*     *          Cram_cc_conjgMul       *
*     *                                 *
*     ***********************************

      subroutine Cram_cc_conjgMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = dconjg(A(i))*B(i)
      end do

      return
      end




*     ***********************************
*     *					*
*     *	         Cram_rc_iMul           *
*     *					*
*     ***********************************

      subroutine Cram_rc_iMul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Cram_rr_Mul            *
*     *					*
*     ***********************************

      subroutine Cram_rr_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      real*8     B(*)
      real*8     C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = A(i)*B(i)
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          Cram_rr_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Cram_rr_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      real*8     B(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        B(i) = B(i)*A(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Cram_zccr_Mulitply2    *
*     *                                 *
*     ***********************************
*
*   This routine used by cpsp force routines  computes 
*     C(i) = Imag (alpha*A(i)*dconjg(B(i)))
*
      subroutine Cram_zccr_Multiply2(nb,alpha,A,B,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = dimag(alpha*A(i)*dconjg(B(i)))
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          Cram_zccr_Mulitply2Add *
*     *                                 *
*     ***********************************
*
*   This routine used by cpsp force routines  computes 
*     C(i) = Imag (alpha*A(i)*dconjg(B(i)))
*
      subroutine Cram_zccr_Multiply2Add(nb,alpha,A,B,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

#include "bafdecls.fh"
#include "cram_common.fh"

*     **** local variables ****
      integer i

      do i=1,(int_mb(nidb2_list(1)+nb))
        C(i) = C(i)+dimag(alpha*A(i)*dconjg(B(i)))
      end do

      return
      end



