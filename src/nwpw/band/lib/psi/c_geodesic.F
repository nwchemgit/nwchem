*
* $Id: c_geodesic.F,v 1.3 2001-12-28 19:04:24 bylaska Exp $
*


*     ***********************************
*     *			                		*
*     *		c_geodesic_init	         	*
*     *					                *
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_init()
      implicit none

#include "mafdecls.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical value
      integer nfft3d,nemax,n,nbrill
      

*     **** external functions **** 
      integer  cpsi_ne,brillioun_nbrillioun
      external cpsi_ne,brillioun_nbrillioun


      call C3dB_nfft3d(1,nfft3d)
      nemax  = cpsi_ne(1)+cpsi_ne(2)
      n      = cpsi_ne(1)
      nbrill = brillioun_nbrillioun()

      value = MA_alloc_get(mt_dcpl,nfft3d*nemax*nbrill,
     >                     'U',U(2),U(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,2*n*n*nbrill,
     >                     'Vt',Vt(2),Vt(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax*nbrill,
     >                     'S',S(2),S(1))
      if (.not. value) 
     > call errquit('c_geodesic_init:out of heap memory',0)

      return
      end

*     ***********************************
*     *									*
*     *		c_geodesic_finalize			*
*     *					                *
*     ***********************************
*
*     Uses - c_geodesic common block
*
      subroutine c_geodesic_finalize()
      implicit none

#include "mafdecls.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      logical value

      value = MA_free_heap(S(2))
      value = value.and.
     >        MA_free_heap(Vt(2))
      value = value.and.
     >        MA_free_heap(U(2))
      if (.not. value) 
     > call errquit(
     >  'c_geodesic_finalize:error freeing of heap memory',0)

      return
      end



*     ***********************************
*     *					                *
*     *		c_geodesic_start		    *
*     *					                *
*     ***********************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_start(A,max_sigma,dE)
      implicit none
      complex*16 A(*)
      real*8     max_sigma,dE

#include "mafdecls.fh"
#include "c_geodesic_common.fh"


*     **** local variables ****
      logical value
      integer nfft3d,neall,nbrill,ispin
      integer nb,ms,i,j,n,nn,ne(2),indx,indxt
      integer V(2)

*     **** external functions **** 
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      real*8   c_electron_eorbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external c_electron_eorbit

      call nwpw_timing_start(10)

      call C3dB_nfft3d(1,nfft3d)
      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)

      neall  = ne(1) + ne(2)
      nbrill = brillioun_nbrillioun()
      n      = ne(1)
      nn     = n*n

*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*n*n*nbrill),'V',V(2),V(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0)

*     **** HomeGrown SVD ****
      call dcopy(nbrill*neall,0.0d0,0,dbl_mb(S(1)),1)
      call Stfl_SVD(ispin,ne,nfft3d,neall,nbrill,n,
     >              A,
     >              dcpl_mb(U(1)),
     >              dbl_mb(S(1)),
     >              dcpl_mb(V(1)))
      
*     **** find max_sigma ****
      max_sigma = 0.0d0
      do i=1,neall*nbrill
        if (dabs(dbl_mb(S(1)+i-1)).gt.max_sigma) 
     >        max_sigma = dbl_mb(S(1)+i-1)
      end do

*     *** calculate Vt ***
      do nb=1,nbrill
         do ms=1,ispin
         do i=1,ne(ms)
         do j=1,ne(ms)
             indx  = (i-1) + (j-1)*n + (ms-1)*nn + (nb-1)*2*nn
             indxt = (j-1) + (i-1)*n + (ms-1)*nn + (nb-1)*2*nn
             dcpl_mb(Vt(1)+indxt) = dconjg(dcpl_mb(V(1)+indx))
         end do
         end do
         end do
      end do

*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*c_electron_eorbit(A)


*     **** deallocate tmp space ****     
      value = value.and.
     >        MA_pop_stack(V(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0)

      call nwpw_timing_end(10)

      return
      end



*     ***************************
*     *		            		*
*     *		c_geodesic_get		*
*     *					        *
*     ***************************
*
*     Uses - c_geodesic common block
*

      subroutine c_geodesic_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical    value
      integer    nfft3d,neall,nbrill,ispin,ne(2)
      integer    nb,ms,j,k,n,nn,shift1,shift2
      real*8     scal1,scal2
      complex*16 zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)


*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      call nwpw_timing_start(10)
      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)
      call C3dB_nfft3d(1,nfft3d)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1) + ne(2)
      n     = ne(1)
      nn    = n*n
      nbrill = brillioun_nbrillioun()

*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*n*n),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0)

      do nb =1,nbrill
      do ms = 1,ispin
        shift1 = (ms-1)*ne(1)
     >         + (nb-1)*neall
        shift2 = (ms-1)*nn
     >         + (nb-1)*2*nn

        do j=1,ne(ms)
          scal1 = dcos(dbl_mb(S(1)+j-1+shift1)*t)
          scal2 = dsin(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+ (ms-1)*nn) 
     >          = scal1 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+ (ms-1)*nn) 
     >          = scal2 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do

        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2),n,
     >              dcpl_mb(tmp1(1)+(ms-1)*nn),n,
     >              zero,
     >              dcpl_mb(tmp2(1)+(ms-1)*nn),n)

        do k=1,ne(ms)
          call zcopy(nfft3d,zero,0,Ynew(1+(k-1+shift1)*nfft3d),1)
          do j=1,ne(ms)
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(tmp2(1)+j-1+(k-1)*n+(ms-1)*nn),
     >                         Yold(1+(j-1+shift1)*nfft3d),
     >                         Ynew(1+(k-1+shift1)*nfft3d))
            call Cram_cc_zaxpy(nb,
     >                         dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn),
     >                         dcpl_mb(U(1)+(j-1+shift1)*nfft3d),
     >                         Ynew(1+(k-1+shift1)*nfft3d))

          end do
        end do

      end do
      end do


*     **** Orthonormality Check ****     

*     call Grsm_gg_trace(npack1,nemax,Ynew,Ynew,sum2)
* 
*     sum1 = dble(nemax)
*     if (dabs(sum2-sum1).gt.1.0d-10) then
*        call Parallel_taskid(taskid)
*         if (taskid.eq.MASTER) 
*    >      write(*,*) 
*    >      "Waring: Gram-Schmidt being performed in geodesic_get"
*        call Grsm_g_MakeOrtho(npack1,psi_ne(1),Ynew)
*        if (psi_ispin().gt.1) then
*          call Grsm_g_MakeOrtho(npack1,psi_ne(2),
*    >                           Ynew(1 + psi_ne(1)*npack1))
*        end if
*     end if

*     **** deallocate tmp space ****     
      value = MA_pop_stack(tmp3(2))
      value = value.and.
     >        MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0)

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
#include "c_geodesic_common.fh"

*     **** local variables ****
      logical    value
      integer    nfft3d,neall,nbrill,ispin,ne(2)
      integer    nb,ms,j,k,n,nn,shift1,shift2
      complex*16 zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      call nwpw_timing_start(10)
      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)

      call C3dB_nfft3d(1,nfft3d)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1) + ne(2)
      n     = ne(1)
      nn    = n*n
      nbrill = brillioun_nbrillioun()

*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*n*n),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*n*n),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0)


      do nb=1,nbrill
      do ms = 1,ispin
        shift1 = (ms-1)*ne(1) 
     >         + (nb-1)*neall
        shift2 = (ms-1)*nn
     >         + (nb-1)*2*nn

        do j=1,ne(ms)
          scal1 = dbl_mb(S(1)+j-1+shift1)
     >            *dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = dbl_mb(S(1)+j-1+shift1)
     >            *dcos(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+(ms-1)*nn) 
     >       = scal1 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn) 
     >       = scal3 * dcpl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do
        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2), n,
     >              dcpl_mb(tmp1(1)+(ms-1)*nn), n,
     >              zero,
     >              dcpl_mb(tmp2(1)+(ms-1)*nn), n)


        do k=1,ne(ms)
          call zcopy(nfft3d,zero,0,Ynew(1+(k-1+shift1)*nfft3d),1)
          do j=1,ne(ms)
            call Cram_cc_zaxpy(nb,
     >              (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n+(ms-1)*nn)),
     >              Yold(1+(j-1+shift1)*nfft3d),
     >              Ynew(1+(k-1+shift1)*nfft3d))
            call Cram_cc_zaxpy(nb,
     >              dcpl_mb(tmp3(1)+j-1+(k-1)*n+(ms-1)*nn),
     >              dcpl_mb(U(1)+(j-1+shift1)*nfft3d),
     >              Ynew(1+(k-1+shift1)*nfft3d))

          end do
        end do

      end do
      end do

*     **** deallocate tmp space ****     
      value = MA_pop_stack(tmp3(2))
      value = value.and.
     >        MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0)

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
#include "c_geodesic_common.fh"


*     **** local variables ****
      logical    value
      integer    nfft3d,neall,nbrill,ne(2),ispin
      integer    nb,ms,j,k,n,nn,shift1,shift2
      complex*16 zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      call nwpw_timing_start(10)

      zero = dcmplx(0.0d0,0.0d0)
      one  = dcmplx(1.0d0,0.0d0)

      call C3dB_nfft3d(1,nfft3d)
      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)
      neall  = ne(1) + ne(2)
      n      = ne(1)
      nn     = n*n
      nbrill = brillioun_nbrillioun()

*     **** allocate tmp space ****
      value = MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(2*nn*nbrill),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0)

*     ***** calculate U^t * G ****
      call Stfl_ggq_dot(ispin,ne,nfft3d,neall,nbrill,n,
     >                  dcpl_mb(U(1)),
     >                  tG,
     >                  dcpl_mb(tmp1(1)))

      call zcopy(2*nn*nbrill,
     >           dcpl_mb(tmp1(1)),1,
     >           dcpl_mb(tmp3(1)),1)

      do nb = 1,nbrill
      do ms = 1,ispin
        shift1 = (ms-1)*ne(1)
     >         + (nb-1)*neall
        shift2 = (ms-1)*nn 
     >         + (nb-1)*2*nn

        do j=1,ne(ms)
          scal1 =       dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = 1.0d0-dcos(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,ne(ms)
           dcpl_mb(tmp1(1)+j-1+(k-1)*n+shift2) 
     >        = scal1 * dcpl_mb(tmp1(1)+j-1+(k-1)*n + shift2)
           dcpl_mb(tmp3(1)+j-1+(k-1)*n+shift2) 
     >        = scal3 * dcpl_mb(tmp3(1)+j-1+(k-1)*n + shift2)
          end do
        end do

        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >              dcpl_mb(  Vt(1)+shift2), n,
     >              dcpl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dcpl_mb(tmp2(1)+shift2), n)

        do k=1,ne(ms)
          do j=1,ne(ms)

            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(tmp2(1)+j-1+(k-1)*n+shift2)),
     >             Yold(1+(j-1+shift1)*nfft3d),
     >               tG(1+(k-1+shift1)*nfft3d))
            call Cram_cc_zaxpy(nb,
     >             (-1.0d0*dcpl_mb(tmp3(1)+j-1+(k-1)*n+shift2)),
     >             dcpl_mb(U(1)+(j-1+shift1)*nfft3d),
     >             tG(1+(k-1+shift1)*nfft3d))
          end do
        end do

      end do
      end do

*     **** deallocate tmp space ****     
      value = MA_pop_stack(tmp3(2))
      value = value.and.
     >        MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0)

      call nwpw_timing_end(10)

      return
      end
