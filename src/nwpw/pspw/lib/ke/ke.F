*
* $Id: ke.F,v 1.4 2001-12-27 20:54:55 bylaska Exp $
*
      subroutine ke_init()
      implicit none 

#include "mafdecls.fh"

*     **** local variables ****
      integer npack1,nfft3d,G(3)
      integer i
      real*8  gg
      logical value
 
      integer tmp1(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      G(1)= G_indx(1)
      G(2)= G_indx(2)
      G(3)= G_indx(3)

      value = MA_alloc_get(mt_dbl,npack1,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) call errquit('out of heap memory',0)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)


      do i = 1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))
         dbl_mb(tmp1(1)+i-1) = -0.5d0*gg
      end do

      call Pack_t_pack(1,dbl_mb(tmp1(1)))
      call Pack_t_Copy(1,dbl_mb(tmp1(1)),dbl_mb(tg_indx))
      

      value = MA_pop_stack(tmp1(2))
      return
      end

      subroutine ke_end()
      implicit none 
#include "mafdecls.fh"
*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

      logical value

      value = MA_free_heap(tg_hndl)
      return
      end


      subroutine ke(ispin,ne,psi1,psi2)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

*     **** local variables ****
      integer npack1
      integer n
      
      call Pack_npack(1,npack1)

      do n=1,(ne(1)+ne(2))
         call Pack_tc_Mul(1,dbl_mb(tg_indx),psi1(1+(n-1)*npack1),
     >                                      psi2(1+(n-1)*npack1))

      end do
      

      return 
      end

      subroutine ke_ave(ispin,ne,psi1,ave)
      implicit none
      integer ispin,ne(2)
      complex*16 psi1(*)
      real*8     ave

#include "mafdecls.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl


*     **** local variables ****
      integer npack1,np
      integer ms,n,n1(2),n2(2)
      real*8  sum

c     complex*16 tmp1(nfft3d)
      integer tmp1(2)
      logical value
     
      call Parallel_np(np)
      call Pack_npack(1,npack1)
      value = MA_push_get(mt_dcpl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      ave = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call Pack_tc_Mul(1,dbl_mb(tg_indx),
     >                       psi1(1+(n-1)*npack1),
     >                       dcpl_mb(tmp1(1)))
            call Pack_cc_idot(1,psi1(1+(n-1)*npack1),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)

            ave = ave + sum
         end do
      end do
      if (np.gt.1) call D3dB_SumAll(ave)
      if (ispin.eq.1) ave = 2.0d0*ave
      ave = -ave

      value = MA_pop_stack(tmp1(2))
      return 
      end 




      subroutine ke_euv(ispin,ne,psi,euv)
*
* $Id: ke.F,v 1.4 2001-12-27 20:54:55 bylaska Exp $
*
      implicit none 
      integer ispin,ne(2)
      complex*16 psi(*)
      real*8 euv(3,3)

#include "mafdecls.fh"

*     **** local variables ****
      integer npack1,nfft3d,G(2,3)
      integer i,j,ms,n,n1(2),n2(2)
      integer u,v,s
      logical value

      real*8 pi,scal,sum
      real*8 hm(3,3),Aus(3,3)
      integer tmp1(2),tmp2(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      real*8   lattice_unitg,lattice_omega,lattice_unita
      external lattice_unitg,lattice_omega,lattice_unita


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      



      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) call errquit('out of stack memory',0)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) call errquit('out of stack memory',0)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack  memory',0)

      value = MA_push_get(mt_dbl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

      value = MA_push_get(mt_dcpl,npack1,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(G(1,1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(G(1,2)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(G(1,3)),1)
      call Pack_t_pack(1,dbl_mb(G(1,1)))
      call Pack_t_pack(1,dbl_mb(G(1,2)))
      call Pack_t_pack(1,dbl_mb(G(1,3)))

*     **** calculate Aus = Sum(n)Sum(G) psi(G,n)**2 G(u)G(s) ****
      call dcopy(9,0.0d0,0,Aus,1)    
      do u=1,3
      do s=u,3
        call Pack_tt_Mul(1,dbl_mb(G(1,u)),
     >                     dbl_mb(G(1,s)),
     >                     dbl_mb(tmp1(1)))
        do ms=1,ispin
        do n=n1(ms),n2(ms)    
            call Pack_tc_Mul(1,dbl_mb(tmp1(1)),
     >                       psi(1+(n-1)*npack1),
     >                       dcpl_mb(tmp2(1)))
            call Pack_cc_dot(1,psi(1+(n-1)*npack1),
     >                       dcpl_mb(tmp2(1)),
     >                       sum)
           
             Aus(u,s) = Aus(u,s) + sum
        end do
        end do
      end do
      end do
      do u=1,3
      do s=u+1,3
         Aus(s,u) = Aus(u,s)
      end do
      end do
      if (ispin.eq.1) call dscal(9,2.0d0,Aus,1)

*     *** calculate euv = -Sum(s) hm(s,v)*Aus(u,s)
      call dcopy(9,0.0d0,0,euv,1)
      do v=1,3
      do u=1,3
         do s=1,3
            euv(u,v) = euv(u,v) - Aus(u,s)*hm(s,v) 
         end do
      end do
      end do

      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      value = value.and.MA_pop_stack(G(2,3)) 
      value = value.and.MA_pop_stack(G(2,2)) 
      value = value.and.MA_pop_stack(G(2,1))
      if (.not. value) call errquit('error poping stack memory',0)
      return
      end



      subroutine ke_Precondition(npack,neall,psi,gradk)
      implicit none
      integer npack,neall
      complex*16   psi(npack,neall)
      complex*16 gradk(npack,neall)

#include "mafdecls.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

*     **** local variables ****
      logical value
      integer k,n
      real*8  sum
      real*8  x,cm,dm
      integer tmp1(2)

      real*8   lattice_ecut
      external lattice_ecut

      value = MA_push_get(mt_dcpl,npack,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

      sum = 2*lattice_ecut()
      do n=1,neall
c        call Pack_tc_Mul(1,dbl_mb(tg_indx),psi(1,n),dcpl_mb(tmp1(1)))
c        call Pack_cc_dot(1,psi(1,n),dcpl_mb(tmp1(1)),sum)
         do k=1,npack
           x = dabs(dbl_mb(tg_indx+k-1)/sum)

           cm = 27.0d0+(18.0d0+(12.0d0+8.0d0*x)*x)*x
           dm = (cm + 16.0d0* x**4)
           cm = cm/dm

c          cm = 1.0d0 + 0.4d0*dexp(-0.1d0*x*x) 

           gradk(k,n) = gradk(k,n)*(cm)
         end do
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack memory',0)

      return
      end







