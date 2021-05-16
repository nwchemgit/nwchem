*
* $Id$
*

*     **************************************
*     *                                    *
*     *             ke_init                *
*     *                                    *
*     **************************************
      subroutine ke_init()
      implicit none 

#include "bafdecls.fh"
#include "errquit.fh"

*     **** parameters - filter ****
c      real*8 eps,ncut
c      parameter (eps=1.0d-6,ncut=151.0d0)
      real*8 ncut

*     **** local variables ****
      integer npack1,nfft3d,G(3)
      integer i
      real*8  gg,g1,ggcut,A,E0,sigma,bb,x
      logical value
 
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      logical  control_smooth_cutoff
      external control_smooth_cutoff
      integer  G_indx
      external G_indx
      real*8   lattice_wcut,control_smooth_cutoff_values,util_erf
      external lattice_wcut,control_smooth_cutoff_values,util_erf
      

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

      logical filter
      integer df(2)
      common / tg2_block / df,filter


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      G(1)= G_indx(1)
      G(2)= G_indx(2)
      G(3)= G_indx(3)

      value = BA_alloc_get(mt_dbl,npack1,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) 
     > call errquit('ke_init:out of heap memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('ke_init:out of stack memory',0, MA_ERR)


      do i=1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))
         dbl_mb(tmp1(1)+i-1) = -0.5d0*gg 
      end do
      call Pack_t_pack(1,dbl_mb(tmp1(1)))
      call Pack_t_Copy(1,dbl_mb(tmp1(1)),dbl_mb(tg_indx))


*     **** Funny decay added to stabalize unitcell optimization and pressures ****
      filter = control_smooth_cutoff()
      if (filter) then
         value = BA_alloc_get(mt_dbl,npack1,'df',df(2),df(1))
         if (.not. value) 
     >      call errquit('ke_init:out of heap memory',3,MA_ERR)

         E0    = lattice_wcut()
         A     = control_smooth_cutoff_values(1)*E0
         sigma = control_smooth_cutoff_values(2)
         bb    = 2.0d0*A/(sigma*dsqrt(4.0d0*datan(1.0d0)))

         do i=1,nfft3d
            gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >            + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >            + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))

            if (gg.gt.(E0-6.0d0*sigma)) then
               x = (0.5d0*gg-E0)/sigma
               dbl_mb(tmp1(1)+i-1) = -A*(1.0d0+util_erf(x))
               dbl_mb(tmp2(1)+i-1) = 1.0d0 + bb*dexp(-x*x)
            else
               dbl_mb(tmp1(1)+i-1) = 0.0d0
               dbl_mb(tmp2(1)+i-1) = 1.0d0
            end if
         end do
         call Pack_t_pack(1,dbl_mb(tmp1(1)))
         call Pack_t_pack(1,dbl_mb(tmp2(1)))
         call Pack_t_Copy(1,dbl_mb(tmp2(1)),dbl_mb(df(1)))
         call Pack_tt_Sum2(1,dbl_mb(tmp1(1)),dbl_mb(tg_indx))
      end if
      
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('ke_init:popping stack memory',0, MA_ERR)
      return
      end


*     ****************************************
*     *                                      *
*     *               ke_end                 *
*     *                                      *
*     ****************************************
      subroutine ke_end()
      implicit none 

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

      logical filter
      integer df(2)
      common / tg2_block / df,filter

      logical value

      value = BA_free_heap(tg_hndl)
      if (filter) value = value.and.BA_free_heap(df(2))
      if (.not. value) 
     >   call errquit('ke_end:error freeing heap',0, MA_ERR)
      return
      end


*     ****************************************
*     *                                      *
*     *               ke                     *
*     *                                      *
*     ****************************************
      subroutine ke(ispin,ne,psi1,psi2)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "bafdecls.fh"

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

*     ****************************************
*     *                                      *
*     *               ke_add                 *
*     *                                      *
*     ****************************************
      subroutine ke_add(ispin,ne,psi1,psi2)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "bafdecls.fh"
cccccccc#include "frac_occ.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

*     **** local variables ****
      integer npack1
      integer n
      
      call Pack_npack(1,npack1)

      do n=1,(ne(1)+ne(2))
         call Pack_tc_MulAdd(1,dbl_mb(tg_indx),psi1(1+(n-1)*npack1),
     >                                         psi2(1+(n-1)*npack1))
      end do
      
      return 
      end


*     ****************************************
*     *                                      *
*     *               ke_ave                 *
*     *                                      *
*     ****************************************
      subroutine ke_ave(ispin,ne,psi1,ave,fractional,occ)
      implicit none
      integer ispin,ne(2)
      complex*16 psi1(*)
      real*8     ave
      logical fractional
      real*8  occ(*)

#include "bafdecls.fh"
#include "errquit.fh"
cccccccc#include "frac_occ.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl


*     **** local variables ****
      integer npack1,np
      integer ms,n,n1(2),n2(2)
      real*8  sum,ave1

      common /eelectron_ejtmp/ sum,ave1


c     complex*16 tmp1(nfft3d)
      integer tmp1(2)
      logical value
     
      call Parallel_np(np)

      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

!$OMP MASTER
      ave1 = 0.0d0
!$OMP END MASTER
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            if (fractional) then
            call Pack_tc_aMul(1,occ(n),
     >                          dbl_mb(tg_indx),
     >                          psi1(1+(n-1)*npack1),
     >                          dcpl_mb(tmp1(1)))
            else
            call Pack_tc_Mul(1,dbl_mb(tg_indx),
     >                       psi1(1+(n-1)*npack1),
     >                       dcpl_mb(tmp1(1)))
            end if
            call Pack_cc_idot(1,psi1(1+(n-1)*npack1),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)

!$OMP MASTER
            ave1 = ave1 + sum
!$OMP END MASTER
         end do
      end do
!$OMP BARRIER
      if (np.gt.1) call Parallel_SumAll(ave1)
      ave = ave1
      if (ispin.eq.1) ave = 2.0d0*ave
      ave = -ave

      value = BA_pop_stack(tmp1(2))
      return 
      end 



*     ****************************************
*     *                                      *
*     *               ke_euv                 *
*     *                                      *
*     ****************************************
      subroutine ke_euv(ispin,ne,psi,euv)
*
* $Id$
*
      implicit none 
      integer ispin,ne(2)
      complex*16 psi(*)
      real*8 euv(3,3)

#include "bafdecls.fh"
#include "errquit.fh"

      logical filter
      integer df(2)
      common / tg2_block / df,filter

*     **** local variables ****
      integer npack1,nfft3d,G(2,3)
      integer i,j,ms,n,n1(2),n2(2),np_i,np_j
      integer u,v,s
      logical value

      real*8 pi,scal,sum,ave
      real*8 hm(3,3),Aus(3,3)
      integer tmp1(2),tmp2(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      real*8   lattice_unitg,lattice_omega,lattice_unita
      external lattice_unitg,lattice_omega,lattice_unita


      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

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

      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack  memory',0, MA_ERR)


      value = BA_push_get(mt_dbl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dcpl,npack1,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


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

        if (filter) call Pack_tt_Mul2(1,dbl_mb(df(1)),dbl_mb(tmp1(1)))
     
        ave = 0.0d0
        do ms=1,ispin
        do n=n1(ms),n2(ms)    
            call Pack_tc_Mul(1,dbl_mb(tmp1(1)),
     >                       psi(1+(n-1)*npack1),
     >                       dcpl_mb(tmp2(1)))
            call Pack_cc_idot(1,psi(1+(n-1)*npack1),
     >                        dcpl_mb(tmp2(1)),
     >                       sum)
             ave = ave + sum
             !Aus(u,s) = Aus(u,s) + sum
        end do
        end do
        if (np_i.gt.1) call D3dB_SumAll(ave)
        if (np_j.gt.1) call D1dB_SumAll(ave)
        Aus(u,s) = Aus(u,s) + ave

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

      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(G(2,3)) 
      value = value.and.BA_pop_stack(G(2,2)) 
      value = value.and.BA_pop_stack(G(2,1))
      if (.not. value) call errquit('error poping stack memory',0,
     &       MA_ERR)
      return
      end


*     ****************************************
*     *                                      *
*     *               ke_Precondition        *
*     *                                      *
*     ****************************************
      subroutine ke_Precondition(npack,neall,psi,gradk)
      implicit none
      integer npack,neall
      complex*16   psi(npack,neall)
      complex*16 gradk(npack,neall)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

      real*8  sum,elocal
      common /eelectron_ejtmp/ sum,elocal

*     **** local variables ****
      logical value
      integer k,n
      real*8  x,cm,dm,Ep,cm2(1)
      integer tmp1(2)

      real*8   lattice_wggcut,control_Ep
      external lattice_wggcut,control_Ep

      integer ispin,ne(2)
c      integer  psi_ispin,psi_ne
c      external psi_ispin,psi_ne

      value = BA_push_get(mt_dcpl,npack,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** My preconditioner ****
      do n=1,neall
         call Pack_tc_Mul(1,dbl_mb(tg_indx),psi(1,n),dcpl_mb(tmp1(1)))
         call Pack_cc_dot(1,psi(1,n),dcpl_mb(tmp1(1)),sum)
CDIR$ NOVECTOR
!$OMP DO
         do k=1,npack
           x = dbl_mb(tg_indx+k-1)
           x = x*dconjg(psi(k,n))*psi(k,n)
           x = x/sum

           cm = 27.0d0+(18.0d0+(12.0d0+8.0d0*x)*x)*x
           dm = (cm + 16.0d0* x**4)
           cm = cm/dm

           gradk(k,n) = gradk(k,n)*(cm)
         end do
!$OMP END DO
      end do


c*     **** Preconditioner of Tassone, Mauri, and Car ****
c      ispin = psi_ispin()
c      ne(1) = psi_ne(1)
c      ne(2) = psi_ne(2)
c      call ke_ave(ispin,ne,gradk,Ep)
c      write(*,*) "E(R):",Ep
c      Ep = control_Ep()-Ep
c      cm = 1.0d0/(Ep)
c      do k=1,npack       
c         x = -dbl_mb(tg_indx+k-1)
c        dm = (x*cm)
c         if (x.gt.Ep) then
c           do n=1,neall
c              gradk(k,n) = gradk(k,n)/dm
c           end do
c         end if
c      end do

 
c*     **** My preconditioner ****
c      ispin = 2
c      ne(1) = 1
c      ne(2) = 0
c      do n=1,neall
cc        call ke_ave(ispin,ne,gradk(1,n),Ep)
cc        write(*,*) "n,E(R)=",n,Ep,control_Ep()-50*Ep
cc        Ep =  control_Ep() - 15*Ep
c        Ep =  control_Ep()
c        cm = 1.0d0/Ep
cCDIR$ NOVECTOR
c!$OMP DO
c        do k=1,npack
c          x = -dbl_mb(tg_indx+k-1)
c          dm = x*cm
c          if (x.gt.Ep) then
c          gradk(k,n) = gradk(k,n)/dm
c          end if
c        end do
c!$OMP END DO
c      end do
c
c*     **** preconditioner #5 ****
c      ispin = 2
c      ne(1) = 1
c      ne(2) = 0
c      do n=1,neall
c        call ke_ave(ispin,ne,gradk(1,n),Ep,.false.,cm2)
c        Ep = 1.5d0*Ep
c        do k=1,npack
c           x = -2.0d0*dbl_mb(tg_indx+k-1)/Ep
c           cm = 27.0d0+(18.0d0+(12.0d0+8.0d0*x)*x)*x
c           dm = (cm + 16.0d0* x**4)
c           cm = (cm/dm)*(2.0d0/Ep)
c           gradk(k,n) = gradk(k,n)*cm
c        end do
c      end do

      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack memory',0)

      return
      end



*     ****************************************
*     *                                      *
*     *               ke_Precondition2       *
*     *                                      *
*     ****************************************
      subroutine ke_Precondition2(npack,neall,residual,Kresidual)
      implicit none
      integer npack,neall
      complex*16 residual(npack,neall)
      complex*16 Kresidual(npack,neall)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for ke.F ****
c     real*8 tg(nfft3d)
      integer tg_indx,tg_hndl
      common / tg_block / tg_indx,tg_hndl

*     **** local variables ****
      logical value
      integer k,n
      real*8  sum
      real*8  x,cm,dm,Ep,cm2(1)
      integer ispin,ne(2)


*     **** preconditioner #5 ****
      ispin = 2
      ne(1) = 1
      ne(2) = 0
      do n=1,neall
        call ke_ave(ispin,ne,residual(1,n),Ep,.false.,cm2)
        Ep = 1.5d0*Ep
        do k=1,npack
           x = -2.0d0*dbl_mb(tg_indx+k-1)/Ep
           cm = 27.0d0+(18.0d0+(12.0d0+8.0d0*x)*x)*x
           dm = (cm + 16.0d0* x**4)
           cm = (cm/dm)*(2.0d0/Ep)
           Kresidual(k,n) = residual(k,n)*cm
        end do
      end do

      return
      end


