*
* $Id$
*
      subroutine cke_init()
      implicit none 

#include "bafdecls.fh"
#include "cke_common.fh"
#include "errquit.fh"

*     **** local variables ****
      integer nfft3d,npack1,G(3)
      integer i,nb,nbrillq
      real*8  gg,gx,gy,gz,kx,ky,kz,E0,A,bb,sigma,x
      logical value
 
      integer tmp1(2),tmp2(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      logical  control_smooth_cutoff
      external control_smooth_cutoff
      real*8   control_smooth_cutoff_values
      external control_smooth_cutoff_values
      integer  c_G_indx, brillioun_nbrillq
      external c_G_indx, brillioun_nbrillq
      real*8   brillioun_k,lattice_wcut,util_erf
      external brillioun_k,lattice_wcut,util_erf

      nbrillq = brillioun_nbrillq()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      G(1)= c_G_indx(1)
      G(2)= c_G_indx(2)
      G(3)= c_G_indx(3)

      value = BA_alloc_get(mt_dbl,npack1*nbrillq,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) 
     > call errquit('cke_init:out of heap memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'tmp1',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('cke_init:out of stack memory',0,MA_ERR)


      do nb=1,nbrillq
         kx = brillioun_k(1,nb)
         ky = brillioun_k(2,nb)
         kz = brillioun_k(3,nb)
         do i = 1,nfft3d
            gx  =  dbl_mb(G(1)+i-1)+kx
            gy  =  dbl_mb(G(2)+i-1)+ky
            gz  =  dbl_mb(G(3)+i-1)+kz
            gg  = (gx*gx + gy*gy + gz*gz)
            dbl_mb(tmp1(1)+i-1) = -0.5d0*gg
         end do

         call Cram_r_pack(nb,dbl_mb(tmp1(1)))
         call Cram_r_Copy(nb,dbl_mb(tmp1(1)),
     >                       dbl_mb(tg_indx+(nb-1)*npack1))
      end do

*     **** Funny decay added to stabalize unitcell optimization and pressures ****
      filter = control_smooth_cutoff()
      if (filter) then
         value = BA_alloc_get(mt_dbl,npack1*nbrillq,'df',df(2),df(1))
         if (.not. value)
     >      call errquit('c_ke_init:out of heap memory',3,MA_ERR)

         E0    = lattice_wcut()
         A     = control_smooth_cutoff_values(1)*E0
         sigma = control_smooth_cutoff_values(2)
         bb    = 2.0d0*A/(sigma*dsqrt(4.0d0*datan(1.0d0)))

         do nb=1,nbrillq
            kx = brillioun_k(1,nb)
            ky = brillioun_k(2,nb)
            kz = brillioun_k(3,nb)
            do i = 1,nfft3d
               gx  =  dbl_mb(G(1)+i-1)+kx
               gy  =  dbl_mb(G(2)+i-1)+ky
               gz  =  dbl_mb(G(3)+i-1)+kz
               gg  = (gx*gx + gy*gy + gz*gz)

               if (gg.gt.(E0-6.0d0*sigma)) then
                  x = (0.5d0*gg-E0)/sigma
                  dbl_mb(tmp1(1)+i-1) = -A*(1.0d0+util_erf(x))
                  dbl_mb(tmp2(1)+i-1) = 1.0d0 + bb*dexp(-x*x)
               else
                  dbl_mb(tmp1(1)+i-1) = 0.0d0
                  dbl_mb(tmp2(1)+i-1) = 1.0d0
               end if
            end do
            call Cram_r_pack(nb,dbl_mb(tmp1(1)))
            call Cram_r_pack(nb,dbl_mb(tmp2(1)))
            call Cram_r_Copy(nb,dbl_mb(tmp2(1)),
     >                          dbl_mb(df(1)+(nb-1)*npack1))
            call Cram_rr_Sum2(nb,dbl_mb(tmp1(1)),
     >                           dbl_mb(tg_indx+(nb-1)*npack1))
         end do
      end if
      
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value)
     > call errquit('cke_init:popping stack memory',0,MA_ERR)
      return
      end

      subroutine cke_end()
      implicit none 
#include "bafdecls.fh"
#include "cke_common.fh"
#include "errquit.fh"

      logical value

      value = BA_free_heap(tg_hndl)
      if (filter) value = value.and.BA_free_heap(df(2))
      if (.not. value)
     >   call errquit('cke_end:error freeing heap',0, MA_ERR)
      return
      end

*     **********************************************
*     *                                            *
*     *             cke                            *
*     *                                            *
*     **********************************************
      subroutine cke(ispin,ne,psi1_tag,psi2_tag)
      implicit none
      integer    ispin,ne(2)
      integer    psi1_tag
      integer    psi2_tag

#include "bafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer npack1,nbrillq,neall
      integer n,shift1,shift2,nshift,nbq

*     **** external functions ****
      integer  Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      external Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      
      call Cram_max_npack(npack1)
      nshift = 2*npack1
      nbrillq = Pneb_nbrillq()
      neall  = ne(1)+ne(2)

      do nbq=1,nbrillq
         shift1 = cpsi_data_get_chnk(psi1_tag,nbq)
         shift2 = cpsi_data_get_chnk(psi2_tag,nbq)
         do n=1,neall
            call Cram_rc_Mul(nbq,dbl_mb(tg_indx+(nbq-1)*npack1),
     >                           dbl_mb(shift1),dbl_mb(shift2))
            shift1 = shift1 + nshift
            shift2 = shift2 + nshift
         end do
      end do
      

      return 
      end

*     **********************************************
*     *                                            *
*     *             cke_orb                        *
*     *                                            *
*     **********************************************
      subroutine cke_orb(nb,orb1,orb2)
      implicit none
      integer    nb
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "bafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer npack1
      integer ne1,shifts
      integer cpsi_ne
      logical cpsi_spin_orbit
      external cpsi_spin_orbit
      external cpsi_ne

      call Cram_max_npack(npack1)
      call Cram_rc_Mul(nb,
     >                 dbl_mb(tg_indx+(nb-1)*npack1),
     >                 orb1,
     >                 orb2)
      if (cpsi_spin_orbit()) then
        ne1=cpsi_ne(1)
        shifts=ne1*npack1
        call Cram_rc_Mul(nb,
     >                 dbl_mb(tg_indx+(nb-1)*npack1),
     >                 orb1(shifts+1),
     >                 orb2(shifts+1))
      end if      
      return 
      end

*     **********************************************
*     *                                            *
*     *             cke_ave                        *
*     *                                            *
*     **********************************************
      subroutine cke_ave(ispin,ne,psi1_tag,ave)
      implicit none
      integer ispin,ne(2)
      integer psi1_tag
      real*8     ave

#include "bafdecls.fh"
#include "cke_common.fh"
#include "errquit.fh"


*     **** local variables ****
      integer npack1,np,nb,nbq,nbrillq,shift1,nshift1
      integer ms,n,n1(2),n2(2),tmp1(2),occ1_tag,shift2
      real*8  sum

*     **** external functions ****
      integer  Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   brillioun_weight
      external Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      external cpsi_data_get_next
      external brillioun_weight
     
      call Parallel_np(np)
      call Cram_max_npack(npack1)
      if (.not. BA_push_get(mt_dcpl,npack1,'tmp1',tmp1(2),tmp1(1)))
     > call errquit('cke_ave:out of stack memory',0,MA_ERR)

      nbrillq = Pneb_nbrillq()
      nshift1 = 2*npack1

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      occ1_tag = cpsi_data_get_next(psi1_tag)

      ave = 0.0d0
      do nbq=1,nbrillq
        shift1 = cpsi_data_get_chnk(psi1_tag,nbq)
        if (occ1_tag.gt.0) then
         shift2 = cpsi_data_get_chnk(occ1_tag,nbq)
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            call Cram_rc_Mul(nbq,
     >                       dbl_mb(tg_indx+(nbq-1)*npack1),
     >                       dbl_mb(shift1),
     >                       dcpl_mb(tmp1(1)))
            call Cram_cc_idot(nbq,dbl_mb(shift1),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)
            ave = ave + sum*brillioun_weight(nbq)*dbl_mb(shift2)
            shift1 = shift1 + nshift1
            shift2 = shift2 + 1
         end do
         end do

        else
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            call Cram_rc_Mul(nbq,
     >                       dbl_mb(tg_indx+(nbq-1)*npack1),
     >                       dbl_mb(shift1),
     >                       dcpl_mb(tmp1(1)))
            call Cram_cc_idot(nbq,dbl_mb(shift1),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)
            ave = ave + sum*brillioun_weight(nbq)
            shift1 = shift1 + nshift1
         end do
         end do
        end if

      end do
      if (np.gt.1) call Parallel_SumAll(ave)
      if (ispin.eq.1) ave = 2.0d0*ave
      ave = -ave

      if (.not.BA_pop_stack(tmp1(2)))
     > call errquit('cke_ave: error popping stack',0,MA_ERR)
      return 
      end 


*     **********************************************
*     *                                            *
*     *             cke_euv                        *
*     *                                            *
*     **********************************************
      subroutine cke_euv(ispin,ne,psi_tag,euv)
      implicit none 
      integer ispin,ne(2)
      integer    psi_tag
      real*8 euv(3,3)

#include "cke_common.fh"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer npack1,nfft3d,nbrillq,G(2,3)
      integer nbq,i,j,ms,n,n1(2),n2(2)
      integer u,v,s
      integer nshift1,shift1,psi_shift
      logical value

      real*8 pi,scal,sum,kx,ky,kz,weight
      real*8 hm(3,3),Aus(3,3)
      integer tmp1(2),tmp2(2)
      integer occ_tag,occ_shift,shift2

*     **** external functions ****
      integer  c_G_indx,Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   brillioun_k,brillioun_weight
      external c_G_indx,Pneb_nbrillq,Pneb_convert_nb,cpsi_data_get_chnk
      external cpsi_data_get_next
      external brillioun_k,brillioun_weight

      real*8   lattice_unitg,lattice_omega,lattice_unita
      external lattice_unitg,lattice_omega,lattice_unita

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      occ_tag = cpsi_data_get_next(psi_tag)

      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      

      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nbrillq = Pneb_nbrillq()
      nshift1  = 2*npack1
      
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)
      value = BA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) 
     > call errquit('cke_euv:out of stack  memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)

      value = BA_push_get(mt_dcpl,npack1,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)

    
*     **** calculate Aus = Sum(k)Sum(n)Sum(G) psi(G+k,n)**2 (G(u)+k(u))(G(s)+k(s)) ****
      call dcopy(9,0.0d0,0,Aus,1)    
      do nbq=1,nbrillq
        shift1 = cpsi_data_get_chnk(psi_tag,nbq)
        weight = brillioun_weight(nbq)
        kx = brillioun_k(1,nbq)
        ky = brillioun_k(2,nbq)
        kz = brillioun_k(3,nbq)
        call dcopy(nfft3d,dbl_mb(c_G_indx(1)),1,dbl_mb(G(1,1)),1)
        call dcopy(nfft3d,dbl_mb(c_G_indx(2)),1,dbl_mb(G(1,2)),1)
        call dcopy(nfft3d,dbl_mb(c_G_indx(3)),1,dbl_mb(G(1,3)),1)
        do i = 1,nfft3d
          dbl_mb(G(1,1)+i-1) =  dbl_mb(G(1,1)+i-1)+kx
          dbl_mb(G(1,2)+i-1) =  dbl_mb(G(1,2)+i-1)+ky
          dbl_mb(G(1,3)+i-1) =  dbl_mb(G(1,3)+i-1)+kz
        end do
        call Cram_r_pack(nbq,dbl_mb(G(1,1)))
        call Cram_r_pack(nbq,dbl_mb(G(1,2)))
        call Cram_r_pack(nbq,dbl_mb(G(1,3)))

*       **** fractional occupation ****
        if (occ_tag.gt.0) then
          shift2 = cpsi_data_get_chnk(occ_tag,nbq)
          do u=1,3
          do s=u,3
            call Cram_rr_Mul(nbq,dbl_mb(G(1,u)),
     >                        dbl_mb(G(1,s)),
     >                        dbl_mb(tmp1(1)))

            if (filter) 
     >         call Cram_rr_Mul2(nbq,dbl_mb(df(1)+(nbq-1)*npack1),
     >                               dbl_mb(tmp1(1)))

            psi_shift = shift1
            occ_shift = shift2
            do ms=1,ispin
            do n=n1(ms),n2(ms)    
              call Cram_rc_Mul(nbq,dbl_mb(tmp1(1)),
     >                         dbl_mb(psi_shift),
     >                         dcpl_mb(tmp2(1)))
              call Cram_cc_dot(nbq,dbl_mb(psi_shift),
     >                         dcpl_mb(tmp2(1)),
     >                         sum)
               Aus(u,s)  = Aus(u,s) + sum*weight*dbl_mb(occ_shift)
               psi_shift = psi_shift + nshift1
               occ_shift = occ_shift + 1
            end do
            end do
          end do
          end do

*       **** no fractional occupation ****
        else
          do u=1,3
          do s=u,3
            call Cram_rr_Mul(nbq,dbl_mb(G(1,u)),
     >                        dbl_mb(G(1,s)),
     >                        dbl_mb(tmp1(1)))

            if (filter) 
     >         call Cram_rr_Mul2(nbq,dbl_mb(df(1)+(nbq-1)*npack1),
     >                               dbl_mb(tmp1(1)))

            psi_shift = shift1
            do ms=1,ispin
            do n=n1(ms),n2(ms)
              call Cram_rc_Mul(nbq,dbl_mb(tmp1(1)),
     >                         dbl_mb(psi_shift),
     >                         dcpl_mb(tmp2(1)))
              call Cram_cc_dot(nbq,dbl_mb(psi_shift),
     >                         dcpl_mb(tmp2(1)),
     >                         sum)
               Aus(u,s)  = Aus(u,s) + sum*weight
               psi_shift = psi_shift + nshift1
            end do
            end do
          end do
          end do
        end if

      end do
      call K1dB_Vector_SumAll(9,Aus)

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
*     *              cke_Precondition        *
*     *                                      *
*     ****************************************
      subroutine cke_Precondition(nbq,npack,neall,psi,gradk)
      implicit none
      integer nbq,npack,neall
      complex*16   psi(npack,neall)
      complex*16 gradk(npack,neall)

#include "bafdecls.fh"
#include "cke_common.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer k,n
      real*8  sum
      real*8  x,cm,dm,Ep,cm2(1)
      integer tmp1(2)

      real*8   lattice_wggcut,control_Ep
      external lattice_wggcut,control_Ep

      integer ispin,ne(2)

      value = BA_push_get(mt_dcpl,npack,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** My preconditioner ****
      ispin = 2
      ne(1) = 1
      ne(2) = 0
      do n=1,neall
         call Cram_cc_dot(nbq,psi(1,n),psi(1,n),x)
         call Cram_rc_Mul(nbq,
     >                    dbl_mb(tg_indx+(nbq-1)*npack),
     >                    psi(1,n),
     >                    dcpl_mb(tmp1(1)))
         call Cram_cc_dot(nbq,psi(1,n),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)
        sum = -sum
        !write(*,*) "x,sum,Ep=",x,sum,Ep,0.5d0*lattice_wggcut()
CDIR$ NOVECTOR
!$OMP DO
        do k=1,npack
           x = -dbl_mb(tg_indx+(nbq-1)*npack+k-1)
           x = x*dconjg(psi(k,n))*psi(k,n)
           x = x/sum

           cm = 27.0d0+(18.0d0+(12.0d0+8.0d0*x)*x)*x
           dm = (cm + 16.0d0* x**4)
           cm = cm/dm

           gradk(k,n) = gradk(k,n)*(cm)
        end do
!$OMP END DO
        !write(*,*) "xmax=",sum
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack memory',0)


      return
      end 


