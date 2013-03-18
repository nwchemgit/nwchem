*
* $Id$
*

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_init             *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_init(nion0,nkatm0,
     >                                nprj,nbasis,psp_type,
     >                                lmax0,sigma,
     >                                nprj_max,l_prj,m_prj,b_prj,
     >                                comp_charge_matrix,
     >                                hartree_matrix)
      implicit none
      integer nion0,nkatm0
      integer nprj(*),nbasis(*),psp_type(*),lmax0(*)
      real*8 sigma(*)
      integer nprj_max
      integer l_prj(nprj_max,*),m_prj(nprj_max,*),b_prj(nprj_max,*)
      integer comp_charge_matrix(*)
      integer hartree_matrix(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      logical value
      integer k,l,m,lm,ii,iii,ia,iia,Gindx,Gall(2),G(3),tsize
      integer iasize,i,j,li,mi,lj,mj,bi,bj,nb
      integer i1,j1,li1,mi1,bi1,lj1,mj1,bj1
      integer indx,comp_charge_tag,comp_indx,matr_ptr
      real*8  scal,gg,fourpioveromega,sumall
      real*8  taunt,rs,pi,w,rcut

*     **** external functions ****
      integer  control_code,G_indx,c_G_indx,ion_katm
      external control_code,G_indx,c_G_indx,ion_katm
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      real*8   lattice_omega,ewald_rcut
      external lattice_omega,ewald_rcut
      integer  psi_data_get_ptr,psi_data_get_chnk
      external psi_data_get_ptr,psi_data_get_chnk
      real*8   nwpw_gaunt,lattice_unita,control_rcut
      external nwpw_gaunt,lattice_unita,control_rcut
      logical  control_use_grid_cmp
      external control_use_grid_cmp

      nion         = nion0
      nkatm        = nkatm0


      pi = 4.0d0*datan(1.0d0)
      rcut = control_rcut()
      if (rcut.le.0.0d0) then
         rs = lattice_unita(1,1)**2 
     >      + lattice_unita(2,1)**2 
     >      + lattice_unita(3,1)**2
         rs = dsqrt(rs)
         rcut=rs/pi

         rs = lattice_unita(1,2)**2 
     >      + lattice_unita(2,2)**2 
     >      + lattice_unita(3,2)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w

         rs = lattice_unita(1,3)**2 
     >      + lattice_unita(2,3)**2 
     >      + lattice_unita(3,3)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w
      end if
      sigma_smooth = rcut


*     **** determine nion_paw, nkatm_paw, katm_paw, ****
*     ****           katm_pawtoion, ion_pawtoion,  ****
*     ****           katm_iontopaw, ion_iontopaw     ****
      nion_paw = 0
      do ii=1,nion
         ia = ion_katm(ii)
         if (psp_type(ia).eq.4) then
            nion_paw = nion_paw + 1
         end if
      end do
      nkatm_paw = 0
      do ia=1,nkatm
         if (psp_type(ia).eq.4) then
            nkatm_paw = nkatm_paw + 1
         end if
      end do
      value =           MA_alloc_get(mt_int,nion,"ion_iontopaw",
     >                               ion_iontopaw(2),ion_iontopaw(1))
      value = value.and.MA_alloc_get(mt_int,nion_paw,"ion_pawtoion",
     >                               ion_pawtoion(2),ion_pawtoion(1))
      value = value.and.MA_alloc_get(mt_int,nkatm,"katm_iontopaw",
     >                               katm_iontopaw(2),katm_iontopaw(1))
      value = value.and.MA_alloc_get(mt_int,nkatm_paw,"katm_pawtoion",
     >                               katm_pawtoion(2),katm_pawtoion(1))
      value = value.and.MA_alloc_get(mt_int,nion_paw,"katm_paw",
     >                               katm_paw(2),katm_paw(1))
      value = value.and.MA_alloc_get(mt_dbl,nkatm_paw,"sigma_paw",
     >                               sigma_paw(2),sigma_paw(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",0,MA_ERR)

      iia = 0
      do ia=1,nkatm
         if (psp_type(ia).eq.4) then
            iia = iia + 1
            int_mb(katm_pawtoion(1)+iia-1) = ia
            int_mb(katm_iontopaw(1)+ia-1)  = iia
         else
            int_mb(katm_iontopaw(1)+ia-1)  = -1
         end if
      end do
      iii = 0
      do ii=1,nion
         ia = ion_katm(ii)
         iia = int_mb(katm_iontopaw(1)+ia-1)
         if (psp_type(ia).eq.4) then
            iii = iii + 1
            int_mb(katm_paw(1)+iii-1)     = iia
            int_mb(ion_pawtoion(1)+iii-1) = ii
            int_mb(ion_iontopaw(1)+ii-1)  = iii
         else
            int_mb(ion_iontopaw(1)+ii-1)  = -1
         end if
      end do
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)
         dbl_mb(sigma_paw(1)+iia-1) = sigma(ia)
         write(*,*) "E3-ion-ion=",
     >      1.0d0/(dsqrt(2*pi)*sigma(ia)),sigma(ia)
      end do


*     **** allocate mult_l and lm_size ****
      value = value.and.MA_alloc_get(mt_int,nkatm_paw,
     >                     "mult_l",mult_l(2),mult_l(1))
      value = value.and.
     >        MA_alloc_get(mt_int,nkatm_paw,
     >                     "lm_size",lm_size(2),lm_size(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

      !**** use_grid_cmp routines ****
      use_grid_cmp = control_use_grid_cmp()

      !**** band structure code ****
      if ((control_code().eq.5).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then
         isgamma = .false.
      else
         isgamma = .true.
      end if

      if (isgamma) then
         call D3dB_nfft3d(1,nfft3d)
         call Pack_npack(0,npack0)
         Gindx = G_indx(1)
      else
         call C3dB_nfft3d(1,nfft3d)
         call Cram_npack(0,npack0)
         Gindx = c_G_indx(1)
      end if

      mult_l_max = 0
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)
         int_mb(mult_l(1)+iia-1)  = 2*lmax0(ia)
         int_mb(lm_size(1)+iia-1) = (2*lmax0(ia)+1)**2
         if (mult_l_max.lt.(2*lmax0(ia))) mult_l_max = 2*lmax0(ia)
      end do

      !*** allocate gk_smooth, gk,and glm ***
      lm_size_max = (mult_l_max+1)**2
      value = MA_alloc_get(mt_dbl,npack0,
     >                     "gk_smooth",gk_smooth(2),gk_smooth(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*nkatm_paw,"gk",gk(2),gk(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*lm_size_max,
     >                     "glm",glm(2),glm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlm",Qlm(2),Qlm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmx",Qlmx(2),Qlmx(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmy",Qlmy(2),Qlmy(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmz",Qlmz(2),Qlmz(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "dEmult_Qlm",dEmult_Qlm(2),dEmult_Qlm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "dE_Qlm",dE_Qlm(2),dE_Qlm(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*nfft3d,"Gall",Gall(2),Gall(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmx(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmy(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmz(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dEmult_Qlm(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dE_Qlm(1)),1)
      call dcopy(3*nfft3d,dbl_mb(Gindx),1,dbl_mb(Gall(1)),1)
      G(1) = Gall(1)
      G(2) = Gall(1)+nfft3d
      G(3) = Gall(1)+2*nfft3d
      if (isgamma) then
         call Pack_t_pack(0,dbl_mb(G(1)))
         call Pack_t_pack(0,dbl_mb(G(2)))
         call Pack_t_pack(0,dbl_mb(G(3)))
      else
         call Cram_r_pack(0,dbl_mb(G(1)))
         call Cram_r_pack(0,dbl_mb(G(2)))
         call Cram_r_pack(0,dbl_mb(G(3)))
      end if

      lm = 0
      do l =0,mult_l_max

*        **** define  |k|**l / (2*l+1)!! ****
         scal = 1.0d0/dble(nwpw_doublefactorial(2*l+1))
         if (l.gt.0) then
            do k=1,npack0
               gg = dbl_mb(G(1)+k-1)**2 
     >            + dbl_mb(G(2)+k-1)**2 
     >            + dbl_mb(G(3)+k-1)**2
               dbl_mb(gk(1)+k-1) = scal*dsqrt(gg)**l
            end do
         else
            call dcopy(npack0,scal,0,dbl_mb(gk(1)),1)
         end if

*        **** define glm =  |k|**l * Tlm(k)/ (2*l+1)!! ****
         do m =-l,l
            call Tesseral3_vector_lm(l,m,npack0,
     >                               dbl_mb(G(1)),
     >                               dbl_mb(G(2)),
     >                               dbl_mb(G(3)),
     >                               dbl_mb(glm(1)+lm*npack0))
            do k=1,npack0
               dbl_mb(glm(1)+lm*npack0+k-1) 
     >            = dbl_mb(glm(1)+lm*npack0+k-1)*dbl_mb(gk(1)+k-1)
            end do

            lm = lm + 1
         end do

      end do

*     **** define gk_smooth(k)  = 4*pi * Exp[-k*k*sigma_smooth**2 / 4] ****
      scal            = 0.25d0*sigma_smooth**2
      fourpioveromega = 16.0d0*datan(1.0d0)/lattice_omega()
      do k=1,npack0
         gg = dbl_mb(G(1)+k-1)**2 
     >      + dbl_mb(G(2)+k-1)**2 
     >      + dbl_mb(G(3)+k-1)**2
         dbl_mb(gk_smooth(1)+k-1) = fourpioveromega*dexp(-gg*scal)
      end do

*     **** define gk(k,iia)  = 4*pi * Exp[-k*k*sigma(iia**2 / 4] ****
      do iia=1,nkatm_paw
         do k=1,npack0
            scal = 0.25d0*dbl_mb(sigma_paw(1)+iia-1)**2
            gg   = dbl_mb(G(1)+k-1)**2 
     >           + dbl_mb(G(2)+k-1)**2 
     >           + dbl_mb(G(3)+k-1)**2
            dbl_mb(gk(1)+(iia-1)*npack0+k-1) 
     >           = fourpioveromega*dexp(-gg*scal)
         end do
      end do
      if (.not.MA_pop_stack(Gall(2)))
     >   call errquit("nwpw_compcharge_init: pop stack",0,MA_ERR)



*     *******************************************************
*     *****  define indexing for compcharge evalulations ****
*     *******************************************************

      if (.not.MA_alloc_get(mt_int,nkatm_paw,"nindx_Tndiff",
     >                      nindx_Tndiff(2),nindx_Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      if (.not.MA_alloc_get(mt_int,nkatm_paw,"shift_Tndiff",
     >                      shift_Tndiff(2),shift_Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      indx = 0
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)
         int_mb(shift_Tndiff(1)+iia-1) = indx

         comp_charge_tag = comp_charge_matrix(ia)

         do l=0,int_mb(mult_l(1)+iia-1)
            comp_indx = psi_data_get_ptr(comp_charge_tag,l+1)
            do m=-l,l
               do j=1,nprj(ia)
                  lj = l_prj(j,ia)
                  mj = m_prj(j,ia)
                  bj = b_prj(j,ia)

                  do i=1,nprj(ia)
                     li = l_prj(i,ia)
                     mi = m_prj(i,ia)
                     bi = b_prj(i,ia)
                     taunt = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                    *dbl_mb(comp_indx+(bi-1)+(bj-1)*nbasis(ia))
                     if (dabs(taunt).gt.1.0d-11) then
                        indx = indx + 1
                     end if
                  end do
               end do
            end do
         end do
         
        int_mb(nindx_Tndiff(1)+iia-1)=indx-int_mb(shift_Tndiff(1)+iia-1)
      end do

      value = MA_alloc_get(mt_int,indx,"lm_Tndiff",
     >                     lm_Tndiff(2),lm_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_int,indx,"iprj_Tndiff",
     >                     iprj_Tndiff(2),iprj_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_int,indx,"jprj_Tndiff",
     >                     jprj_Tndiff(2),jprj_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,indx,"coeff_Tndiff",
     >                     coeff_Tndiff(2),coeff_Tndiff(1))
      if (.not.value)
     > call errquit(
     >  "nwpw_compcharge_init: error allocating work arrays",4,MA_ERR)

      indx = 0
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)
        
         comp_charge_tag = comp_charge_matrix(ia)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            comp_indx = psi_data_get_ptr(comp_charge_tag,l+1)
            do m=-l,l

               do j=1,nprj(ia)
                  lj = l_prj(j,ia)
                  mj = m_prj(j,ia)
                  bj = b_prj(j,ia)

                  do i=1,nprj(ia)
                     li = l_prj(i,ia)
                     mi = m_prj(i,ia)
                     bi = b_prj(i,ia)
                     taunt = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                    *dbl_mb(comp_indx+(bi-1)+(bj-1)*nbasis(ia))
                     if (dabs(taunt).gt.1.0d-11) then
                        int_mb(lm_Tndiff(1)+indx)   = lm
                        int_mb(iprj_Tndiff(1)+indx) = i
                        int_mb(jprj_Tndiff(1)+indx) = j
                        dbl_mb(coeff_Tndiff(1)+indx) = taunt
                        indx = indx + 1
                     end if
                  end do
               end do
               lm = lm + 1
            end do
         end do
      end do



*     ************************************************************
*     *****  define indexing for hartree matrix evalulations  ****
*     ************************************************************

      if (.not.MA_alloc_get(mt_int,nkatm_paw,"nindx_hartree",
     >                      nindx_hartree(2),nindx_hartree(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",5,MA_ERR)

      if (.not.MA_alloc_get(mt_int,nkatm_paw,"shift_hartree",
     >                      shift_hartree(2),shift_hartree(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",6,MA_ERR)

      indx = 0
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)
         int_mb(shift_hartree(1)+iia-1) = indx

         nb = nbasis(ia)
         matr_ptr = psi_data_get_chnk(hartree_matrix(ia))

         do j = 1,nprj(ia)
            lj = l_prj(j,ia)
            mj = m_prj(j,ia)
            bj = b_prj(j,ia)
            do i = 1,nprj(ia)
               li = l_prj(i,ia)
               mi = m_prj(i,ia)
               bi = b_prj(i,ia)
               do j1 = 1,nprj(ia)
                  lj1 = l_prj(j1,ia)
                  mj1 = m_prj(j1,ia)
                  bj1 = b_prj(j1,ia)
                  do i1 = 1,nprj(ia)
                     li1 = l_prj(i1,ia)
                     mi1 = m_prj(i1,ia)
                     bi1 = b_prj(i1,ia)
                     do l=0,int_mb(mult_l(1)+iia-1)
                        do m=-l,l
                           taunt =
     >                        nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *nwpw_gaunt(.false.,l,m,li1,mi1,lj1,mj1)
     >                       *dbl_mb(matr_ptr
     >                              + (bi-1)
     >                              + (bj-1)*nb
     >                              + (bi1-1)*nb*nb
     >                              + (bj1-1)*nb*nb*nb
     +                              + l*nb*nb*nb*nb)
                           if (dabs(taunt).gt.1.0d-11) then
                              indx = indx + 1
                           end if
                        end do
                     end do

                  end do
               end do

            end do
         end do
         int_mb(nindx_hartree(1)+iia-1) = indx
     >                                  - int_mb(shift_hartree(1)+iia-1)
      end do


      value = MA_alloc_get(mt_int,indx,"iprj_hartree",
     >                 iprj_hartree(2),iprj_hartree(1))
      value = value.and.
     >     MA_alloc_get(mt_int,indx,"jprj_hartree",
     >                 jprj_hartree(2),jprj_hartree(1))
      value = value.and.
     >     MA_alloc_get(mt_int,indx,"iprj1_hartree",
     >                 iprj1_hartree(2),iprj1_hartree(1))
      value = value.and.
     >     MA_alloc_get(mt_int,indx,"jprj1_hartree",
     >                 jprj1_hartree(2),jprj1_hartree(1))
      value = value.and.
     >     MA_alloc_get(mt_dbl,indx,"coeff_hartree",
     >                 coeff_hartree(2),coeff_hartree(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_init:allocate fail",7,MA_ERR)


      indx = 0
      do iia=1,nkatm_paw
         ia = int_mb(katm_pawtoion(1)+iia-1)

         nb = nbasis(ia)
         matr_ptr = psi_data_get_chnk(hartree_matrix(ia))

         do j = 1,nprj(ia)
            lj = l_prj(j,ia)
            mj = m_prj(j,ia)
            bj = b_prj(j,ia)
            do i = 1,nprj(ia)
               li = l_prj(i,ia)
               mi = m_prj(i,ia)
               bi = b_prj(i,ia)
               do j1 = 1,nprj(ia)
                  lj1 = l_prj(j1,ia)
                  mj1 = m_prj(j1,ia)
                  bj1 = b_prj(j1,ia)
                  do i1 = 1,nprj(ia)
                     li1 = l_prj(i1,ia)
                     mi1 = m_prj(i1,ia)
                     bi1 = b_prj(i1,ia)
                     do l=0,int_mb(mult_l(1)+iia-1)
                        do m=-l,l
                           taunt =
     >                        nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *nwpw_gaunt(.false.,l,m,li1,mi1,lj1,mj1)
     >                       *dbl_mb(matr_ptr
     >                              + (bi-1)
     >                              + (bj-1)*nb
     >                              + (bi1-1)*nb*nb
     >                              + (bj1-1)*nb*nb*nb
     +                              + l*nb*nb*nb*nb)
                           if (dabs(taunt).gt.1.0d-11) then
                              int_mb(iprj_hartree(1)+indx)  = i
                              int_mb(jprj_hartree(1)+indx)  = j
                              int_mb(iprj1_hartree(1)+indx) = i1
                              int_mb(jprj1_hartree(1)+indx) = j1
                              dbl_mb(coeff_hartree(1)+indx) = taunt
                              indx = indx + 1
                           end if
                        end do
                     end do

                  end do
               end do

            end do
         end do
      end do


*     **** initialize the gaussian integrals ****

      return
      end


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_end              *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      logical value
      value =           MA_free_heap(ion_iontopaw(2))
      value = value.and.MA_free_heap(ion_pawtoion(2))
      value = value.and.MA_free_heap(katm_iontopaw(2))
      value = value.and.MA_free_heap(katm_pawtoion(2))
      value = value.and.MA_free_heap(katm_paw(2))
      value = value.and.MA_free_heap(sigma_paw(2))
      value = value.and.MA_free_heap(mult_l(2))
      value = value.and.MA_free_heap(lm_size(2))
      value = value.and.MA_free_heap(gk_smooth(2))
      value = value.and.MA_free_heap(gk(2))
      value = value.and.MA_free_heap(glm(2))
      value = value.and.MA_free_heap(Qlm(2))
      value = value.and.MA_free_heap(Qlmx(2))
      value = value.and.MA_free_heap(Qlmy(2))
      value = value.and.MA_free_heap(Qlmz(2))
      value = value.and.MA_free_heap(dEmult_Qlm(2))
      value = value.and.MA_free_heap(dE_Qlm(2))
      value = value.and.MA_free_heap(shift_Tndiff(2))
      value = value.and.MA_free_heap(nindx_Tndiff(2))
      value = value.and.MA_free_heap(lm_Tndiff(2))
      value = value.and.MA_free_heap(iprj_Tndiff(2))
      value = value.and.MA_free_heap(jprj_Tndiff(2))
      value = value.and.MA_free_heap(coeff_Tndiff(2))

      value = value.and.MA_free_heap(shift_hartree(2))
      value = value.and.MA_free_heap(nindx_hartree(2))
      value = value.and.MA_free_heap(iprj_hartree(2))
      value = value.and.MA_free_heap(jprj_hartree(2))
      value = value.and.MA_free_heap(iprj1_hartree(2))
      value = value.and.MA_free_heap(jprj1_hartree(2))
      value = value.and.MA_free_heap(coeff_hartree(2))

      if (.not.value)
     >   call errquit("nwpw_compcharge_end: freeing heap",0,MA_ERR)
      return
      end 

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_lm_size_max          *
*     *                                               *
*     *************************************************
      integer function nwpw_compcharge_lm_size_max()
      implicit none
#include "nwpw_compcharge.fh"

      nwpw_compcharge_lm_size_max = lm_size_max
      return
      end

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_mult_l               *
*     *                                               *
*     *************************************************
      integer function nwpw_compcharge_mult_l(ia)
      implicit none
      integer ia

#include "mafdecls.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iia

      iia = int_mb(katm_iontopaw(1)+ia-1)

      if (iia.eq.-1) then
         nwpw_compcharge_mult_l = -1
      else
         nwpw_compcharge_mult_l = int_mb(mult_l(1)+iia-1)
      end if
      return
      end

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_Qlm                  *
*     *                                               *
*     *************************************************
      real*8 function nwpw_compcharge_Qlm(ms,ii,l,m)
      implicit none
      integer ms,ii,l,m

#include "mafdecls.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer lm,iii,iia,indx
      real*8 tmp

      tmp = 0.0d0
      iii = int_mb(ion_iontopaw(1)+ii-1)
      if (iii.ne.-1) then
         iia = int_mb(katm_paw(1)+iii-1)

         if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
            lm = l*(l+1) + m
            indx  = (iii-1)*2*lm_size_max + (ms-1)*lm_size_max + lm
            tmp = dbl_mb(Qlm(1)+indx)
         end if
      end if

      nwpw_compcharge_Qlm = tmp
      return
      end


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_gen_Qlm          *
*     *                                               *
*     *************************************************

*     This routine computes Qlm for atom ii.

      subroutine nwpw_compcharge_gen_Qlm(ii,ia,ispin,nprj,wmatrix)
      implicit none
      integer ii,ia
      integer ispin,nprj
      real*8  wmatrix(nprj,nprj,ispin)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iii,iia,n,k,ms,lm,iprj,jprj,shift,indx,indx1
      real*8  scal,coeff,w

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      iii = int_mb(ion_iontopaw(1) +ii-1)
      iia = int_mb(katm_iontopaw(1)+ia-1)

      indx = (iii-1)*2*lm_size_max
      call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)+indx),1)
      scal = 1.0d0/lattice_omega()

      shift = int_mb(shift_Tndiff(1)+iia-1)
      do k=1,int_mb(nindx_Tndiff(1)+iia-1)
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = scal*dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
            indx1 = indx + (ms-1)*lm_size_max + lm
            dbl_mb(Qlm(1)+indx1) = dbl_mb(Qlm(1)+indx1) 
     >                           + coeff*wmatrix(iprj,jprj,ms)
         end do
      end do
      call D1dB_Vector_SumAll(ispin*lm_size_max,dbl_mb(Qlm(1)+indx))

      return
      end 

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_add_dE_Qlm           *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_add_dE_Qlm(ispin,ii,l,m,deqlm)
      implicit none
      integer ispin,ii,l,m
      real*8 deqlm

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer ms,iii,iia,indx,lm

      iii = int_mb(ion_iontopaw(1)+ii-1)
      if (iii.ne.-1) then
         iia = int_mb(katm_paw(1)+iii-1)
         if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
            lm = l*(l+1) + m
            do ms =1,ispin
               indx  = (iii-1)*2*lm_size_max + (ms-1)*lm_size_max + lm
               dbl_mb(dE_Qlm(1)+indx) = deqlm
            end do
         end if
      end if
      return
      end

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_reduce_dE_Qlm        *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_reduce_dE_Qlm()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dE_Qlm(1)))

      return
      end

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_zero_dE_Qlm          *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_zero_dE_Qlm()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dE_Qlm(1)),1)
      return
      end


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_gen_sw2          *
*     *                                               *
*     *************************************************

*     This routine computes sw2 from Qlm for atom ii. Where sw2
*
*      dE_Qlm/dpsi_i(r) = Sum(I,lm) dE/dQlm * dQlm/dpsi_i(r) = sum(iprj) |iprj> * sw2(iprj,i)
*     
      subroutine nwpw_compcharge_gen_sw2(ii,ia,ispin,ne,nprj,sw1,sw2)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iii,iia,n,k,ms,lm,iprj,jprj,shift,indx,indx1,n1(2),n2(2)
      real*8  scal,coeff,w

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      iii = int_mb(ion_iontopaw(1) +ii-1)
      iia = int_mb(katm_iontopaw(1)+ia-1)

      indx = (iii-1)*2*lm_size_max

      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      !scal = 1.0d0/lattice_omega()

      shift = int_mb(shift_Tndiff(1)+iia-1)
      do k=1,int_mb(nindx_Tndiff(1)+iia-1)
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
            indx1 = indx + (ms-1)*lm_size_max + lm
            !w     = dbl_mb(dEmult_Qlm(1)+indx1)*coeff
            w     = dbl_mb(dE_Qlm(1)+indx1)*coeff

            do n=n1(ms),n2(ms)
               sw2(n,iprj) = sw2(n,iprj) + sw1(n,jprj)*w
            end do

         end do
      end do

      return
      end 


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_gen_Qlmxyz       *
*     *                                               *
*     *************************************************

*     This routine computes Qlmx,Qlmy,Qlmz for atom ii.

      subroutine nwpw_compcharge_gen_Qlmxyz(ii,ia,ispin,ne,nprj,sw1,
     >                                      sw1x,sw1y,sw1z)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw1x(ne(1)+ne(2),nprj)
      real*8  sw1y(ne(1)+ne(2),nprj)
      real*8  sw1z(ne(1)+ne(2),nprj)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iii,iia,n,k,ms,lm,iprj,jprj,shift,indx,indx1,n1(2),n2(2)
      real*8  scal,coeff,wx,wy,wz

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      iii = int_mb(ion_iontopaw(1) +ii-1)
      iia = int_mb(katm_iontopaw(1)+ia-1)

      indx = (iii-1)*2*lm_size_max
      call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlmx(1)+indx),1)
      call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlmy(1)+indx),1)
      call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlmz(1)+indx),1)
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      scal = 1.0d0/lattice_omega()

      shift = int_mb(shift_Tndiff(1)+iia-1)
      do k=1,int_mb(nindx_Tndiff(1)+iia-1)
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = scal*dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
            wx = 0.0d0
            wy = 0.0d0
            wz = 0.0d0
            do n=n1(ms),n2(ms)
               wx = wx + sw1x(n,iprj)*sw1(n,jprj) 
     >                 + sw1(n,iprj)*sw1x(n,jprj)
               wy = wy + sw1y(n,iprj)*sw1(n,jprj) 
     >                 + sw1(n,iprj)*sw1y(n,jprj)
               wz = wz + sw1z(n,iprj)*sw1(n,jprj) 
     >                 + sw1(n,iprj)*sw1z(n,jprj)
            end do
c            indx1 = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max + lm
            indx1 = indx + (ms-1)*lm_size_max + lm
            dbl_mb(Qlmx(1)+indx1) = dbl_mb(Qlmx(1)+indx1) + coeff*wx
            dbl_mb(Qlmy(1)+indx1) = dbl_mb(Qlmy(1)+indx1) + coeff*wy
            dbl_mb(Qlmz(1)+indx1) = dbl_mb(Qlmz(1)+indx1) + coeff*wz
         end do
      end do
      call D1dB_Vector_SumAll(2*lm_size_max,dbl_mb(Qlmx(1)+indx))
      call D1dB_Vector_SumAll(2*lm_size_max,dbl_mb(Qlmy(1)+indx))
      call D1dB_Vector_SumAll(2*lm_size_max,dbl_mb(Qlmz(1)+indx))

      return
      end


c*     *************************************************
c*     *                                               *
c*     *              nwpw_compcharge_Qlm_sub          *
c*     *                                               *
c*     *************************************************
c      subroutine nwpw_compcharge_Qlm_sub(ispin,ne,nprj,Tndiff,sw1,qlm)
c      implicit none
c      integer ispin,ne(2),nprj
c      real*8  Tndiff(nprj,nprj)
c      real*8  sw1(ne(1)+ne(2),nprj)
c      real*8  qlm(*)
c
c*     **** local variables ****
c      integer ms,i,j,n,n1(2),n2(2)
c      real*8 psum
c
c      n1(1) = 1
c      n1(2) = ne(1)+1
c      n2(1) = ne(1)
c      n2(2) = ne(1)+ne(2)
c      do ms=1,ispin
c         qlm(ms) = 0.0d0
c      end do
c      do j=1,nprj
c         do i=1,nprj
c            do ms=1,ispin
c               psum = 0.0d0
c               do n=n1(ms),n2(ms)
c                  psum = psum + sw1(n,i)*sw1(n,j)
c               end do
c               qlm(ms) = qlm(ms) + Tndiff(i,j)*psum
c            end do
c         end do
c      end do
c      return
c      end


*     *************************************************
*     *                                               *
*     *         nwpw_compcharge_gen_rgaussian         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_rgaussian(l,sigma,nr,r,gl)
      implicit none
      integer  l
      integer  nr
      real*8 sigma
      double precision r(*)
      double precision gl(*)

      !*** local variables ***
      integer i
      real*8  c,pi

*     ***** external functions ****
      real*8   nwpw_double_factorial
      external nwpw_double_factorial

      pi = 4.0d0*datan(1.0d0)
      c  = 2.0d0**(l+2)/(dsqrt(pi)
     >    *(nwpw_double_factorial(2*l+1))*sigma**(2*l+3))

      !*** this fixes possible underflow error ***
      do i=1,nr
         gl(i) = 0.0d0
      end do

      do i=1,nr
c        gl(i) = c*(r(i)**l)*dexp(-(r(i)/sigma)**2)
        if (dabs(r(i)) .lt. (8.0d0*sigma))
     >     gl(i) = c*r(i)**l*dexp(-(r(i)/sigma)**2)
      end do

      return
      end

*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dn_cmp_smooth         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_dn_cmp_smooth(ms,dng_smooth)
      implicit none
      integer ms
      complex*16 dng_smooth(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall


      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp_smooth:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_smooth,1)
      if (isgamma) then
      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l
               
               indx = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx)
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)

      end do
      call Pack_tc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      else

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l

               indx = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx)
               if (mod(l,2).eq.0) then
                  call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)

      end do
      call Cram_rc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      end if


      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp_smooth:stack",1,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_dn_cmp2         *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transforms of dn_cmp and dn_cmp_smooth
 
      subroutine nwpw_compcharge_gen_dn_cmp2(ispin,dng_cmp,dng_smooth)
      implicit none
      integer    ispin
      complex*16 dng_cmp(*)
      complex*16 dng_smooth(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall

      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_cmp,1)
      call dcopy(2*npack0,0.0d0,0,dng_smooth,1)

      if (isgamma) then

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l
               
               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)
         call Pack_tc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do
      call Pack_tc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      else

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l

               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
               if (mod(l,2).eq.0) then
                  call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)
         call Cram_rc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do
      call Cram_rc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      end if

      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",1,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_dn_cmp          *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transform of dn_cmp.
 
      subroutine nwpw_compcharge_gen_dn_cmp(ispin,dng_cmp)
      implicit none
      integer    ispin
      complex*16 dng_cmp(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall


      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_cmp,1)

      if (isgamma) then

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l
               
               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_tc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do

      else

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            do m=-l,l

               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
               if (mod(l,2).eq.0) then
                  call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Cram_rc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)
      end do

      end if

      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",1,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_dn_cmp2_zv      *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transforms of dn_cmp and dn_cmp_smooth
 
      subroutine nwpw_compcharge_gen_dn_cmp2_zv(ispin,zv,
     >                                          dng_cmp,dng_smooth)
      implicit none
      integer    ispin
      real*8 zv(*)
      complex*16 dng_cmp(*)
      complex*16 dng_smooth(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_cmp,1)
      call dcopy(2*npack0,0.0d0,0,dng_smooth,1)

      if (isgamma) then

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)

         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            if (l.eq.0) then
               q00 = -zv(ia)/dsqrt(fourpi)
            else
               q00 = 0.0d0
            end if
            do m=-l,l
               
               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2) + q00
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)
         call Pack_tc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do
      call Pack_tc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      else

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia  = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            if (l.eq.0) then
               q00 = -zv(ia)/dsqrt(fourpi)
            else
               q00 = 0.0d0
            end if
            do m=-l,l

               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2) + q00
               if (mod(l,2).eq.0) then
                  call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)
         call Cram_rc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do
      call Cram_rc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      end if

      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",1,MA_ERR)

      return
      end




*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_dn_cmp_zv       *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transform of dn_cmp.
 
      subroutine nwpw_compcharge_gen_dn_cmp_zv(ispin,zv,dng_cmp)
      implicit none
      integer    ispin
      real*8     zv(*)
      complex*16 dng_cmp(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_cmp,1)

      if (isgamma) then

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            if (l.eq.0) then
               q00 = -zv(ia)/dsqrt(fourpi)
            else
               q00 = 0.0d0
            end if
            do m=-l,l
               
               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2) + q00
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_tc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do

      else

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
            if (l.eq.0) then
               q00 = -zv(ia)/dsqrt(fourpi)
            else
               q00 = 0.0d0
            end if
            do m=-l,l

               indx1 = (iii-1)*2*lm_size_max+ lm
               indx2 = (iii-1)*2*lm_size_max+(ispin-1)*lm_size_max + lm
               q = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2) + q00
               if (mod(l,2).eq.0) then
                  call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Cram_rc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)
      end do

      end if

      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",1,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_zv_cmp          *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transform of dn_cmp.
 
      subroutine nwpw_compcharge_gen_zv_cmp(zv,dng_cmp)
      implicit none
      real*8     zv(*)
      complex*16 dng_cmp(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_zv_cmp:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_cmp,1)

      if (isgamma) then

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         l = 0
         m = 0
         q00 = -zv(ia)/dsqrt(fourpi)

         q = q00
         call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))

         call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                          dcpl_mb(QTlm(1)))

         call Pack_tc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

      end do

      else

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         l = 0
         m = 0
         q00 = -zv(ia)/dsqrt(fourpi)

         q = q00
         call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))

         call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                          dcpl_mb(QTlm(1)))

         call Cram_rc_Mul2(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                       dcpl_mb(QTlm(1)))
         call Cram_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)
      end do

      end if

      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_zv_cmp:stack",1,MA_ERR)

      return
      end





*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_glm             *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transforms of dn_cmp

      subroutine nwpw_compcharge_gen_glm(ii,l,m,glm_out)
      implicit none
      integer ii,l,m
      complex*16 glm_out(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,lm,tmp(2),exi(2)

      iii = int_mb(ion_iontopaw(1)+ii-1)
      iia  = int_mb(katm_paw(1)+iii-1)
      if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
         lm = l*(l+1) + m

         value = MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
         value = value.and.
     >           MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
         if (.not.value)
     >     call errquit("nwpw_compcharge_gen_glm:stack",0,MA_ERR)

         if (isgamma) then
            call strfac_pack(0,ii,dcpl_mb(exi(1)))
            if (mod(l,2).eq.0) then
               call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                          dcpl_mb(exi(1)),
     >                          dcpl_mb(tmp(1)))
            else
               call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
            end if
            call Pack_tc_Mul(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                         dcpl_mb(tmp(1)),
     >                         glm_out)
         else
            call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
            if (mod(l,2).eq.0) then
               call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(tmp(1)))
            else
               call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                             dcpl_mb(exi(1)),
     >                             dcpl_mb(tmp(1)))
            end if
            call Cram_rc_Mul(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                         dcpl_mb(tmp(1)),
     >                         glm_out)

         end if
         value =           MA_pop_stack(exi(2))
         value = value.and.MA_pop_stack(tmp(2))
         if (.not.value)
     >     call errquit("nwpw_compcharge_gen_glm:stack",1,MA_ERR)


      else
         call errquit("nwpw_compcharge_gen_glm:bad l,m,",l,m)
      end if
      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_gen_glm2            *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transforms of dn_cmp and dn_cmp_smooth
 
      subroutine nwpw_compcharge_gen_glm2(ii,l,m,glm_out,glm_smooth_out)
      implicit none
      integer ii,l,m
      complex*16 glm_out(*)
      complex*16 glm_smooth_out(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,lm,tmp(2),exi(2)

      iii = int_mb(ion_iontopaw(1)+ii-1)
      iia  = int_mb(katm_paw(1)+iii-1)
      if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
         lm = l*(l+1) + m

         value = MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
         value = value.and.
     >           MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
         if (.not.value)
     >     call errquit("nwpw_compcharge_gen_glm2:stack",0,MA_ERR)

         if (isgamma) then
            call strfac_pack(0,ii,dcpl_mb(exi(1)))
            if (mod(l,2).eq.0) then
               call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                          dcpl_mb(exi(1)),
     >                          dcpl_mb(tmp(1)))
            else
               call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
            end if
            call Pack_tc_Mul(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                         dcpl_mb(tmp(1)),
     >                         glm_out)
            call Pack_tc_Mul(0,dbl_mb(gk_smooth(1)),
     >                         dcpl_mb(tmp(1)),
     >                         glm_smooth_out)
         else
            call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
            if (mod(l,2).eq.0) then
               call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(tmp(1)))
            else
               call Cram_rc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                             dcpl_mb(exi(1)),
     >                             dcpl_mb(tmp(1)))
            end if
            call Cram_rc_Mul(0,dbl_mb(gk(1)+(iia-1)*npack0),
     >                         dcpl_mb(tmp(1)),
     >                         glm_out)
            call Cram_rc_Mul(0,dbl_mb(gk_smooth(1)),
     >                         dcpl_mb(tmp(1)),
     >                         glm_smooth_out)
         end if

         value =           MA_pop_stack(exi(2))
         value = value.and.MA_pop_stack(tmp(2))
         if (.not.value)
     >     call errquit("nwpw_compcharge_gen_glm2:stack",1,MA_ERR)


      else
         call errquit("nwpw_compcharge_gen_glm2:bad l,m,",l,m)
      end if
      return
      end







*     *************************************************
*     *                                               *
*     *         nwpw_gintegrals_set_gcount            *
*     *                                               *
*     *************************************************

      subroutine nwpw_gintegrals_set_gcount()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer taskid,np,pcount,gcount,nshl3d
      integer l1,m1,l2,m2,iii,jjj,iia,jja

*     **** external functions ****
      integer  control_version,ewald_nshl3d
      external control_version,ewald_nshl3d

      call Parallel_taskid(taskid)
      call Parallel_np(np)

      periodic = (control_version().eq.3)

      if (periodic) then
         nshl3d = ewald_nshl3d()
      else
         nshl3d = 1
      end if

      pcount = 0
      gcount = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)

*        **** calculate on-site integrals ****
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
               gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end if
         end do
         end do

*        **** calculate IJ integrals ****
         do jjj=iii+1,nion_paw
            jja = int_mb(katm_paw(1)+jjj-1)

            do l1=0,int_mb(mult_l(1)+iia-1)
            do m1=-l1,l1
               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end do
            end do
         end do
      end do
      ngauss_max = gcount
      ngauss     = gcount
      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_gintegrals_init                *
*     *                                               *
*     *************************************************

      subroutine nwpw_gintegrals_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value


      call nwpw_gintegrals_set_gcount()

      value =           MA_alloc_get(mt_int,ngauss_max,"lm1_gauss",
     >                              lm1_gauss(2),lm1_gauss(1))
      value = value.and.MA_alloc_get(mt_int,ngauss_max,"lm2_gauss",
     >                              lm2_gauss(2),lm2_gauss(1))
      value = value.and.MA_alloc_get(mt_int,ngauss_max,"iii1_gauss",
     >                              iii1_gauss(2),iii1_gauss(1))
      value = value.and.MA_alloc_get(mt_int,ngauss_max,"iii2_gauss",
     >                              iii2_gauss(2),iii2_gauss(1))
      value = value.and.MA_alloc_get(mt_dbl,ngauss_max,"e_gauss",
     >                              e_gauss(2),e_gauss(1))
      value = value.and.MA_alloc_get(mt_dbl,3*ngauss_max,"f_gauss",
     >                              f_gauss(2),f_gauss(1))
      if (.not.value)
     > call errquit("nwpw_gintegrals_init:cannot allocate memory",
     >             0,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_gintegrals_end                 *
*     *                                               *
*     *************************************************
      subroutine nwpw_gintegrals_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value

      value =           MA_free_heap(lm1_gauss(2))
      value = value.and.MA_free_heap(lm2_gauss(2))
      value = value.and.MA_free_heap(iii1_gauss(2))
      value = value.and.MA_free_heap(iii2_gauss(2))
      value = value.and.MA_free_heap(e_gauss(2))
      value = value.and.MA_free_heap(f_gauss(2))
      if (.not.value)
     > call errquit("nwpw_gintegrals_end:cannot allocate memory",
     >             0,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *             nwpw_gintegrals_set               *
*     *                                               *
*     *************************************************
      subroutine nwpw_gintegrals_set(move)
      implicit none
      logical move

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      real*8 tole
      parameter (tole=1.0d-15)

      integer taskid,np,pcount
      integer ii,jj,ia,ja,indx
      integer iii,jjj,iia,jja
      integer lm1,l1,m1,lm2,l2,m2
      integer l,nshl3d,rcell,rcell_hndl
      real*8 R1(3),R12(3),s1,s2,Rab(3),R
      real*8 W1,W2,W3,W4,dW1(3),dW2(3),dW3(3),dW4(3)
      real*8 e1,de1(3)

*     **** external functions ****
      real*8   ion_rion,nwpw_WGaussian,nwpw_UGaussian
      external ion_rion,nwpw_WGaussian,nwpw_UGaussian
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      integer  ewald_nshl3d,ewald_rcell_ptr
      external ewald_nshl3d,ewald_rcell_ptr


      call Parallel_taskid(taskid)
      call Parallel_np(np)


      if (periodic) then
         nshl3d = ewald_nshl3d()
         rcell  = ewald_rcell_ptr()
      else
         if (.not. MA_push_get(mt_dbl,3,"rcellflm",rcell_hndl,rcell))
     >   call errquit("nwpw_compcharge_set_gintegrals:stack",1,MA_ERR)

         nshl3d = 1
         dbl_mb(rcell)   = 0.0d0
         dbl_mb(rcell+1) = 0.0d0
         dbl_mb(rcell+2) = 0.0d0
      end if
      call dcopy(ngauss_max,0.0d0,0,dbl_mb(e_gauss(1)),1)
      call dcopy(3*ngauss_max,0.0d0,0,dbl_mb(f_gauss(1)),1)

      pcount = 0
      indx   = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)
         s1  = dbl_mb(sigma_paw(1)+iia-1)

*        **** calculate on-site integrals ****
         lm1 = 0
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
                W1=nwpw_UGaussian(l1,m1,s1,l1,m1,s1)
                W2=nwpw_UGaussian(l1,m1,s1,l1,m1,sigma_smooth)
                W4=nwpw_UGaussian(l1,m1,sigma_smooth,l1,m1,sigma_smooth)
                e1 = 0.5d0*W1 + 0.5d0*W4 - W2
                if (dabs(e1).gt.tole) then
                  dbl_mb(e_gauss(1)+indx) = e1
                  int_mb(lm1_gauss(1)+indx) = (iii-1)*2*lm_size_max+ lm1
                  int_mb(lm2_gauss(1)+indx) = (iii-1)*2*lm_size_max+ lm1
                  int_mb(iii1_gauss(1)+indx) = iii
                  int_mb(iii2_gauss(1)+indx) = iii

                  indx = indx + 1
                end if
               
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               lm2 = 0
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     e1 = 0.0d0
                     do l=2,nshl3d
                        Rab(1) = dbl_mb(rcell+l-1)
                        Rab(2) = dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                           W1 = nwpw_WGaussian(l1,m1,s1,l2,m2,s1,Rab)
                           W2 = nwpw_WGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
c                           W3 = nwpw_WGaussian(l1,m1,sigma_smooth,
c     >                                         l2,m2,s1,Rab)
                           W4 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           e1=e1+0.5d0*W1+0.5d0*W4-W2
                        end if
                     end do
                     if (dabs(e1).gt.tole) then
                        dbl_mb(e_gauss(1)+indx) = e1
                        int_mb(lm1_gauss(1)+indx) 
     >                     = (iii-1)*2*lm_size_max+lm1
                        int_mb(lm2_gauss(1)+indx)
     >                     = (iii-1)*2*lm_size_max+lm2
                        int_mb(iii1_gauss(1)+indx) = iii
                        int_mb(iii2_gauss(1)+indx) = iii

                        indx = indx + 1
                     end if
                  end if
                  pcount = pcount + 1
                  lm2 = lm2 + 1
               end do
               end do
            end if

            lm1 = lm1 + 1
         end do
         end do


*        **** calculate IJ integrals ****
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         R1(1) = ion_rion(1,ii)
         R1(2) = ion_rion(2,ii)
         R1(3) = ion_rion(3,ii)
         do jjj=iii+1,nion_paw
            jja = int_mb(katm_paw(1)+jjj-1)
            s2  = dbl_mb(sigma_paw(1)+jja-1)

            jj  = int_mb(ion_pawtoion(1)+jjj-1)
            R12(1) = R1(1) - ion_rion(1,jj)
            R12(2) = R1(2) - ion_rion(2,jj)
            R12(3) = R1(3) - ion_rion(3,jj)

            lm1 = 0
            do l1=0,int_mb(mult_l(1)+iia-1)
            do m1=-l1,l1

               lm2 = 0
               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     e1 = 0.0d0
                     de1(1) = 0.0d0
                     de1(2) = 0.0d0
                     de1(3) = 0.0d0
                     do l=1,nshl3d
                        Rab(1) = R12(1) + dbl_mb(rcell+l-1)
                        Rab(2) = R12(2) + dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = R12(3) + dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                        if (move) then
                           call nwpw_dWGaussian(l1,m1,s1,l2,m2,s2,Rab,
     >                                          W1,dW1)
                           call nwpw_dWGaussian(l1,m1,s1,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W2,dW2)
                           call nwpw_dWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,s2,Rab,
     >                                          W3,dW3)
                           call nwpw_dWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W4,dW4)
                           e1 = e1 + W1 + W4 - W2 - W3
                           de1(1)=de1(1)+dW1(1)+dW4(1)-dW2(1)-dW3(1)
                           de1(2)=de1(2)+dW1(2)+dW4(2)-dW2(2)-dW3(2)
                           de1(3)=de1(3)+dW1(3)+dW4(3)-dW2(3)-dW3(3)

                        else
                           W1 = nwpw_WGaussian(l1,m1,s1,l2,m2,s2,Rab)
                           W2 = nwpw_WGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
                           W3 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                      l2,m2,s2,Rab)
                           W4 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           e1 = e1 + W1 + W4 - W2 - W3
                        end if
                        end if
                     end do
                     if (dabs(e1).gt.tole) then
                        dbl_mb(e_gauss(1)+indx) = e1
                        if (move) then
                           dbl_mb(f_gauss(1)+3*indx)   = de1(1)
                           dbl_mb(f_gauss(1)+3*indx+1) = de1(2)
                           dbl_mb(f_gauss(1)+3*indx+2) = de1(3)
                        end if
                        int_mb(lm1_gauss(1)+indx)
     >                     = (iii-1)*2*lm_size_max+lm1
                        int_mb(lm2_gauss(1)+indx)
     >                     = (jjj-1)*2*lm_size_max+lm2
                        int_mb(iii1_gauss(1)+indx) = iii
                        int_mb(iii2_gauss(1)+indx) = jjj

                        indx = indx + 1
                     end if

                  end if
                  pcount = pcount + 1
                  lm2 = lm2 + 1
               end do
               end do

               lm1 = lm1 + 1
            end do
            end do

         end do
      end do
      ngauss = indx

      if (.not.periodic) then
        if (.not.MA_pop_stack(rcell_hndl))
     >   call errquit("nwpw_compcharge_set_gintegrals:stack",2,MA_ERR)
      end if

      return
      end

*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_E_multipole         *
*     *                                               *
*     *************************************************
*
*    This routine calulates the multipole energy defined by
* 
*                         //
*    Emultipole  = 0.5 * || (ncmp(r)-ncmp_smooth(r))*(ncmp(r')-ncmp_smooth(r'))
*                        || ---------------------------------------------------  dr dr'
*                        //                       |r-r'|
*
*    using two electron two center Gaussian Coulomb integrals.
*
*     Uses - compensation charge Qlm, and a table of Gaussian integrals stored in e_gauss
*
      real*8 function nwpw_compcharge_E_multipole(ispin)
      implicit none
      integer ispin

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 E,q1,q2
      integer iii,jjj,l1,m1,l2,m2


      E = 0.0d0
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max

         call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
         call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)



         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)
         E = E + q1*q2*dbl_mb(e_gauss(1)+k-1)

      end do
      call Parallel_SumAll(E)

      nwpw_compcharge_E_multipole = E
      return
      end

*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_F_multipole         *
*     *                                               *
*     *************************************************
*
      subroutinen nwpw_compcharge_F_multipole(ispin,fion)
      implicit none
      integer ispin
      real*8 fion(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      integer iii1,iii2,ii1,ii2
      integer fmult(2)
      real*8 q1,q1x,q1y,q1z,q2,q2x,q2y,q2z
      real*8 tmp,tmpx,tmpy,tmpz,R,R12(3)

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      if (.not.MA_push_get(mt_dbl,3*nion_paw,"fmult",fmult(2),fmult(1)))
     >   call errquit("nwpw_compcharge_F_multipole:stack",1,MA_ERR)
      call dcopy(3*nion_paw,0.0d0,0,dbl_mb(fmult(1)))

      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max
         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)

         q1x = dbl_mb(Qlmx(1)+indx1up) + dbl_mb(Qlmx(1)+indx1dn)
         q1y = dbl_mb(Qlmy(1)+indx1up) + dbl_mb(Qlmy(1)+indx1dn)
         q1z = dbl_mb(Qlmz(1)+indx1up) + dbl_mb(Qlmz(1)+indx1dn)

         q2x = dbl_mb(Qlmx(1)+indx2up) + dbl_mb(Qlmx(1)+indx2dn)
         q2y = dbl_mb(Qlmy(1)+indx2up) + dbl_mb(Qlmy(1)+indx2dn)
         q2z = dbl_mb(Qlmz(1)+indx2up) + dbl_mb(Qlmz(1)+indx2dn)

         iii1 = int_mb(iii1_gauss(1)+k-1)
         iii2 = int_mb(iii2_gauss(1)+k-1)
         ii1 = int_mb(ion_pawtoion(1)+iii1-1)
         ii2 = int_mb(ion_pawtoion(1)+iii2-1)
         R12(1) = ion_rion(1,ii1) - ion_rion(1,ii2)
         R12(2) = ion_rion(2,ii1) - ion_rion(2,ii2)
         R12(3) = ion_rion(3,ii1) - ion_rion(3,ii2)
         R = dsqrt(R12(1)*R12(1)+R12(2)*R12(2)+R12(3)*R12(3))
         R12(1) = R12(1)/R
         R12(2) = R12(2)/R
         R12(3) = R12(3)/R
        
         tmp  =       dbl_mb(e_gauss(1)+3*(k-1))
         tmpx = q1*q2*dbl_mb(f_gauss(1)+3*(k-1))  *R12(1)
         tmpy = q1*q2*dbl_mb(f_gauss(1)+3*(k-1)+1)*R12(2)
         tmpz = q1*q2*dbl_mb(f_gauss(1)+3*(k-1)+2)*R12(3)

         dbl_mb(fmult(1)+3*(iii1-1))   = dbl_mb(fmult(1)+3*(iii1-1))
     >                                 - tmpx - q1x*q2*tmp
         dbl_mb(fmult(1)+3*(iii1-1)+1) = dbl_mb(fmult(1)+3*(iii1-1)+1)
     >                                 - tmpy - q1y*q2*tmp
         dbl_mb(fmult(1)+3*(iii1-1)+2) = dbl_mb(fmult(1)+3*(iii1-1)+2) 
     >                                 - tmpz - q1z*q2*tmp

         dbl_mb(fmult(1)+3*(iii2-1))   = dbl_mb(fmult(1)+3*(iii2-1))
     >                                 + tmpx - q1*q2x*tmp
         dbl_mb(fmult(1)+3*(iii2-1)+1) = dbl_mb(fmult(1)+3*(iii2-1)+1)
     >                                 + tmpy - q1*q2y*tmp
         dbl_mb(fmult(1)+3*(iii2-1)+2) = dbl_mb(fmult(1)+3*(iii2-1)+2) 
     >                                 + tmpz - q1*q2z*tmp
      end do
      call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(fmult(1)))


      if (.not.MA_pop_stack(fmult(2)))
     >  call errquit("nwpw_compcharge_F_multipole:stack",2,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *        nwpw_compcharge_add_dEmult_Qlm         *
*     *                                               *
*     *************************************************
*
*    This routine calculates the derivative of the multipole
*    energy wrt to Qlm
*
      subroutine nwpw_compcharge_add_dEmult_Qlm(ispin)
      implicit none
      integer ispin

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 q1,q2

c      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dEmult_Qlm(1)),1)
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max
         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)

         dbl_mb(dE_Qlm(1)+indx1up) = dbl_mb(dE_Qlm(1)+indx1up) 
     >                                 + q2*dbl_mb(e_gauss(1)+k-1)
         dbl_mb(dE_Qlm(1)+indx1dn) = dbl_mb(dE_Qlm(1)+indx1dn) 
     >                                 + q2*dbl_mb(e_gauss(1)+k-1)

         dbl_mb(dE_Qlm(1)+indx2up) = dbl_mb(dE_Qlm(1)+indx2up) 
     >                                 + q1*dbl_mb(e_gauss(1)+k-1)
         dbl_mb(dE_Qlm(1)+indx2dn) = dbl_mb(dE_Qlm(1)+indx2dn) 
     >                                 + q1*dbl_mb(e_gauss(1)+k-1)

      end do
c      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
c     >                            dbl_mb(dE_Qlm(1)))
      
      return
      end


      real*8 function nwpw_compcharge_E_multipole_zv(ispin,zv)
      implicit none
      integer ispin
      real*8 zv(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 E,q1,q2,q11,q22,fourpi
      integer iii,jjj,l1,m1,l2,m2,ia,iia,ja,jja

      fourpi = 16.0d0*datan(1.0d0)

      E = 0.0d0
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max

         call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
         call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)

         if (l1.eq.0) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            q11 = -zv(ia)/dsqrt(fourpi)
         else
            q11 = 0.0d0
         end if

         if (l2.eq.0) then
            jja = int_mb(katm_paw(1)+jjj-1)
            ja  = int_mb(katm_pawtoion(1)+jja-1)
            q22 = -zv(ia)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if


         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)+q11
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)+q22
         E = E + q1*q2*dbl_mb(e_gauss(1)+k-1)

c         if (dabs(q1*q2*dbl_mb(e_gauss(1)+k-1)).gt.1.0d-15) then
c         write(*,'(A,6I4,E20.12,2E15.6)') "ii,li,mi,jj,lj,mj:",
c     >             iii,l1,m1,jjj,l2,m2,
c     >             dbl_mb(e_gauss(1)+k-1),q1,q2
c         end if 
      end do
      call Parallel_SumAll(E)

      nwpw_compcharge_E_multipole_zv = E
      return
      end





*     *************************************************
*     *                                               *
*     *         nwpw_compcharge_coulomb_atom          *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_coulomb_atom(ii,ia,ispin,ne,nprj,
     >                                        wmatrix,sw1,sw2,eh_atom)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  wmatrix(nprj,nprj,ispin)
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      real*8  eh_atom

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iia,shift

      call nwpw_timing_start(4)
      call nwpw_timing_start(22)

      iia = int_mb(katm_iontopaw(1)+ia-1)

      shift = int_mb(shift_hartree(1)+iia-1)
      call nwpw_compcharge_coulomb_sub(ispin,ne,nprj,wmatrix,sw1,sw2,
     >                           eh_atom,
     >                           int_mb(nindx_hartree(1)+iia-1),
     >                           int_mb(iprj_hartree(1)+shift),
     >                           int_mb(jprj_hartree(1)+shift),
     >                           int_mb(iprj1_hartree(1)+shift),
     >                           int_mb(jprj1_hartree(1)+shift),
     >                           dbl_mb(coeff_hartree(1)+shift))


      call nwpw_timing_end(4)
      call nwpw_timing_end(22)
      return
      end

c     *********************************************
c     *                                           *
c     *           nwpw_compcharge_coulomb_sub     *
c     *                                           *
c     *********************************************
      subroutine nwpw_compcharge_coulomb_sub(ispin,ne,nprj,wmatrix,
     >                                 sw1,sw2,eh,
     >                                 nindx,
     >                                 iprj_hartree, jprj_hartree,
     >                                 iprj1_hartree,jprj1_hartree,
     >                                 coeff_hartree)
      implicit none
      integer ispin,ne(2),nprj
      real*8  wmatrix(nprj,nprj,ispin)
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      real*8  eh
      integer n1dgrid,nbasis,lmax2
      integer nindx,iprj_hartree(*),jprj_hartree(*)
      integer iprj1_hartree(*),jprj1_hartree(*)
      real*8  coeff_hartree(*)

      integer n,i,iprj,jprj,iprj1,jprj1
      real*8  coeff,w,w1,scal

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      scal = 1.0d0/lattice_omega()

      do i=1,nindx
         iprj  = iprj_hartree(i)
         jprj  = jprj_hartree(i)
         iprj1 = iprj1_hartree(i)
         jprj1 = jprj1_hartree(i)
         coeff = coeff_hartree(i)
         w  = wmatrix(iprj,jprj,1)   + wmatrix(iprj,jprj,ispin)
         w1 = wmatrix(iprj1,jprj1,1) + wmatrix(iprj1,jprj1,ispin)
         eh = eh + coeff*w*w1*scal*scal*0.5d0
         do n=1,ne(1)+ne(2)
            sw2(n,iprj) = sw2(n,iprj) 
     >                  + scal*coeff*sw1(n,jprj)*w1
         end do
      end do
   
      call nwpw_timing_end(21)
      return
      end


*     *************************************************
*     *                                               *
*     *         nwpw_compcharge_coulomb_e_atom        *
*     *                                               *
*     *************************************************
      real*8 function nwpw_compcharge_coulomb_e_atom(ii,ia,ispin,nprj,
     >                                               wmatrix)
      implicit none
      integer ii,ia
      integer ispin,nprj
      real*8  wmatrix(nprj,nprj,ispin)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iia,shift
      real*8  eh_atom

*     **** external functions ****
      real*8   nwpw_compcharge_coulomb_e_sub
      external nwpw_compcharge_coulomb_e_sub

      call nwpw_timing_start(4)
      call nwpw_timing_start(22)

      iia = int_mb(katm_iontopaw(1)+ia-1)

      shift = int_mb(shift_hartree(1)+iia-1)
      eh_atom = nwpw_compcharge_coulomb_e_sub(ispin,nprj,wmatrix,
     >                           int_mb(nindx_hartree(1)+iia-1),
     >                           int_mb(iprj_hartree(1)+shift),
     >                           int_mb(jprj_hartree(1)+shift),
     >                           int_mb(iprj1_hartree(1)+shift),
     >                           int_mb(jprj1_hartree(1)+shift),
     >                           dbl_mb(coeff_hartree(1)+shift))
      call nwpw_timing_end(4)
      call nwpw_timing_end(22)

      nwpw_compcharge_coulomb_e_atom = eh_atom
      return
      end



c     *********************************************
c     *                                           *
c     *         nwpw_compcharge_coulomb_e_sub     *
c     *                                           *
c     *********************************************
      real*8 function nwpw_compcharge_coulomb_e_sub(ispin,nprj,wmatrix,
     >                                 nindx,
     >                                 iprj_hartree, jprj_hartree,
     >                                 iprj1_hartree,jprj1_hartree,
     >                                 coeff_hartree)
      implicit none
      integer ispin,nprj
      real*8  wmatrix(nprj,nprj,ispin)
      integer n1dgrid,nbasis,lmax2
      integer nindx,iprj_hartree(*),jprj_hartree(*)
      integer iprj1_hartree(*),jprj1_hartree(*)
      real*8  coeff_hartree(*)

      integer n,i,iprj,jprj,iprj1,jprj1
      real*8  coeff,w,w1,scal
      real*8  eh

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      eh = 0.0d0
      scal = 1.0d0/lattice_omega()
      do i=1,nindx
         iprj  = iprj_hartree(i)
         jprj  = jprj_hartree(i)
         iprj1 = iprj1_hartree(i)
         jprj1 = jprj1_hartree(i)
         coeff = coeff_hartree(i)
         w  = wmatrix(iprj,jprj,1)   + wmatrix(iprj,jprj,ispin)
         w1 = wmatrix(iprj1,jprj1,1) + wmatrix(iprj1,jprj1,ispin)
         eh = eh + coeff*w*w1*scal*scal*0.5d0
      end do

      nwpw_compcharge_coulomb_e_sub = eh
      return
      end



      subroutine nwpw_compcharge_indxiiilm(indx,iii,l,m)
      implicit none
      integer indx
      integer iii,l,m

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      integer lm

      iii = indx/(2*lm_size_max) + 1
      lm = indx - (iii-1)*2*lm_size_max
      l = 0
      do while ((lm+1) .gt. (l+1)**2)
         l = l + 1
      end do
      m = lm - l*(l+1)
      return
      end





