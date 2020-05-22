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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      logical value
      integer k,l,m,lm,ii,iii,ia,iia,Gindx,Gall(2),G(3),tsize
      integer iasize,i,j,li,mi,lj,mj,bi,bj,nb
      integer i1,j1,li1,mi1,bi1,lj1,mj1,bj1
      integer indx,comp_charge_tag,comp_indx,matr_ptr
      real*8  scal,gg,fourpioveromega,sumall
      real*8  taunt,rs,pi,w,rcut,phase

*     **** external functions ****
      integer  control_code,G_indx,c_G_indx,ion_katm
      external control_code,G_indx,c_G_indx,ion_katm
      integer  nwpw_doublefactorial,control_version
      external nwpw_doublefactorial,control_version
      real*8   lattice_omega
      external lattice_omega
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
      if ((control_version().eq.4).and.(rcut.le.0.0d0)) rcut = 1.0d0
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
      value =           BA_alloc_get(mt_int,nion,"ion_iontopaw",
     >                               ion_iontopaw(2),ion_iontopaw(1))
      value = value.and.BA_alloc_get(mt_int,nion_paw,"ion_pawtoion",
     >                               ion_pawtoion(2),ion_pawtoion(1))
      value = value.and.BA_alloc_get(mt_int,nkatm,"katm_iontopaw",
     >                               katm_iontopaw(2),katm_iontopaw(1))
      value = value.and.BA_alloc_get(mt_int,nkatm_paw,"katm_pawtoion",
     >                               katm_pawtoion(2),katm_pawtoion(1))
      value = value.and.BA_alloc_get(mt_int,nion_paw,"katm_paw",
     >                               katm_paw(2),katm_paw(1))
      value = value.and.BA_alloc_get(mt_dbl,nkatm_paw,"sigma_paw",
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
      end do


*     **** allocate mult_l and lm_size ****
      value = value.and.BA_alloc_get(mt_int,nkatm_paw,
     >                     "mult_l",mult_l(2),mult_l(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nkatm_paw,
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
      call nwpw_gaunt_init(.false.,2*mult_l_max)

      !*** allocate gk_smooth, gk,and glm ***
      lm_size_max = (mult_l_max+1)**2
      value = BA_alloc_get(mt_dbl,npack0,
     >                     "vk_smooth",vk_smooth(2),vk_smooth(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,npack0,
     >                     "gk_smooth",gk_smooth(2),gk_smooth(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,npack0*nkatm_paw,"gk",gk(2),gk(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,npack0*lm_size_max,
     >                     "glm",glm(2),glm(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlm",Qlm(2),Qlm(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmx",Qlmx(2),Qlmx(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmy",Qlmy(2),Qlmy(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "Qlmz",Qlmz(2),Qlmz(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "dEmult_Qlm",dEmult_Qlm(2),dEmult_Qlm(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "dElocal_Qlm",dElocal_Qlm(2),dElocal_Qlm(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nion_paw*lm_size_max,
     >                     "dE_Qlm",dE_Qlm(2),dE_Qlm(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nfft3d,"Gall",Gall(2),Gall(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmx(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmy(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(Qlmz(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dEmult_Qlm(1)),1)
      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,
     >           dbl_mb(dElocal_Qlm(1)),1)
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

         phase = 1.0d0 
         if (mod(l,4).eq.0) then
            phase = 1.0d0
         else if (mod(l,4).eq.1) then
            phase = -1.0d0
         else if (mod(l,4).eq.2) then
            phase = -1.0d0
         else if (mod(l,4).eq.3) then
            phase = 1.0d0
         end if

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

*        **** define glm = (-i)**l *  |k|**l * Tlm(k)/ (2*l+1)!! ****
         do m =-l,l
            call Tesseral3_vector_lm(l,m,npack0,
     >                               dbl_mb(G(1)),
     >                               dbl_mb(G(2)),
     >                               dbl_mb(G(3)),
     >                               dbl_mb(glm(1)+lm*npack0))
            do k=1,npack0
               dbl_mb(glm(1)+lm*npack0+k-1) 
     >          = phase*dbl_mb(glm(1)+lm*npack0+k-1)*dbl_mb(gk(1)+k-1)
            end do

            lm = lm + 1
         end do

      end do

*     **** define vk_smooth(k) ****
      call nwpw_compcharge_gen_vk_smooth(isgamma,sigma_smooth,npack0,
     >                                   dbl_mb(G(1)),
     >                                   dbl_mb(G(2)),
     *                                   dbl_mb(G(3)),
     >                                   dbl_mb(vk_smooth(1)))

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
      if (.not.BA_pop_stack(Gall(2)))
     >   call errquit("nwpw_compcharge_init: pop stack",0,MA_ERR)



*     *******************************************************
*     *****  define indexing for compcharge evalulations ****
*     *******************************************************

      if (.not.BA_alloc_get(mt_int,nkatm_paw,"nindx_Tndiff",
     >                      nindx_Tndiff(2),nindx_Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      if (.not.BA_alloc_get(mt_int,nkatm_paw,"shift_Tndiff",
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
                     if (dabs(taunt).gt.1.0d-15) then
                        indx = indx + 1
                     end if
                  end do
               end do
            end do
         end do
         
        int_mb(nindx_Tndiff(1)+iia-1)=indx-int_mb(shift_Tndiff(1)+iia-1)
      end do

      value = BA_alloc_get(mt_int,indx,"lm_Tndiff",
     >                     lm_Tndiff(2),lm_Tndiff(1))
      value = value.and.
     >        BA_alloc_get(mt_int,indx,"iprj_Tndiff",
     >                     iprj_Tndiff(2),iprj_Tndiff(1))
      value = value.and.
     >        BA_alloc_get(mt_int,indx,"jprj_Tndiff",
     >                     jprj_Tndiff(2),jprj_Tndiff(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,indx,"coeff_Tndiff",
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

                     if (dabs(taunt).gt.1.0d-15) then
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

      if (.not.BA_alloc_get(mt_int,nkatm_paw,"nindx_hartree",
     >                      nindx_hartree(2),nindx_hartree(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",5,MA_ERR)

      if (.not.BA_alloc_get(mt_int,nkatm_paw,"shift_hartree",
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
                           if (dabs(taunt).gt.1.0d-15) then
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


      value = BA_alloc_get(mt_int,indx,"iprj_hartree",
     >                 iprj_hartree(2),iprj_hartree(1))
      value = value.and.
     >     BA_alloc_get(mt_int,indx,"jprj_hartree",
     >                 jprj_hartree(2),jprj_hartree(1))
      value = value.and.
     >     BA_alloc_get(mt_int,indx,"iprj1_hartree",
     >                 iprj1_hartree(2),iprj1_hartree(1))
      value = value.and.
     >     BA_alloc_get(mt_int,indx,"jprj1_hartree",
     >                 jprj1_hartree(2),jprj1_hartree(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,indx,"coeff_hartree",
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
                           if (dabs(taunt).gt.1.0d-15) then
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      logical value

      call nwpw_gaunt_end()

      value =           BA_free_heap(ion_iontopaw(2))
      value = value.and.BA_free_heap(ion_pawtoion(2))
      value = value.and.BA_free_heap(katm_iontopaw(2))
      value = value.and.BA_free_heap(katm_pawtoion(2))
      value = value.and.BA_free_heap(katm_paw(2))
      value = value.and.BA_free_heap(sigma_paw(2))
      value = value.and.BA_free_heap(mult_l(2))
      value = value.and.BA_free_heap(lm_size(2))
      value = value.and.BA_free_heap(vk_smooth(2))
      value = value.and.BA_free_heap(gk_smooth(2))
      value = value.and.BA_free_heap(gk(2))
      value = value.and.BA_free_heap(glm(2))
      value = value.and.BA_free_heap(Qlm(2))
      value = value.and.BA_free_heap(Qlmx(2))
      value = value.and.BA_free_heap(Qlmy(2))
      value = value.and.BA_free_heap(Qlmz(2))
      value = value.and.BA_free_heap(dEmult_Qlm(2))
      value = value.and.BA_free_heap(dElocal_Qlm(2))
      value = value.and.BA_free_heap(dE_Qlm(2))
      value = value.and.BA_free_heap(shift_Tndiff(2))
      value = value.and.BA_free_heap(nindx_Tndiff(2))
      value = value.and.BA_free_heap(lm_Tndiff(2))
      value = value.and.BA_free_heap(iprj_Tndiff(2))
      value = value.and.BA_free_heap(jprj_Tndiff(2))
      value = value.and.BA_free_heap(coeff_Tndiff(2))

      value = value.and.BA_free_heap(shift_hartree(2))
      value = value.and.BA_free_heap(nindx_hartree(2))
      value = value.and.BA_free_heap(iprj_hartree(2))
      value = value.and.BA_free_heap(jprj_hartree(2))
      value = value.and.BA_free_heap(iprj1_hartree(2))
      value = value.and.BA_free_heap(jprj1_hartree(2))
      value = value.and.BA_free_heap(coeff_hartree(2))

      if (.not.value)
     >   call errquit("nwpw_compcharge_end: freeing heap",0,MA_ERR)
      return
      end 


*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_gen_vk_smooth        *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_vk_smooth(isgamma,sigma_smooth,
     >                                         npack0,Gx,Gy,Gz,vk)
      implicit none
      logical isgamma
      real*8  sigma_smooth
      integer npack0
      real*8  Gx(*),Gy(*),Gz(*),vk(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,nrho,nray,nx,ny,nz
      real*8  unita(3,3),bmesh,log_bmesh
      integer Gray(2),vlray(2),tmpray(2),f(2),rho(2)

*     **** external functions ****
      integer  nwpw_kbpp_calc_nray
      external nwpw_kbpp_calc_nray
      real*8   lattice_unita
      external lattice_unita

*     **** define extra local smooth psp ****
      if (isgamma) then
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
      else
         call C3dB_nx(1,nx)
         call C3dB_ny(1,ny)
         call C3dB_nz(1,nz)
      end if

      do j=1,3
         do i=1,3
            unita(i,j) = lattice_unita(i,j)
         end do
      end do

      nray = nwpw_kbpp_calc_nray(nx,ny,nz,unita)
      bmesh     = 1.005d0
      log_bmesh = dlog(bmesh)
      nrho = int(dlog(25.0d0/0.00025d0)/log_bmesh) + 1

      !*** make sure loggrid is odd ***
      if (mod(nrho,2).eq.0) nrho = nrho + 1

      value = BA_push_get(mt_dbl,nray,'Gray',Gray(2),Gray(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nray,'vlray',vlray(2),vlray(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nray,'tmpray',tmpray(2),tmpray(1))

      value = value.and.
     >        BA_push_get(mt_dbl,nrho,'fjky99',f(2),f(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nrho,'rhoy99',rho(2),rho(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_vk_smooth:stack",4,MA_ERR)

      call nwpw_kbpp_generate_G_ray(nx,ny,nz,unita,dbl_mb(Gray(1)))

      dbl_mb(rho(1)) = 0.00025d0
      do i=2,nrho
         dbl_mb(rho(1)+i-1) = bmesh*dbl_mb(rho(1)+i-2)
      end do


      call nwpw_compcharge_gen_vlray(sigma_smooth,
     >                             log_bmesh,nrho,
     >                             dbl_mb(rho(1)),dbl_mb(f(1)),
     >                             nray,
     >                             dbl_mb(Gray(1)),
     >                             dbl_mb(vlray(1)),
     >                             dbl_mb(tmpray(1)))

      call nwpw_compcharge_gen_smoothpsp(nray,
     >                                   dbl_mb(Gray(1)),
     >                                   dbl_mb(vlray(1)),
     >                                   npack0,Gx,Gy,Gz,
     >                                   vk)

      value =           BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(f(2))
      value = value.and.BA_pop_stack(tmpray(2))
      value = value.and.BA_pop_stack(vlray(2))
      value = value.and.BA_pop_stack(Gray(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_vk_smooth:pop stack",4,MA_ERR)

      return
      end 

*     ***********************************
*     *                                 *
*     *   nwpw_compcharge_gen_vlray     *
*     *                                 *
*     ***********************************

      subroutine nwpw_compcharge_gen_vlray(sigma_smooth,
     >                                  log_amesh,nrho,rho,f,
     >                                  nray,Gray,vlray,tmpray)
      implicit none
      real*8  sigma_smooth,log_amesh
      integer nrho
      real*8  rho(nrho),f(nrho)
      integer nray
      real*8  Gray(nray)
      real*8  vlray(nray,2)
      real*8  tmpray(nray)

*     **** local variables ****
      logical periodic
      integer i,k
      real*8 q,sn,cs,xerf1,yerf1,xerf,yerf,yp1,dG,rlocal
      real*8 twopi,fourpi,ecut

*     **** external functions ****
      logical  control_kbpp_filter
      external control_kbpp_filter
      integer  control_version
      external control_version
      real*8   util_erf,log_integrate_def,control_ecut
      external util_erf,log_integrate_def,control_ecut

      periodic = (control_version().eq.3)

      twopi  = 8.0d0*datan(1.0d0)
      fourpi = 2.0d0*twopi
      ecut   = control_ecut()
      rlocal = 1.0d0

      call dcopy(nray,0.0d0,0,vlray,1)
      do k=2,nray
         q = Gray(k)
c         do i=1,nrho
c            sn = dsin(q*rho(i))
c            xerf1 = rho(i)/sigma_smooth
c            yerf1 = util_erf(xerf1)
c            f(i)=(-yerf1)*sn
c         end do
c         cs = dcos(q*rho(nrho))

c   aperiodic should not use this?
c         if (periodic) then
c            vlray(k,1)=fourpi
c     >            *log_integrate_def(0,f,0,rho,log_amesh,nrho)/q
c     >           -fourpi*cs/(q*q)
            vlray(k,1)=-(fourpi/(q*q))*dexp(-0.25d0*(sigma_smooth*q)**2)
c         else
c            vlray(k,1)=fourpi
c     >            *log_integrate_def(0,f,0,rho,log_amesh,nrho)/q
c         end if
      end do

c      do i=1,nrho
c         xerf1 = rho(i)/sigma_smooth
c         yerf1 = util_erf(xerf1)
c         !f(i)=(-yerf1)*rho(i)
c         f(i)=(-yerf1)
c      end do

c   aperiodic shou not use this?
c      if (periodic) then
c         vlray(1,1) = fourpi*log_integrate_def(1,f,1,rho,log_amesh,nrho)
c     >              + twopi*rho(nrho)**2
         vlray(1,1) = 0.5d0*twopi*sigma_smooth**2
c         write(*,*) "Vksmooth(G=0)=",vlray(1,1),
c     >              0.5d0*twopi*sigma_smooth**2
c      else
c         vlray(1,1) = fourpi*log_integrate_def(0,f,0,rho,log_amesh,nrho)
c      end if

      if (control_kbpp_filter())
     >  call nwpw_kbpp_filter_ray(nray,Gray,ecut,vlray)

      dG = Gray(3)-Gray(2)

      !**** five point formula ***
      yp1 = ( -50.0d0*vlray(2,1)
     >       + 96.0d0*vlray(3,1)
     >       - 72.0d0*vlray(4,1)
     >       + 32.0d0*vlray(5,1)
     >       -  6.0d0*vlray(6,1))/(24.0d0*dG)
      call nwpw_spline(Gray(2),vlray(2,1),nray-1,yp1,0.0d0,
     >                 vlray(2,2),tmpray)

      return
      end

*     ***********************************
*     *                                 *
*     *  nwpw_compcharge_gen_smoothpsp  *
*     *                                 *
*     ***********************************
      subroutine nwpw_compcharge_gen_smoothpsp(nray,Gray,vlray,
     >                                      npack0,Gx,Gy,Gz,
     >                                      vlsmooth)
      implicit none
      integer nray
      real*8  Gray(nray),vlray(nray,2)
      integer npack0
      real*8  Gx(npack0),Gy(npack0),Gz(npack0)
      real*8  vlsmooth(npack0)

*     **** local variables ****
      integer k,nxray
      real*8  Q,P,dG

*     **** external functions ****
      real*8   nwpw_splint
      external nwpw_splint

      dG = Gray(3)-Gray(2)
      do k=1,npack0
         Q = dsqrt(Gx(k)**2 + Gy(k)**2 + Gz(k)**2)
         if (Q.gt.1.0d-9) then
            nxray = (Q/dG) + 1.0d0
            P = nwpw_splint(Gray(2),vlray(2,1),vlray(2,2),
     >                      nray-1,nxray-1,Q)
            vlsmooth(k)=P
         else
            vlsmooth(k)=vlray(1,1)
         end if
      end do

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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,n,k,ms,lm,iprj,jprj,shift,indx,indx1
      integer tmp(2)
      integer tid,nthr
      real*8  scal,coeff,w

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

!     *** allocate tmp space ***
      value = BA_push_get(mt_dbl,2*lm_size_max*nthr,'tmp3',
     >                    tmp(2),tmp(1))
      if (.not.value) 
     >   call errquit("nwpw_compcharge_gen_Qlm:cannot push stack",
     >                0,MA_ERR)

      call Parallel_shared_vector_zero(.false.,2*lm_size_max*nthr,
     >                                 dbl_mb(tmp(1)))

      iii = int_mb(ion_iontopaw(1) +ii-1)
      iia = int_mb(katm_iontopaw(1)+ia-1)
      indx = (iii-1)*2*lm_size_max

      !call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)+indx),1)
      call Parallel_shared_vector_zero(.true.,2*lm_size_max,
     >                                 dbl_mb(Qlm(1)+indx))

      scal = 1.0d0/lattice_omega()
      shift = int_mb(shift_Tndiff(1)+iia-1)
      do k=tid+1,int_mb(nindx_Tndiff(1)+iia-1),nthr
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = scal*dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
c            indx1 = indx + (ms-1)*lm_size_max + lm
c            dbl_mb(Qlm(1)+indx1) = dbl_mb(Qlm(1)+indx1) 
c     >                           + coeff*wmatrix(iprj,jprj,ms)
            indx1 = (ms-1)*lm_size_max + lm + tid*2*lm_size_max
            dbl_mb(tmp(1)+indx1) = dbl_mb(tmp(1)+indx1) 
     >                           + coeff*wmatrix(iprj,jprj,ms)
         end do
      end do
!$OMP BARRIER
      call Parallel_shared_vector_reduction(2*lm_size_max,nthr,
     >                                      dbl_mb(tmp(1)),
     >                                      dbl_mb(Qlm(1)+indx))
      !call D1dB_Vector_SumAll(ispin*lm_size_max,dbl_mb(Qlm(1)+indx))

    
!     *** deallocate tmp space ***
      value = BA_pop_stack(tmp(2))
      if (.not.value) 
     >   call errquit("nwpw_compcharge_gen_Qlm:pop stack",1,MA_ERR)
      return
      end 

      subroutine Parallel_shared_vector_reduction(nsize,nthr,A,B)
      implicit none
      integer nsize,nthr
      real*8 A(nsize,nthr)
      real*8 B(nsize)
      integer k

!$OMP DO REDUCTION(+:B)
      do k=1,nthr
         call daxpy(nsize,1.0d0,A(1,k),1,B,1)
      end do
!$OMP END DO
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dE_Qlm(1)),1)
      return
      end


*     *************************************************
*     *                                               *
*     *         nwpw_compcharge_gen_dE_Qlm            *
*     *                                               *
*     *************************************************
*
*   This routine generates dE_Qlm for the non-multipole part of Ecmp_cmp and Ecmp_pw
*
*   dE_Qlm is defined as follows:
*   dE/dQlm =  Sum(G) dcongj(glm(G,ii))*vcmp(G) + dconjg(glm_smooth(G,ii))*vcmp_smooth(G)
*
*   where Qlm = is a function if lm, and ii
*
*   The force wrt to glm and glm_smooth is computed as follows:
*    fion(xyz,ii) = Sum(G) G(xyz)*Qlm(ii)*dimag(glm(G,ii)*dcongj(vcmp(G)))
*                        + G(xyz)*Qlm(ii)*dimag(glm_smooth(G,ii)*dconjg(vcmp_smooth(G)))
*
*    if option%2==1 then include docmp
*    if option>1    then include docmp_smooth
*
      subroutine nwpw_compcharge_gen_dE_Qlm(option,
     >                                      ispin,vcmp,vcmp_smooth,
     >                                      move,fion)
      implicit none
      integer option
      integer ispin
      complex*16 vcmp(*)
      complex*16 vcmp_smooth(*)
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     *** local variables ****
      logical value,docmp,docmp_smooth
      integer np,np_j,np_k,taskid,taskid_j,taskid_k,pcount
      integer iii,iia,ii,ia,lm,l,m,k,jjj,ja,jja,l2,m2,ms
      integer indx1up,indx2up,indx,indx1,indx2
      integer glm_cmp(2),glm_cmp_smooth(2),gvtmp(2),ftmp(2),tmp(2)
      integer Gx,Gy,Gz
      real*8  q22,e1,e2,fourpi
      real*8 e1x,e1y,e1z,e2x,e2y,e2z,qq

*     **** external functions ****
      integer  Pack_G_indx,c_G_indx
      external Pack_G_indx,c_G_indx

      docmp        = mod(option,2).eq.1
      docmp_smooth = option.gt.1
c      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,
c     >           dbl_mb(dE_Qlm(1)),1)
      call Parallel_shared_vector_zero(.true.,2*nion_paw*lm_size_max,
     >                                 dbl_mb(dE_Qlm(1)))

*     **** allocating stack memory ****
      value = BA_push_get(mt_dcpl,npack0,'glm_cmp',
     >                    glm_cmp(2),glm_cmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'glm_cmp_smooth',
     >                    glm_cmp_smooth(2),glm_cmp_smooth(1))
      value = value.and.BA_push_get(mt_dbl,8,'tmpdot',tmp(2),tmp(1))
      call Parallel_shared_vector_zero(.true.,8,dbl_mb(tmp(1)))
      if (move) then
         value = value.and.
     >           BA_push_get(mt_dbl,npack0,'gvtmp',
     >                       gvtmp(2),gvtmp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,3*nion_paw,'ftmp',
     >                       ftmp(2),ftmp(1))
         !call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)
         call Parallel_shared_vector_zero(.true.,3*nion_paw,
     >                                    dbl_mb(ftmp(1)))
      end if
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_dE_Qlm:out of stack',
     >               0,MA_ERR)

      if (isgamma) then

      if (move) then
         Gx = Pack_G_indx(0,1)
         Gy = Pack_G_indx(0,2)
         Gz = Pack_G_indx(0,3)
      end if
     
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         !if (move) call dcopy(npack0,0.0d0,0,dbl_mb(gvtmp(1)),1)
         if (move) call Parallel_shared_vector_zero(.true.,npack0,
     >                                              dbl_mb(gvtmp(1)))

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np_j).eq.taskid_j) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               if (docmp)
     >            call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           dbl_mb(tmp(1)))
               if (docmp_smooth)
     >            call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           dbl_mb(tmp(1)+1))
!$OMP MASTER
               do ms=1,ispin
                  indx  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dE_Qlm(1)+indx) = dbl_mb(tmp(1))
     >                                   + dbl_mb(tmp(1)+1)
               end do
!$OMP END MASTER
!$OMP BARRIER

c               e1 = 0.0d0
c               e2 = 0.0d0
c               if (docmp)
c     >            call Pack_cc_idot(0,
c     >                           dcpl_mb(glm_cmp(1)),
c     >                           vcmp,
c     >                           e1)
c               if (docmp_smooth)
c     >            call Pack_cc_idot(0,
c     >                           dcpl_mb(glm_cmp_smooth(1)),
c     >                           vcmp_smooth,
c     >                           e2)
c               !lm = l*(l+1)+m
c               do ms=1,ispin
c                  indx  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
c                  dbl_mb(dE_Qlm(1)+indx) = (e1+e2)
c               end do

               if (move) then
                  indx1  = (iii-1)*2*lm_size_max+lm
                  indx2  = (iii-1)*2*lm_size_max
     >                   + (ispin-1)*lm_size_max+lm
                  qq = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
                  if (docmp)
     >               call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp(1)),
     >                                       vcmp,
     >                                       dbl_mb(gvtmp(1)))
                  if (docmp_smooth)
     >               call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp_smooth(1)),
     >                                       vcmp_smooth,
     >                                       dbl_mb(gvtmp(1)))
               end if
            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
         if (move) then
            call Pack_tt_idot(0,dbl_mb(Gx),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)))
            call Pack_tt_idot(0,dbl_mb(Gy),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)+1))
            call Pack_tt_idot(0,dbl_mb(Gz),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)+2))

!$OMP MASTER
            dbl_mb(ftmp(1)+3*(iii-1))   = dbl_mb(ftmp(1)+3*(iii-1))  
     >                                  - dbl_mb(tmp(1))
            dbl_mb(ftmp(1)+3*(iii-1)+1) = dbl_mb(ftmp(1)+3*(iii-1)+1)
     >                                  - dbl_mb(tmp(1)+1)
            dbl_mb(ftmp(1)+3*(iii-1)+2) = dbl_mb(ftmp(1)+3*(iii-1)+2)
     >                                  - dbl_mb(tmp(1)+2)
!$OMP END MASTER
!$OMP BARRIER
         end if
      end do



      else

      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)
      taskid = taskid_j + np_j*taskid_k
      np     = np_j+np_k
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np).eq.taskid) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               e1 = 0.0d0
               e2 = 0.0d0
               if (docmp)
     >            call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           e1)
               if (docmp_smooth)
     >            call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           e2)
               !lm = l*(l+1)+m
               do ms=1,ispin
                  indx = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dE_Qlm(1)+indx) = (e1+e2)
               end do

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
      end do

      end if

      if (move) then
         call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
!$OMP DO
         do iii=1,nion_paw
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
            fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
            fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
         end do
!$OMP END DO
      end if

*     **** popping stack memory ****
      value = .true.
      if (move) then
         value = value.and.BA_pop_stack(ftmp(2))
         value = value.and.BA_pop_stack(gvtmp(2))
      end if
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(glm_cmp_smooth(2))
      value = value.and.BA_pop_stack(glm_cmp(2))
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_dE_Qlm:popping stack',
     >               1,MA_ERR)


      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dE_Qlm(1)))

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

#include "bafdecls.fh"
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
            w     = dbl_mb(dE_Qlm(1)+indx1)*coeff
!$OMP DO
            do n=n1(ms),n2(ms)
               sw2(n,iprj) = sw2(n,iprj) + sw1(n,jprj)*w
            end do
!$OMP END DO
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

#include "bafdecls.fh"
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
*     *      nwpw_compcharge_gen_v_cmp_smooth         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_v_cmp_smooth(zv,v_smooth)
      implicit none
      real*8     zv(*)
      complex*16 v_smooth(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,tmp(2),exi(2)


      value = BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value) 
     >  call errquit("nwpw_compcharge_gen_v_cmp_smooth:stack",0,MA_ERR)


      !call dcopy(2*npack0,0.0d0,0,v_smooth,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,v_smooth)
      if (isgamma) then
      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ia   = int_mb(katm_pawtoion(1)+iia-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_tc_Mul(0,dbl_mb(vk_smooth(1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))

         call Pack_cc_daxpy(0,zv(ia),dcpl_mb(tmp(1)),v_smooth)
      end do

      else

      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
         call Cram_rc_Mul(0,dbl_mb(vk_smooth(1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))
         call Cram_cc_daxpy(0,zv(ia),dcpl_mb(tmp(1)),v_smooth)

      end do

      end if

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_v_cmp_smooth:stack",1,MA_ERR)

      return
      end

*     *************************************************
*     *                                               *
*     *      nwpw_compcharge_gen_f_cmp_smooth         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_f_cmp_smooth(zv,dng,fion)
      implicit none
      real*8     zv(*)
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer np_j,taskid_j,G(3)
      integer i,iii,iia,ii,ia,tmp(2),exi(2),xtmp(2),ftmp(2),tmp1(2)
      real*8 fx,fy,fz

      integer  Pack_G_indx
      external Pack_G_indx

      value = BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,npack0,"xtmp",xtmp(2),xtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nion_paw,"ftmp",ftmp(2),ftmp(1))
      value = value.and.BA_push_get(mt_dbl,3,"tmp1",tmp1(2),tmp1(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_f_cmp_smooth:stack",0,MA_ERR)

      !call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion_paw,
     >                                 dbl_mb(ftmp(1)))

      if (isgamma) then
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      G(1)  = Pack_G_indx(0,1)
      G(2)  = Pack_G_indx(0,2)
      G(3)  = Pack_G_indx(0,3)

      do iii=1,nion_paw
       if (mod(iii-1,np_j).eq.taskid_j) then
         iia  = int_mb(katm_paw(1)+iii-1)
         ia   = int_mb(katm_pawtoion(1)+iia-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_tc_Mul(0,dbl_mb(vk_smooth(1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))

         call Pack_cct_iconjgMulb(0,
     >                            dng,
     >                            dcpl_mb(tmp(1)),
     >                            dbl_mb(xtmp(1)))

         call Pack_tt_idot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                     dbl_mb(tmp1(1)))
         call Pack_tt_idot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                     dbl_mb(tmp1(1)+1))
         call Pack_tt_idot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                     dbl_mb(tmp1(1)+2))
!$OMP MASTER
         dbl_mb(ftmp(1)+3*(iii-1))   = dbl_mb(ftmp(1)+3*(iii-1))
     >                               + zv(ia)*dbl_mb(tmp1(1))
         dbl_mb(ftmp(1)+3*(iii-1)+1) = dbl_mb(ftmp(1)+3*(iii-1)+1)
     >                               + zv(ia)*dbl_mb(tmp1(1)+1)
         dbl_mb(ftmp(1)+3*(iii-1)+2) = dbl_mb(ftmp(1)+3*(iii-1)+2)
     >                               + zv(ia)*dbl_mb(tmp1(1)+2)
!$OMP END MASTER
!$OMP BARRIER

       end if
      end do
      call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
!$OMP DO
      do iii=1,nion_paw
         ii = int_mb(ion_pawtoion(1)+iii-1)
         fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
         fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
         fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
      end do
!$OMP END DO

      else
      G(1)  = 0  !*** need to define :(
      G(2)  = 0
      G(3)  = 0

      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
         call Cram_rc_Mul(0,dbl_mb(vk_smooth(1)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))

         do i=1,npack0
           dbl_mb(xtmp(1)+i-1)
     >        = dimag(dng(i))* dble(dcpl_mb(tmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(tmp(1)+i-1))
         end do
         call Cram_rr_dot(0,dbl_mb(G(1)),dbl_mb(xtmp(1)),fx)
         call Cram_rr_dot(0,dbl_mb(G(2)),dbl_mb(xtmp(1)),fy)
         call Cram_rr_dot(0,dbl_mb(G(3)),dbl_mb(xtmp(1)),fz)

         fion(1,ii) = fion(1,ii) + zv(ia)*fx
         fion(2,ii) = fion(2,ii) + zv(ia)*fy
         fion(3,ii) = fion(3,ii) + zv(ia)*fz

      end do

      end if

      value =           BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(ftmp(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_f_cmp_smooth:stack",1,MA_ERR)

      return
      end




*     *************************************************
*     *                                               *
*     *      nwpw_compcharge_gen_vlr_cmp_smooth       *
*     *                                               *
*     *************************************************
*
*  This routine calculates the long-range part of the smooth
*  compensation pseudopotential (used by version4)
*
*   Note - will currently only work with pspw, since aperiodic 
*          boundary conditions are not implemented in BAND
*
      subroutine nwpw_compcharge_gen_vlr_cmp_smooth(zv,
     >                                              rgrid,
     >                                              vlr_smooth)
      implicit none
      real*8 zv(*)
      real*8 rgrid(3,*)
      real*8 vlr_smooth(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer np_j,taskid_j,n2ft3d,n2ft3d_map
      integer i,ii,ia,iia,iii
      real*8  x,y,z,r,b,c

*     **** external functions ****
      logical  control_fast_erf
      external control_fast_erf
      real*8   util_erf,ion_rion
      external util_erf,ion_rion

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)

      !call dcopy(n2ft3d,0.0d0,0,vlr_smooth,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,vlr_smooth)

      b = 1.0d0/dsqrt(4.0d0*datan(1.0d0))   !*** 1/sqrt(pi)
      c = 1.0d0/sigma_smooth

      if (control_fast_erf()) then
         do iii=1,nion_paw
         if (mod(iii-1,np_j).eq.taskid_j) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            ii  = int_mb(ion_pawtoion(1)+iii-1)

            x = ion_rion(1,ii)
            y = ion_rion(2,ii)
            z = ion_rion(3,ii)
!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt((rgrid(1,i)-x)**2
     >                  +(rgrid(2,i)-y)**2
     >                  +(rgrid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                  xerf = r*c
                  yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                  yerf = (1.0d0 - 1.0d0/yerf**4)
c                  yerf = util_erf(xerf)
                  vlr_smooth(i) = vlr_smooth(i) - (zv(ia)/r)*yerf
               else
                  vlr_smooth(i) = vlr_smooth(i) - 2.0d0*zv(ia)*b*c
               end if

            end do
!$OMP END DO
         end if
         end do
      else
         do iii=1,nion_paw
         if (mod(iii-1,np_j).eq.taskid_j) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            ii  = int_mb(ion_pawtoion(1)+iii-1)

            x = ion_rion(1,ii)
            y = ion_rion(2,ii)
            z = ion_rion(3,ii)
!$OMP DO
            do i=1,n2ft3d_map
               r = dsqrt((rgrid(1,i)-x)**2
     >                  +(rgrid(2,i)-y)**2
     >                  +(rgrid(3,i)-z)**2)
               if (r.gt.1.0d-15) then
                  xerf = r*c
                  yerf = util_erf(xerf)
                  vlr_smooth(i) = vlr_smooth(i) - (zv(ia)/r)*yerf
               else
                  vlr_smooth(i) = vlr_smooth(i) - 2.0d0*zv(ia)*b*c
               end if

            end do
!$OMP END DO
         end if
         end do
      end if
      if (np_j.gt.1) call D1dB_Vector_SumAll(n2ft3d_map,vlr_smooth)

      call nwpw_timing_end(5)
      return
      end 

*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_f_lr_cmp_smooth       *
*     *                                               *
*     *************************************************
*   Note - will currently only work with pspw, since aperiodic 
*          boundary conditions are not implemented in BAND
*
      subroutine nwpw_compcharge_gen_f_lr_cmp_smooth(zv,
     >                                               rgrid,
     >                                               rho,fion)
      implicit none
      real*8 zv(*)
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** Error function parameters ****
      real*8 xerf,yerf,fterf,verf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer np_j,taskid_j,n2ft3d,n2ft3d_map
      integer i,ii,ia,iia,iii,np1,np2,np3
      integer ftmp(2)
      real*8  x,y,z,r,b,c,v,dv,q,sqrt_pi
      real*8  rx,ry,rz,fx,fy,fz

*     **** external functions ****
      logical  control_fast_erf
      external control_fast_erf
      real*8   util_erf,ion_rion,lattice_omega
      external util_erf,ion_rion,lattice_omega


*     **** constants ****
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      call nwpw_timing_start(5)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      dv = lattice_omega()/dble(np1*np2*np3)

      if (.not.BA_push_get(mt_dbl,3*nion_paw,'ftmp',ftmp(2),ftmp(1)))
     > call errquit('nwpw_compcharge_gen_f_lr_cmp_smooth:stack',
     >              0,MA_ERR)
      call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)

      b = 1.0d0/dsqrt(4.0d0*datan(1.0d0))   !*** 1/sqrt(pi)
      c = 1.0d0/sigma_smooth

      if (control_fast_erf()) then
!$OMP MASTER
         do iii=1,nion_paw
         if (mod(iii-1,np_j).eq.taskid_j) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            q   = -zv(ia)
            x = ion_rion(1,ii)
            y = ion_rion(2,ii)
            z = ion_rion(3,ii)
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            do i=1,n2ft3d_map
               rx = x - rgrid(1,i)
               ry = y - rgrid(2,i)
               rz = z - rgrid(3,i)
               r  = dsqrt( rx**2 + ry**2 + rz**2)

               if (r .gt. 1.0d-8) then
                 yerf=r*c
                 fterf = (1.0d0
     >                 + yerf*(c1 + yerf*(c2
     >                 + yerf*(c3 + yerf*(c4
     >                 + yerf*(c5 + yerf*c6))))))**4
                 verf = (1.0d0 - 1.0d0/fterf**4)
c                 verf = util_erf(yerf)
                 v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                    - verf)/r**3
               else
                 v = 0.0d0
               end if
               fx = fx + rho(i)*rx*v
               fy = fy + rho(i)*ry*v
               fz = fz + rho(i)*rz*v

            end do
            dbl_mb(ftmp(1)+3*(iii-1))   = -fx*dv
            dbl_mb(ftmp(1)+3*(iii-1)+1) = -fy*dv
            dbl_mb(ftmp(1)+3*(iii-1)+2) = -fz*dv
         end if
         end do
!$OMP END MASTER
!$OMP BARRIER
      else
!$OMP MASTER
         do iii=1,nion_paw
         if (mod(iii-1,np_j).eq.taskid_j) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            q   = -zv(ia)

            x = ion_rion(1,ii)
            y = ion_rion(2,ii)
            z = ion_rion(3,ii)
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            do i=1,n2ft3d_map
               rx = x - rgrid(1,i)
               ry = y - rgrid(2,i)
               rz = z - rgrid(3,i)
               r  = dsqrt( rx**2 + ry**2 + rz**2)

               if (r .gt. 1.0d-8) then
                 yerf=r*c
                 verf = util_erf(yerf)
                 v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                    - verf)/r**3
               else
                 v = 0.0d0
               end if
               fx = fx + rho(i)*rx*v
               fy = fy + rho(i)*ry*v
               fz = fz + rho(i)*rz*v
            end do
            dbl_mb(ftmp(1)+3*(iii-1))   = -fx*dv
            dbl_mb(ftmp(1)+3*(iii-1)+1) = -fy*dv
            dbl_mb(ftmp(1)+3*(iii-1)+2) = -fz*dv
         end if
         end do
!$OMP END MASTER
!$OMP BARRIER
      end if
      call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
!$OMP DO
      do iii=1,nion_paw
         ii = int_mb(ion_pawtoion(1)+iii-1)
         fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
         fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
         fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
      end do
!$OMP END DO

      if (.not.BA_pop_stack(ftmp(2)))
     >call errquit('nwpw_compcharge_gen_f_lr_cmp_smooth:popping stack',
     >             1,MA_ERR)

      call nwpw_timing_end(5)
      return
      end 



*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dn_cmp_smooth_ms      *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_dn_cmp_smooth_ms(ms,dng_smooth)
      implicit none
      integer ms
      complex*16 dng_smooth(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall


      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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


      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp_smooth:stack",1,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dn_cmp_smooth         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_dn_cmp_smooth(ispin,dng_smooth)
      implicit none
      integer ispin
      complex*16 dng_smooth(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall


      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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

      end do
      call Cram_rc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)

      end if


      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall

      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      !call dcopy(2*npack0,0.0d0,0,dng_cmp,1)
      !call dcopy(2*npack0,0.0d0,0,dng_smooth,1)
      call Parallel_shared_vector_zero(.false.,2*npack0,dng_cmp)
      call Parallel_shared_vector_zero(.true., 2*npack0,dng_smooth)

      if (isgamma) then

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         !call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         call Parallel_shared_vector_zero(.true.,2*npack0,
     >                                    dcpl_mb(QTlm(1)))
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

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall


      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp:stack",0,MA_ERR)

      !call dcopy(2*npack0,0.0d0,0,dng_cmp,1)
      call Parallel_shared_vector_zero(.true.,2*npack0,dng_cmp)

      if (isgamma) then

      do iii=1,nion_paw
         iia  = int_mb(katm_paw(1)+iii-1)
         ii = int_mb(ion_pawtoion(1)+iii-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         !call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         call Parallel_shared_vector_zero(.true.,2*npack0,
     >                                    dcpl_mb(QTlm(1)))
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

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = BA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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

c         call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
c     >                          dcpl_mb(QTlm(1)))
c         call Pack_tc_Mul(0,dbl_mb(gk(1)+(iia-1)*npack0),
c     >                       dcpl_mb(tmp1(1)),
c     >                       dcpl_mb(QTlm(1)))
c         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_cmp)

         call Pack_tc_aMulAdd(0,q00,
     >                        dbl_mb(gk(1)+(iia-1)*npack0),
     >                        dcpl_mb(tmp(1)),
     >                        dng_cmp)

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

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_zv_cmp:stack",1,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *       nwpw_compcharge_gen_zv_cmp_smooth       *
*     *                                               *
*     *************************************************
*
*   This routine returns the fourier transform of dn_cmp.

      subroutine nwpw_compcharge_gen_zv_cmp_smooth(zv,cmp_smooth)
      implicit none
      real*8     zv(*)
      complex*16 cmp_smooth(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,ii,ia,lm,l,m,indx1,indx2,tmp(2),exi(2)
      real*8  q,sumall,q00,fourpi

      fourpi = 16.0d0*datan(1.0d0)

      value = BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_zv_cmp_smooth:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,cmp_smooth,1)

      if (isgamma) then

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         lm = 0
         l = 0
         m = 0
         q00 = -zv(ia)/dsqrt(fourpi)

         q = q00
         call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))
         call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),cmp_smooth)
      end do
      call Pack_tc_Mul2(0,dbl_mb(gk_smooth(1)),cmp_smooth)

      else

      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

         lm = 0
         l = 0
         m = 0
         q00 = -zv(ia)/dsqrt(fourpi)

         q = q00
         call Cram_rc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(tmp(1)))
         call Cram_cc_daxpy(0,q,dcpl_mb(tmp(1)),cmp_smooth)
      end do
      call Cram_rc_Mul2(0,dbl_mb(gk_smooth(1)),cmp_smooth)

      end if

      value =           BA_pop_stack(exi(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_zv_cmp_smooth:stack",1,MA_ERR)

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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,lm,tmp(2),exi(2)

      iii = int_mb(ion_iontopaw(1)+ii-1)
      iia  = int_mb(katm_paw(1)+iii-1)
      if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
         lm = l*(l+1) + m

         value = BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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
         value =           BA_pop_stack(exi(2))
         value = value.and.BA_pop_stack(tmp(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer iii,iia,lm,tmp(2),exi(2)

      iii = int_mb(ion_iontopaw(1)+ii-1)
      iia  = int_mb(katm_paw(1)+iii-1)
      if ((l.le.int_mb(mult_l(1)+iia-1)).and.(abs(m).le.l)) then
         lm = l*(l+1) + m

         value = BA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
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

         value =           BA_pop_stack(exi(2))
         value = value.and.BA_pop_stack(tmp(2))
         if (.not.value)
     >     call errquit("nwpw_compcharge_gen_glm2:stack",1,MA_ERR)


      else
         call errquit("nwpw_compcharge_gen_glm2:bad l,m,",l,m)
      end if
      return
      end




********************************* gintegrals tester **************************************


*     ************************************************
*     *                                              *
*     *           nwpw_gintegrals_tester             *
*     *                                              *
*     ************************************************

      subroutine nwpw_gintegrals_tester()
      implicit none

      integer l1,m1,l2,m2,i,j,k
      real*8 W4,sigma
      real*8 R1(3),R2(3),Rab(3)
      complex*16 CW4,CWpp,CW4p,CW4m
      complex*16 CW4pp,CW4pm,CW4mp,CW4mm,CWW

*     **** external functions ****
      real*8   nwpw_WGaussian,nwpw_UGaussian
      external nwpw_WGaussian,nwpw_UGaussian
      complex*16 paw_mult_interaction
      external   paw_mult_interaction
      complex*16 nwpw_CWGaussian
      external   nwpw_CWGaussian

      l1 = 1
      m1 = 0
      l2 = 2
      m2 = 1

      do l1=0,4
      do m1=-l1,l1
      do l2=0,4
      do m2=-l2,l2

c      l1 = 0
c      m1 = 0
c      l2 = 1
c      m2 = -1

      sigma = 6.0d0

      do i= 1,1
      do j= 1,1
      do k= 1,1
      if ((i.ne.0).or.(j.ne.0).or.(k.ne.0)) then

      R1(1) = 0.0d0
      R1(2) = 0.0d0
      R1(3) = 0.0d0

      R2(1) = 5.0d0*i
      R2(2) = 5.0d0*j
      R2(3) = 5.0d0*k

      Rab(1) = R1(1) - R2(1)
      Rab(2) = R1(2) - R2(2)
      Rab(3) = R1(3) - R2(3)

      W4   = nwpw_WGaussian(l1,m1,sigma,l2,m2,sigma,Rab)
      CWpp = nwpw_CWGaussian(l1,m1,sigma,l2,m2,sigma,Rab)
      CWW  = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,R2)

      if ((m1.eq.0).and.(m2.eq.0)) then
         !CW4 = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,Rab)
         CW4 = nwpw_CWGaussian(l1,m1,sigma,l2, m2,sigma,Rab)

      else if ((m1.eq.0).and.(m2.lt.0)) then
         !CW4p = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,R2)
         !CW4m = paw_mult_interaction(l1,m1,sigma,R1,l2,-m2,sigma,R2)
         CW4p = nwpw_CWGaussian(l1,m1,sigma,l2, m2,sigma,Rab)
         CW4m = nwpw_CWGaussian(l1,m1,sigma,l2,-m2,sigma,Rab)
         CW4 = (CW4p - (-1.0d0)**m2 * CW4m)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))

      else if ((m1.eq.0).and.(m2.gt.0)) then
         !CW4p = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,R2)
         !CW4m = paw_mult_interaction(l1,m1,sigma,R1,l2,-m2,sigma,R2)
         CW4p = nwpw_CWGaussian(l1,m1,sigma,l2, m2,sigma,Rab)
         CW4m = nwpw_CWGaussian(l1,m1,sigma,l2,-m2,sigma,Rab)
         CW4 = (CW4m + (-1.0d0)**m2 * CW4p)/dsqrt(2.0d0)

      else if ((m1.lt.0).and.(m2.eq.0)) then
         !CW4p = paw_mult_interaction(l1, m1,sigma,R1,l2,m2,sigma,R2)
         !CW4m = paw_mult_interaction(l1,-m1,sigma,R1,l2,m2,sigma,R2)
         CW4p = nwpw_CWGaussian(l1, m1,sigma,l2,m2,sigma,Rab)
         CW4m = nwpw_CWGaussian(l1,-m1,sigma,l2,m2,sigma,Rab)
         CW4 = (CW4p - (-1.0d0)**m1 * CW4m)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))

      else if ((m1.gt.0).and.(m2.eq.0)) then
         !CW4p = paw_mult_interaction(l1, m1,sigma,R1,l2,m2,sigma,R2)
         !CW4m = paw_mult_interaction(l1,-m1,sigma,R1,l2,m2,sigma,R2)
         CW4p = nwpw_CWGaussian(l1, m1,sigma,l2,m2,sigma,Rab)
         CW4m = nwpw_CWGaussian(l1,-m1,sigma,l2,m2,sigma,Rab)
         CW4 = (CW4m + (-1.0d0)**m1 * CW4p)/dsqrt(2.0d0)

      else
         !CW4pp = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,R2)
         !CW4pm = paw_mult_interaction(l1,m1,sigma,R1,l2,-m2,sigma,R2)
         !CW4mp = paw_mult_interaction(l1,-m1,sigma,R1,l2, m2,sigma,R2)
         !CW4mm = paw_mult_interaction(l1,-m1,sigma,R1,l2,-m2,sigma,R2)
         CW4pp = nwpw_CWGaussian(l1,m1,sigma,l2, m2,sigma,Rab)
         CW4pm = nwpw_CWGaussian(l1,m1,sigma,l2,-m2,sigma,Rab)
         CW4mp = nwpw_CWGaussian(l1,-m1,sigma,l2, m2,sigma,Rab)
         CW4mm = nwpw_CWGaussian(l1,-m1,sigma,l2,-m2,sigma,Rab)
         if ((m1.lt.0).and.(m2.lt.0)) then
            CW4 = -(CW4pp 
     >             + (-1.0d0)**(m1+m2) * CW4mm 
     >             - (-1.0d0)**m2      * CW4pm 
     >             - (-1.0d0)**m1      * CW4mp)/2.0d0

         else if ((m1.lt.0).and.(m2.gt.0)) then
            CW4 = (CW4pm 
     >            - (-1.0d0)**(m1+m2) * CW4mp 
     >            - (-1.0d0)**m1      * CW4mm 
     >            + (-1.0d0)**m2      * CW4pp)*dcmplx(0.0d0,1.0d0/2.0d0)

         else if ((m1.gt.0).and.(m2.lt.0)) then
            CW4 = (CW4mp 
     >            - (-1.0d0)**(m1+m2) * CW4pm 
     >            + (-1.0d0)**m1      * CW4pp 
     >            - (-1.0d0)**m2      * CW4mm)*dcmplx(0.0d0,1.0d0/2.0d0)

         else if ((m1.gt.0).and.(m2.gt.0)) then
            CW4 =  (CW4mm 
     >             + (-1.0d0)**(m1+m2) * CW4pp 
     >             + (-1.0d0)**m1      * CW4pm 
     >             + (-1.0d0)**m2      * CW4mp)/2.0d0
         end if
      end if
      write(*,*) "gintegral tester=",l1,m1,l2,m2,i,j,k 
      write(*,*) "    - W4, CW4   = ",W4,CW4 ,W4-dble(CW4)
      !write(*,*) "    - CWpp,CWW  = ",CWpp,CWW,CWpp-CWW
      write(*,*)


c      CW4p = paw_mult_interaction(l1,m1,sigma,R1,l2, m2,sigma,R2)
c      CW4m = paw_mult_interaction(l1,m1,sigma,R1,l2,-m2,sigma,R2)
c      CW4 = (CW4m + (-1.0d0)**m2 * CW4p) / dsqrt(2.0d0)
c
c      write(*,*) "gintegral tester=",l1,m1,l2,m2,i,j,k 
c      write(*,*) "    - W4, CWpp = ",W4,CWpp,W4-dble(CWpp)
c      write(*,*) "    - W4, CW4p  = ",W4,CW4p
c      write(*,*) "    - W4, CW4m  = ",W4,CW4m
c      write(*,*) "    - W4, CW4   = ",W4,CW4 ,W4-dble(CW4)
c      write(*,*)

      end if
      end do
      end do
      end do

      end do
      end do
      end do
      end do

      return
      end
********************************* gintegrals tester **************************************


********************************* real Gaussian Integrals **************************************
********************************* real Gaussian Integrals **************************************


********************************* complex Gaussian integrals ***********************************

*     ************************************************ 
*     *                                              *
*     *           nwpw_cgintegrals_set_gcount        *
*     *                                              *
*     ************************************************ 

      subroutine nwpw_cgintegrals_set_gcount()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer taskid,np,pcount,gcount,nshl3d,tcount
      integer l1,m1,l2,m2,iii,jjj,iia,jja
      integer tid,nthr

*     **** external functions ****
      integer  control_version,ewald_nshl3d
      external control_version,ewald_nshl3d
      integer  Parallel_maxthreads
      external Parallel_maxthreads

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      nthr = Parallel_maxthreads()
      do l1 = 1,nthr
         int_mb(tcgauss(1)+l1-1) = 0
      end do

      cperiodic = (control_version().eq.3)

      if (cperiodic) then
         nshl3d = ewald_nshl3d()
      else
         nshl3d = 1
      end if

      pcount = 0
      gcount = 0
      tcount = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)

*        **** calculate on-site integrals ****
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
               tid = mod(gcount,nthr)
               int_mb(tcgauss(1)+tid) = int_mb(tcgauss(1)+tid) + 1
               gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     tid = mod(gcount,nthr)
                     int_mb(tcgauss(1)+tid) = int_mb(tcgauss(1)+tid) + 1
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
                     tid = mod(gcount,nthr)
                     int_mb(tcgauss(1)+tid) = int_mb(tcgauss(1)+tid) + 1
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end do
            end do
         end do
      end do
      ncgauss_max = gcount
      ncgauss     = gcount

      tcount = 0
      do l1 = 1,nthr
         int_mb(tgauss_shift(1)+l1-1) = tcount
         tcount = tcount + int_mb(tcgauss(1)+l1-1)
      end do

      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_cgintegrals_init               *
*     *                                               *
*     *************************************************

      subroutine nwpw_cgintegrals_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer nthr

*     **** external functions ****
      integer  Parallel_maxthreads
      external Parallel_maxthreads

      nthr = Parallel_maxthreads()
      value =           BA_alloc_get(mt_int,nthr,"tcgauss",
     >                               tcgauss(2),tcgauss(1))
      value = value.and.BA_alloc_get(mt_int,nthr,"tcgauss_shift",
     >                               tcgauss_shift(2),tcgauss_shift(1))

      call nwpw_cgintegrals_set_gcount()

      value = value.and.BA_alloc_get(mt_int,ncgauss_max,"lm1_cgauss",
     >                              lm1_cgauss(2),lm1_cgauss(1))
      value = value.and.BA_alloc_get(mt_int,ncgauss_max,"lm2_cgauss",
     >                              lm2_cgauss(2),lm2_cgauss(1))
      value = value.and.BA_alloc_get(mt_int,ncgauss_max,"iii1_cgauss",
     >                              iii1_cgauss(2),iii1_cgauss(1))
      value = value.and.BA_alloc_get(mt_int,ncgauss_max,"iii2_cgauss",
     >                              iii2_cgauss(2),iii2_cgauss(1))
      value = value.and.BA_alloc_get(mt_dcpl,ncgauss_max,"e_cgauss",
     >                              e_cgauss(2),e_cgauss(1))
      value = value.and.BA_alloc_get(mt_dcpl,3*ncgauss_max,"f_cgauss",
     >                              f_cgauss(2),f_cgauss(1))
      if (.not.value)
     > call errquit("nwpw_cgintegrals_init:cannot allocate memory",
     >             0,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_cgintegrals_end                *
*     *                                               *
*     *************************************************
      subroutine nwpw_cgintegrals_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value

      value =           BA_free_heap(tcgauss(2))
      value = value.and.BA_free_heap(tcgauss_shift(2))
      value = value.and.BA_free_heap(lm1_cgauss(2))
      value = value.and.BA_free_heap(lm2_cgauss(2))
      value = value.and.BA_free_heap(iii1_cgauss(2))
      value = value.and.BA_free_heap(iii2_cgauss(2))
      value = value.and.BA_free_heap(e_cgauss(2))
      value = value.and.BA_free_heap(f_cgauss(2))
      if (.not.value)
     > call errquit("nwpw_cgintegrals_end:cannot allocate memory",
     >             0,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *             nwpw_cgintegrals_set              *
*     *                                               *
*     *************************************************
c
c  The logic of this routine needs to be completely reworked for threading.
c  It's well designed for MPI parallelism, so one option is to expand all the
c  data structures over tasks and threads instead of just tasks.
c  Another option is to define thread shifts for indx.... However the threshold
c  check with tole would have to be eliminated.
c
      subroutine nwpw_cgintegrals_set(move)
      implicit none
      logical move

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      real*8 tole
      parameter (tole=1.0d-25)

      logical value
      integer taskid,np,pcount,gcount
      integer ii,jj,ia,ja,indx,shft,inds
      integer iii,jjj,iia,jja
      integer lm1,l1,m1,lm2,l2,m2
      integer l,nshl3d,rcell,rcell_hndl
      real*8 R1(3),R12(3),s1,s2,Rab(3),Rba(3),R,ss
      complex*16 W1,W2,W3,W4,dW1(3),dW2(3),dW3(3),dW4(3)
      complex*16 e1,de1(3)

      integer tid,nthr
      integer lm1_tauss(2),lm2_tauss(2),iii1_tauss(2),iii2_tauss(2)
      integer e_tauss(2),f_tauss(2)

*     **** external functions ****
      real*8   ion_rion,nwpw_UGaussian
      external ion_rion,nwpw_UGaussian
      complex*16 nwpw_CWGaussian
      external   nwpw_CWGaussian
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      integer  ewald_nshl3d,ewald_rcell_ptr
      external ewald_nshl3d,ewald_rcell_ptr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(64)
      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      shft = int_mb(tcgauss_shift(1)+tid)

c     **** allocate temporary memory ****
      value =           BA_push_get(mt_int,ncgauss_max,"lm1_tauss",
     >                              lm1_tauss(2),lm1_tauss(1))
      value = value.and.BA_push_get(mt_int,ncgauss_max,"lm2_tauss",
     >                              lm2_tauss(2),lm2_tauss(1))
      value = value.and.BA_push_get(mt_int,ncgauss_max,"iii1_tauss",
     >                              iii1_tauss(2),iii1_tauss(1))
      value = value.and.BA_push_get(mt_int,ncgauss_max,"iii2_tauss",
     >                              iii2_tauss(2),iii2_tauss(1))
      value = value.and.BA_push_get(mt_dcpl,ncgauss_max,"e_tauss",
     >                              e_tauss(2),e_tauss(1))
      value = value.and.BA_push_get(mt_dcpl,3*ncgauss_max,"f_tauss",
     >                              f_tauss(2),f_tauss(1))

      if (cperiodic) then
         nshl3d = ewald_nshl3d()
         rcell  = ewald_rcell_ptr()
      else
         if (.not. BA_push_get(mt_dbl,3,"rcellflm",rcell_hndl,rcell))
     >   call errquit("nwpw_compcharge_set_cgintegrals:stack",1,MA_ERR)

         nshl3d = 1
!$OMP MASTER
         dbl_mb(rcell)   = 0.0d0
         dbl_mb(rcell+1) = 0.0d0
         dbl_mb(rcell+2) = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      end if


      !**** these should not need to be called!!!! ****
!$OMP MASTER
      call dcopy(2*ncgauss_max,0.0d0,0,dcpl_mb(e_tauss(1)),1)
      call dcopy(2*3*ncgauss_max,0.0d0,0,dcpl_mb(f_tauss(1)),1)
      call dcopy(2*ncgauss_max,0.0d0,0,dcpl_mb(e_cgauss(1)),1)
      call dcopy(2*3*ncgauss_max,0.0d0,0,dcpl_mb(f_cgauss(1)),1)
!$OMP END MASTER
!$OMP BARRIER


      pcount = 0
      gcount = 0
      indx   = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)
         s1  = dbl_mb(sigma_paw(1)+iia-1)

*        **** calculate on-site integrals ****
         lm1 = 0
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
            if (mod(gcount,nthr).eq.tid) then
                W1=nwpw_UGaussian(l1,m1,s1,l1,m1,s1)
                W2=nwpw_UGaussian(l1,m1,s1,l1,m1,sigma_smooth)
                W4=nwpw_UGaussian(l1,m1,sigma_smooth,l1,m1,sigma_smooth)
                e1 = 0.5d0*W1 + 0.5d0*W4 - W2
                inds = indx + shft
                dcpl_mb(e_tauss(1)+inds) = e1
                int_mb(lm1_tauss(1)+inds) = (iii-1)*2*lm_size_max+lm1
                int_mb(lm2_tauss(1)+inds) = (iii-1)*2*lm_size_max+lm1
                int_mb(iii1_tauss(1)+inds) = iii
                int_mb(iii2_tauss(1)+inds) = iii

                indx = indx + 1
c               !end if
            end if
            gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               lm2 = 0
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                  if (mod(gcount,nthr).eq.tid) then
                     e1 = dcmplx(0.0d0,0.0d0)
                     do l=2,nshl3d
                        Rab(1) = dbl_mb(rcell+l-1)
                        Rab(2) = dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                           W1 = nwpw_CWGaussian(l1,m1,s1,l2,m2,s1,Rab)
                           W2 = nwpw_CWGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
c                           W3 = nwpw_CWGaussian(l1,m1,sigma_smooth,
c     >                                         l2,m2,s1,Rab)
                           W4 = nwpw_CWGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           !e1=e1+0.5d0*W1+0.5d0*W4-0.5d0*W2-0.5d0*W3
                           !e1=e1+(W1+W4-2.0d0*W2)
                           !e1=e1+(W1+W4-W2-W3)
                           e1=e1+0.5d0*W1+0.5d0*W4-W2
                        end if
                     end do

                     inds = indx + shft
                     dcpl_mb(e_tauss(1)+inds) = e1
                     int_mb(lm1_tauss(1)+inds) 
     >                     = (iii-1)*2*lm_size_max+lm1
                     int_mb(lm2_tauss(1)+inds)
     >                     = (iii-1)*2*lm_size_max+lm2
                     int_mb(iii1_tauss(1)+inds) = iii
                     int_mb(iii2_tauss(1)+inds) = iii

                     indx = indx + 1
                     !end if
                  end if
                  gcount = gcount + 1
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
                  if (mod(gcount,nthr).eq.tid) then
                     e1 = dcmplx(0.0d0,0.0d0)
                     de1(1) = dcmplx(0.0d0,0.0d0)
                     de1(2) = dcmplx(0.0d0,0.0d0)
                     de1(3) = dcmplx(0.0d0,0.0d0)
                     do l=1,nshl3d
                        Rab(1) = R12(1) + dbl_mb(rcell+l-1)
                        Rab(2) = R12(2) + dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = R12(3) + dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                        if (move) then
                           call nwpw_dCWGaussian(l1,m1,s1,l2,m2,s2,Rab,
     >                                          W1,dW1)
                           call nwpw_dCWGaussian(l1,m1,s1,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W2,dW2)
                           call nwpw_dCWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,s2,Rab,
     >                                          W3,dW3)
                           call nwpw_dCWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W4,dW4)
                           e1 = e1 + (W1 + W4 - W2 - W3)
                           de1(1)=de1(1)+(dW1(1)+dW4(1)-dW2(1)-dW3(1))
                           de1(2)=de1(2)+(dW1(2)+dW4(2)-dW2(2)-dW3(2))
                           de1(3)=de1(3)+(dW1(3)+dW4(3)-dW2(3)-dW3(3))

                        else
                           W1 = nwpw_CWGaussian(l1,m1,s1,l2,m2,s2,Rab)
                           W2 = nwpw_CWGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
                           W3 = nwpw_CWGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,s2,Rab)
                           W4 = nwpw_CWGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           e1 = e1 + (W1 + W4 - W2 - W3)
                        end if
                        end if
                     end do
                     
                     inds = indx + shft
                     dcpl_mb(e_tauss(1)+inds) = e1
                     if (move) then
                        dcpl_mb(f_tauss(1)+3*inds)   = de1(1)
                        dcpl_mb(f_tauss(1)+3*inds+1) = de1(2)
                        dcpl_mb(f_tauss(1)+3*inds+2) = de1(3)
                     end if
                     int_mb(lm1_tauss(1)+inds)
     >                  = (iii-1)*2*lm_size_max+lm1
                     int_mb(lm2_tauss(1)+inds)
     >                  = (jjj-1)*2*lm_size_max+lm2
                     int_mb(iii1_tauss(1)+inds) = iii
                     int_mb(iii2_tauss(1)+inds) = jjj

                     indx = indx + 1
                     !end if
                  end if
                  gcount = gcount + 1
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

!$OMP BARRIER
!$OMP MASTER
      call nwpw_cgintegral_stripper(ncgauss_max,
     >                             int_mb(iii1_tauss(1)),
     >                             int_mb(iii2_tauss(1)),
     >                             int_mb(lm1_tauss(1)),
     >                             int_mb(lm2_tauss(1)),
     >                             dcpl_mb(e_tauss(1)),
     >                             dcpl_mb(f_tauss(1)),
     >                             ncgauss,
     >                             int_mb(iii1_cgauss(1)),
     >                             int_mb(iii2_cgauss(1)),
     >                             int_mb(lm1_cgauss(1)),
     >                             int_mb(lm2_cgauss(1)),
     >                             dcpl_mb(e_cgauss(1)),
     >                             dcpl_mb(f_cgauss(1)))
!$OMP END MASTER
!$OMP BARRIER
c Need to have barrier before deallocation below.  This barrier could be removed
c if the tauss variables were on the heap and not deallocated rather than stack.


c     **** deallocate stack memory ****
      if (.not.periodic) then
        if (.not.BA_pop_stack(rcell_hndl))
     >   call errquit("nwpw_compcharge_set_cgintegrals:stack",2,MA_ERR)
      end if
      value =           BA_pop_stack(f_tauss(2))
      value = value.and.BA_pop_stack(e_tauss(2))
      value = value.and.BA_pop_stack(iii2_tauss(2))
      value = value.and.BA_pop_stack(iii1_tauss(2))
      value = value.and.BA_pop_stack(lm2_tauss(2))
      value = value.and.BA_pop_stack(lm1_tauss(2))
        if (.not.value)
     >   call errquit("nwpw_compcharge_set_cgintegrals:stack",3,MA_ERR)

      call nwpw_timing_end(64)
      return
      end

*     *************************************************
*     *                                               *
*     *            nwpw_cgintegrals_stripper          *
*     *                                               *
*     *************************************************
c
c  This routine is used to remove unecessary integrals
c
      subroutine nwpw_cgintegral_stripper(ng_in,
     >                                   iii1_in,iii2_in,
     >                                   lm1_in,lm2_in,e_in,f_in,
     >                                   ng_out,
     >                                   iii1_out,iii2_out,
     >                                   lm1_out,lm2_out,e_out,f_out)
      implicit none
      integer    ng_in,iii1_in(*),iii2_in(*),lm1_in(*),lm2_in(*)
      complex*16 e_in(*),f_in(*)
      integer    ng_out,iii1_out(*),iii2_out(*),lm1_out(*),lm2_out(*)
      complex*16 e_out(*),f_out(*)

c     **** local variables ****
      real*8 tole
      parameter (tole=1.0d-25)
      integer i

      ng_out = 0
      do i=1,ng_in
         if (Abs(e_in(i)).gt.tole) then
            ng_out = ng_out + 1
            iii1_out(ng_out) = iii1_in(i)
            iii2_out(ng_out) = iii2_in(i)
            lm1_out(ng_out)  = lm1_in(i)
            lm2_out(ng_out)  = lm2_in(i)
            e_out(ng_out)    = e_in(i)
            f_out(3*(ng_out-1)+1) = f_in(3*(i-1)+1)
            f_out(3*(ng_out-1)+2) = f_in(3*(i-1)+2)
            f_out(3*(ng_out-1)+3) = f_in(3*(i-1)+3)
         end if
      end do
      return
      end

*     *************************************************
*     *                                               *
*     *            nwpw_cgintegrals_unscramble        *
*     *                                               *
*     *************************************************

      subroutine nwpw_cintegrals_unscramble(ng_in,
     >                                   uindex,ucoeff,
     >                                   iii1_in,iii2_in,
     >                                   lm1_in,lm2_in,e_in,f_in,
     >                                   ng_out,
     >                                   iii1_out,iii2_out,
     >                                   lm1_out,lm2_out,e_out,f_out)
      implicit none
      integer    ng_in
      integer    uindex(4,*)
      real*8     ucoeff(4,*)
      integer iii1_in(*),iii2_in(*),lm1_in(*),lm2_in(*)
      complex*16 e_in(*),f_in(*)
      integer    ng_out,iii1_out(*),iii2_out(*),lm1_out(*),lm2_out(*)
      complex*16 e_out(*),f_out(*)

c     **** local variables ****
      real*8 tole
      parameter (tole=1.0d-25)
      integer i
      real*8 etmp,ftmpx,ftmpy,ftmpz

      ng_out = 0
      do i=1,ng_in
         etmp = ucoeff(1,i)*e_in(uindex(1,i))
     >        + ucoeff(2,i)*e_in(uindex(2,i))
     >        + ucoeff(3,i)*e_in(uindex(3,i))
     >        + ucoeff(4,i)*e_in(uindex(4,i))
         ftmpx = ucoeff(1,i)*f_in(3*(uindex(1,i)-1)+1)
     >         + ucoeff(2,i)*f_in(3*(uindex(2,i)-1)+1)
     >         + ucoeff(3,i)*f_in(3*(uindex(3,i)-1)+1)
     >         + ucoeff(4,i)*f_in(3*(uindex(4,i)-1)+1)
         ftmpy = ucoeff(1,i)*f_in(3*(uindex(1,i)-1)+2)
     >         + ucoeff(2,i)*f_in(3*(uindex(2,i)-1)+2)
     >         + ucoeff(3,i)*f_in(3*(uindex(3,i)-1)+2)
     >         + ucoeff(4,i)*f_in(3*(uindex(4,i)-1)+2)
         ftmpz = ucoeff(1,i)*f_in(3*(uindex(1,i)-1)+3)
     >         + ucoeff(2,i)*f_in(3*(uindex(2,i)-1)+3)
     >         + ucoeff(3,i)*f_in(3*(uindex(3,i)-1)+3)
     >         + ucoeff(4,i)*f_in(3*(uindex(4,i)-1)+3)
         if (Abs(etmp).gt.tole) then
            ng_out = ng_out + 1
            iii1_out(ng_out) = iii1_in(i)
            iii2_out(ng_out) = iii2_in(i)
            lm1_out(ng_out)  = lm1_in(i)
            lm2_out(ng_out)  = lm2_in(i)
            e_out(ng_out)    = etmp
            f_out(3*(ng_out-1)+1) = ftmpx
            f_out(3*(ng_out-1)+2) = ftmpy
            f_out(3*(ng_out-1)+3) = ftmpz
         end if
      end do

      return
      end 

********************************* complex Gaussian integrals **************************************


*     *****************************************************
*     *                                                   *
*     *        nwpw_compcharge_E_multipole_Zv_ee          *
*     *                                                   *
*     *****************************************************
*
*    This routine calulates the multipole energy defined by 
* 
*                    //
*    Emultipole  =  || (nelccmp(r)-nelccmp_smooth(r))*(nioncmp(r')-nioncmp_smooth(r'))
*                   || ---------------------------------------------------------------  dr dr'
*                   //                       |r-r'|
*
      real*8 function nwpw_compcharge_E_multipole_zv_ee(ispin,zv)
      implicit none
      integer ispin
      real*8 zv(*)


#include "bafdecls.fh"
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
            q11 = -1.0d0*zv(ia)/dsqrt(fourpi)
         else
            q11 = 0.0d0
         end if

         if (l2.eq.0) then
            jja = int_mb(katm_paw(1)+jjj-1)
            ja  = int_mb(katm_pawtoion(1)+jja-1)
            q22 = -1.0d0*zv(ja)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if

         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)
         E = E + (q11*q2 + q1*q22)*dbl_mb(e_gauss(1)+k-1)

      end do
      call Parallel_SumAll(E)


      nwpw_compcharge_E_multipole_Zv_ee = E
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 E,q1,q2
      integer iii,jjj,l1,m1,l2,m2,lm

      integer iic,jjc,lic,ljc
c      real*8  etest,etest2
      common /nwpw_compcharge_E_mul_tmp/ E

c      do iii=1,nion_paw
c         lm = 1
c         do l1=0,int_mb(mult_l(1))
c         do m1=-l1,l1
c          write(*,*) "ii,lm,Qlm=",iii,lm,
c     >    dbl_mb(Qlm(1)+(iii-1)*2*lm_size_max+lm-1)
c     >  +dbl_mb(Qlm(1)+(iii-1)*2*lm_size_max+(ispin-1)*lm_size_max+lm-1)
c          lm = lm + 1
c         end do
c         end do
c      end do
c      do iic=1,nion_paw
c      do jjc=1,nion_paw
c      etest2=0.0d0
c      do lic=0,2
c      do ljc=0,2
c      etest = 0.0d0
c      !iic = 1
c      !jjc = 1
c      !lic = 1
c      !ljc = 1

!$OMP MASTER
      E = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:E)
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max

         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)
         E = E + q1*q2*dbl_mb(e_gauss(1)+k-1)

c         call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
c         call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)
c         etest= q1*q2*dbl_mb(e_gauss(1)+k-1)
c
cc         if (dabs(etest).gt.1.0d-9) 
c     >    write(24,'(A,6I3,4E14.6)') "ii,li,mi,jj,lj,mj=", 
c     >      iii,l1,m1,jjj,l2,m2,q1,q2,etest,
c     >      dbl_mb(e_gauss(1)+k-1)
c
cc         if ((iic.eq.iii).and.(lic.eq.l1).and.
cc     >       (jjc.eq.jjj).and.(ljc.eq.l2)) then
cc            etest = etest + q1*q2*dbl_mb(e_gauss(1)+k-1)
cc         end if

      end do
!$OMP END DO


      call Parallel_SumAll(E)


      nwpw_compcharge_E_multipole = E
      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_F_multipole_all     *
*     *                                               *
*     *************************************************
*
*  This routine computes the multiple only contribution to the force.
*  Note that the components of the force due to the variation wrt ncmp
*  are included in the standard non-local psp forces.
*
      subroutine nwpw_compcharge_F_multipole_all(ispin,zv,fion)
      implicit none
      integer ispin
      real*8 zv(*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      integer iii1,iii2,ii1,ii2
      integer iii,l1,m1,iia,ia
      integer jjj,l2,m2,jja,ja
      integer fmult(2)
      real*8 q1,q1x,q1y,q1z,q2,q2x,q2y,q2z
      real*8 q11,q22,fourpi
      real*8 tmp,tmpx,tmpy,tmpz,R,R12(3)

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      fourpi = 16.0d0*datan(1.0d0)

      if (.not.BA_push_get(mt_dbl,3*nion_paw,"fmult",fmult(2),fmult(1)))
     >   call errquit("nwpw_compcharge_F_multipole:stack",1,MA_ERR)
      !call dcopy(3*nion_paw,0.0d0,0,dbl_mb(fmult(1)),1)
      call Parallel_shared_vector_zero(.true.,3*nion_paw,
     >                                 dbl_mb(fmult(1)))

!$OMP MASTER
      do k=1,ngauss
         iii1 = int_mb(iii1_gauss(1)+k-1)
         iii2 = int_mb(iii2_gauss(1)+k-1)
         if (iii1.ne.iii2) then
            indx1up = int_mb(lm1_gauss(1)+k-1)
            indx1dn = indx1up + (ispin-1)*lm_size_max
            indx2up = int_mb(lm2_gauss(1)+k-1)
            indx2dn = indx2up + (ispin-1)*lm_size_max

            call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
            call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)

            if (l1.eq.0) then
               iia = int_mb(katm_paw(1)+iii-1)
               ia  = int_mb(katm_pawtoion(1)+iia-1)
               q11 = -1.0d0*zv(ia)/dsqrt(fourpi)
            else
               q11 = 0.0d0
            end if
            if (l2.eq.0) then
               jja = int_mb(katm_paw(1)+jjj-1)
               ja  = int_mb(katm_pawtoion(1)+jja-1)
               q22 = -1.0d0*zv(ja)/dsqrt(fourpi)
            else
               q22 = 0.0d0
            end if

            q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
            q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)

c            ii1 = int_mb(ion_pawtoion(1)+iii1-1)
c            ii2 = int_mb(ion_pawtoion(1)+iii2-1)
c            R12(1) = ion_rion(1,ii1) - ion_rion(1,ii2)
c            R12(2) = ion_rion(2,ii1) - ion_rion(2,ii2)
c            R12(3) = ion_rion(3,ii1) - ion_rion(3,ii2)
c            R = dsqrt(R12(1)*R12(1)+R12(2)*R12(2)+R12(3)*R12(3))
c            R12(1) = R12(1)/R
c            R12(2) = R12(2)/R
c            R12(3) = R12(3)/R

            tmp  = (q1*q2 + q11*q2 + q1*q22)
            !tmp  = (q11*q2 + q1*q22)

            call daxpy(3,-tmp,dbl_mb(f_gauss(1)+3*(k-1)), 1,
     >                        dbl_mb(fmult(1)+3*(iii1-1)),1)
            call daxpy(3, tmp,dbl_mb(f_gauss(1)+3*(k-1)), 1,
     >                        dbl_mb(fmult(1)+3*(iii2-1)),1)

c            tmpx = tmp*dbl_mb(f_gauss(1)+3*(k-1))
c            tmpy = tmp*dbl_mb(f_gauss(1)+3*(k-1)+1)
c            tmpz = tmp*dbl_mb(f_gauss(1)+3*(k-1)+2)
c            dbl_mb(fmult(1)+3*(iii1-1))
c     >         = dbl_mb(fmult(1)+3*(iii1-1))   - tmpx 
c            dbl_mb(fmult(1)+3*(iii1-1)+1)
c     >         = dbl_mb(fmult(1)+3*(iii1-1)+1) - tmpy 
c            dbl_mb(fmult(1)+3*(iii1-1)+2)
c     >         = dbl_mb(fmult(1)+3*(iii1-1)+2) - tmpz 
c
c            dbl_mb(fmult(1)+3*(iii2-1))  =dbl_mb(fmult(1)+3*(iii2-1))
c     >                                   +tmpx
c            dbl_mb(fmult(1)+3*(iii2-1)+1)=dbl_mb(fmult(1)+3*(iii2-1)+1)
c     >                                   +tmpy
c            dbl_mb(fmult(1)+3*(iii2-1)+2)=dbl_mb(fmult(1)+3*(iii2-1)+2)
c     >                                   +tmpz 
         end if
      end do
!$OMP END MASTER
!$OMP BARRIER
      call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(fmult(1)))

!$OMP DO
      do iii1=1,nion_paw
         ii1 = int_mb(ion_pawtoion(1)+iii1-1)
         fion(1,ii1) = fion(1,ii1) + dbl_mb(fmult(1)+3*(iii1-1))
         fion(2,ii1) = fion(2,ii1) + dbl_mb(fmult(1)+3*(iii1-1)+1)
         fion(3,ii1) = fion(3,ii1) + dbl_mb(fmult(1)+3*(iii1-1)+2)
      end do
!$OMP END DO

      if (.not.BA_pop_stack(fmult(2)))
     >  call errquit("nwpw_compcharge_F_multipole_all:stack",2,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_F_multipole         *
*     *                                               *
*     *************************************************
*
      subroutine nwpw_compcharge_F_multipole(ispin,fion)
      implicit none
      integer ispin
      real*8 fion(*)

#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_dbl,3*nion_paw,"fmult",fmult(2),fmult(1)))
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


      if (.not.BA_pop_stack(fmult(2)))
     >  call errquit("nwpw_compcharge_F_multipole:stack",2,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dE_ncmp_vloc_Qlm      *
*     *                                               *
*     *************************************************
*
*    This routine calulates dE_ncmp_vloc/dQlm where
*  
*                  /
*   E_ncmp_vloc  = | (ncmp(r)*Vl2(r) + ncmp_smooth(r)*(Vl1(r)-Vl2(r))) dr - plane-wave integrals
*                  /
*
*                  /
*                + | (ncmp(r)-ncmp_smooth(r))*(Vl1(r)-Vl2(r)) dr  - Gaussian two-center integrals
*                  /
*  
*   vcmp = V2l and vcmp_smooth=Vl1-Vl2
*

      subroutine nwpw_compcharge_gen_dE_ncmp_vloc_Qlm(ispin,zv,
     >                                           vcmp,vcmp_smooth,
     >                                           move,fion)
      implicit none
      integer ispin
      real*8     zv(*)
      complex*16 vcmp(*)
      complex*16 vcmp_smooth(*)
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     *** local variables ****
      logical value
      integer np,np_j,np_k,taskid,taskid_j,taskid_k,pcount
      integer iii,iia,ii,ia,lm,l,m,k,jjj,ja,jja,l1,m1,l2,m2,ms
      integer indx1up,indx2up,indx1,indx2,Gx,Gy,Gz
      integer glm_cmp(2),glm_cmp_smooth(2),gvtmp(2),ftmp(2),tmp(2)
      real*8  q11,q22,e1,e2,fourpi
      real*8   qq,e1x,e1y,e1z

      integer  Pack_G_indx,c_G_indx
      external Pack_G_indx,c_G_indx

      call nwpw_timing_start(38)

c      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,
c     >           dbl_mb(dElocal_Qlm(1)),1)
      call Parallel_shared_vector_zero(.true.,2*nion_paw*lm_size_max,
     >                                 dbl_mb(dElocal_Qlm(1)))

*     **** allocating stack memory ****
      value = BA_push_get(mt_dcpl,npack0,'glm_cmp',
     >                    glm_cmp(2),glm_cmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'glm_cmp_smooth',
     >                    glm_cmp_smooth(2),glm_cmp_smooth(1))
      value = value.and.BA_push_get(mt_dbl,30,'tmpdot',tmp(2),tmp(1))
      if (move) then
         value = value.and.
     >           BA_push_get(mt_dbl,npack0,'gvtmp',
     >                       gvtmp(2),gvtmp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,3*nion_paw,'ftmp',
     >                       ftmp(2),ftmp(1))
         !call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)
         call Parallel_shared_vector_zero(.true.,3*nion_paw,
     >                                    dbl_mb(ftmp(1)))
      end if
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:out of stack',
     >               0,MA_ERR)

      if (isgamma) then
      if (move) then
         Gx = Pack_G_indx(0,1)
         Gy = Pack_G_indx(0,2)
         Gz = Pack_G_indx(0,3)
      end if
     
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         !if (move) call dcopy(npack0,0.0d0,0,dbl_mb(gvtmp(1)),1)
         if (move) 
     >      call Parallel_shared_vector_zero(.true.,npack0,
     >                                       dbl_mb(gvtmp(1)))

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np_j).eq.taskid_j) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           dbl_mb(tmp(1)))
               call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           dbl_mb(tmp(1)+1))
               !lm = l*(l+1)+m
!$OMP MASTER
               do ms=1,ispin
                  indx1  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx1) 
     >               = (dbl_mb(tmp(1))+dbl_mb(tmp(1)+1))
               end do
!$OMP END MASTER
!$OMP BARRIER

               if (move) then
                  indx1  = (iii-1)*2*lm_size_max+lm
                  indx2  = (iii-1)*2*lm_size_max
     >                   + (ispin-1)*lm_size_max+lm
                  qq = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
                  call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp(1)),
     >                                       vcmp,
     >                                       dbl_mb(gvtmp(1)))
                  call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp_smooth(1)),
     >                                       vcmp_smooth,
     >                                       dbl_mb(gvtmp(1)))
               end if

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
         if (move) then
            call Pack_tt_idot(0,dbl_mb(Gx),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)))
            call Pack_tt_idot(0,dbl_mb(Gy),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)+1))
            call Pack_tt_idot(0,dbl_mb(Gz),dbl_mb(gvtmp(1)),
     >                        dbl_mb(tmp(1)+2))
!$OMP MASTER
            dbl_mb(ftmp(1)+3*(iii-1))   = dbl_mb(ftmp(1)+3*(iii-1))  
     >                                  - dbl_mb(tmp(1))
            dbl_mb(ftmp(1)+3*(iii-1)+1) = dbl_mb(ftmp(1)+3*(iii-1)+1)
     >                                  - dbl_mb(tmp(1)+1)
            dbl_mb(ftmp(1)+3*(iii-1)+2) = dbl_mb(ftmp(1)+3*(iii-1)+2)
     >                                  - dbl_mb(tmp(1)+2)
!$OMP END MASTER
!$OMP BARRIER


         end if
      end do !** iii **


      else

      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)
      taskid = taskid_j + np_j*taskid_k
      np     = np_j+np_k
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np).eq.taskid) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           e1)
               call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           e2)
               !lm = l*(l+1)+m
               do ms=1,ispin
                  indx1 = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx1) = (e1+e2)
               end do

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
      end do !** iii **

      end if

      if (move) then
         call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
!$OMP DO
         do iii=1,nion_paw
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
            fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
            fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
         end do
!$OMP END DO
      end if

*     **** popping stack memory ****
      value = .true.
      if (move) then
         value = value.and.BA_pop_stack(ftmp(2))
         value = value.and.BA_pop_stack(gvtmp(2))
      end if
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(glm_cmp_smooth(2))
      value = value.and.BA_pop_stack(glm_cmp(2))
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:popping stack',
     >               1,MA_ERR)

      fourpi = 16.0d0*datan(1.0d0)
!$OMP MASTER
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx2up = int_mb(lm2_gauss(1)+k-1)

         call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
         call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)

         if (l1.eq.0) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            q11 = -1.0d0*zv(ia)/dsqrt(fourpi)
         else
            q11 = 0.0d0
         end if

         if (l2.eq.0) then
            jja = int_mb(katm_paw(1)+jjj-1)
            ja  = int_mb(katm_pawtoion(1)+jja-1)
            q22 = -1.0d0*zv(ja)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if

         !q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         !q2 = q22
         !E = E + q1*q2*dbl_mb(e_gauss(1)+k-1)
      
         do ms=1,ispin
            indx1 = indx1up + (ms-1)*lm_size_max
            indx2 = indx2up + (ms-1)*lm_size_max
            dbl_mb(dElocal_Qlm(1)+indx1) = dbl_mb(dElocal_Qlm(1)+indx1)
     >                                   + q22*dbl_mb(e_gauss(1)+k-1)
            dbl_mb(dElocal_Qlm(1)+indx2) = dbl_mb(dElocal_Qlm(1)+indx2)
     >                                   + q11*dbl_mb(e_gauss(1)+k-1)
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER
      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dElocal_Qlm(1)))

      call nwpw_timing_end(38)
      return
      end


*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw   *
*     *                                               *
*     *************************************************
*
*    This routine calulates dE_ncmp_vloc/dQlm where
*  
*                  /
*   E_ncmp_vloc  = | ncmp(r)*Vl2(r) dr - plane-wave integrals
*                  /
*
*   vcmp = V2l 
*

      subroutine nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_pw(ispin,vcmp,
     >                                                   move,fion)
      implicit none
      integer ispin
      complex*16 vcmp(*)
      logical move
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     *** local variables ****
      logical value
      integer np,np_j,np_k,taskid,taskid_j,taskid_k,pcount
      integer iii,iia,ii,ia,lm,l,m,k,jjj,ja,jja,l2,m2,ms
      integer indx,indx1,indx2,Gx,Gy,Gz
      integer glm_cmp(2),glmvcmp(2),ftmp(2),tmp(2)
      real*8  e1,qq,e1x,e1y,e1z
      integer  Pack_G_indx,c_G_indx
      external Pack_G_indx,c_G_indx

      call nwpw_timing_start(38)

c      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,
c     >           dbl_mb(dElocal_Qlm(1)),1)
      call Parallel_shared_vector_zero(.true.,2*nion_paw*lm_size_max,
     >                                 dbl_mb(dElocal_Qlm(1)))

*     **** allocating stack memory ****
      value = BA_push_get(mt_dcpl,npack0,'glm_cmp',
     >                    glm_cmp(2),glm_cmp(1))
      value = value.and.BA_push_get(mt_dbl,3,'tmpdot',tmp(2),tmp(1))
      if (move) then
         value = value.and.BA_push_get(mt_dbl,npack0,'glmvcmp',
     >                                 glmvcmp(2),glmvcmp(1))
         value = value.and.BA_push_get(mt_dbl,3*nion_paw,'ftmp',
     >                                 ftmp(2),ftmp(1))
c         call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)
         call Parallel_shared_vector_zero(.true.,3*nion_paw,
     >                                    dbl_mb(ftmp(1)))
      end if
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:out of stack',
     >               0,MA_ERR)

      if (isgamma) then
      if (move) then
         Gx = Pack_G_indx(0,1)
         Gy = Pack_G_indx(0,2)
         Gz = Pack_G_indx(0,3)
      end if
     
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         !if (move) call dcopy(npack0,0.0d0,0,dbl_mb(glmvcmp(1)),1)
         if (move) call Parallel_shared_vector_zero(.true.,npack0,
     >                                              dbl_mb(glmvcmp(1)))

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np_j).eq.taskid_j) then
               call nwpw_compcharge_gen_glm(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)))
               call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           dbl_mb(tmp(1)))
               !lm = l*(l+1)+m
!$OMP MASTER
               do ms=1,ispin
                  indx  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx) = dbl_mb(tmp(1))
               end do
!$OMP END MASTER
!$OMP BARRIER

               if (move) then
                  indx1  = (iii-1)*2*lm_size_max+lm
                  indx2  = (iii-1)*2*lm_size_max
     >                   + (ispin-1)*lm_size_max+lm
                  qq = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
                  call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp(1)),
     >                                       vcmp,
     >                                       dbl_mb(glmvcmp(1)))
               end if
            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
         if (move) then
            call Pack_tt_idot(0,dbl_mb(Gx),dbl_mb(glmvcmp(1)),
     >                        dbl_mb(tmp(1)))
            call Pack_tt_idot(0,dbl_mb(Gy),dbl_mb(glmvcmp(1)),
     >                        dbl_mb(tmp(1)+1))
            call Pack_tt_idot(0,dbl_mb(Gz),dbl_mb(glmvcmp(1)),
     >                        dbl_mb(tmp(1)+2))
!$OMP MASTER
            dbl_mb(ftmp(1)+3*(iii-1))   = dbl_mb(ftmp(1)+3*(iii-1))  
     >                                  - dbl_mb(tmp(1))
            dbl_mb(ftmp(1)+3*(iii-1)+1) = dbl_mb(ftmp(1)+3*(iii-1)+1)
     >                                  - dbl_mb(tmp(1)+1)
            dbl_mb(ftmp(1)+3*(iii-1)+2) = dbl_mb(ftmp(1)+3*(iii-1)+2)
     >                                  - dbl_mb(tmp(1)+2)
!$OMP END MASTER
!$OMP BARRIER

         end if
      end do


      else

      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)
      taskid = taskid_j + np_j*taskid_k
      np     = np_j+np_k
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np).eq.taskid) then
               call nwpw_compcharge_gen_glm(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)))
               call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           e1)
               !lm = l*(l+1)+m
               do ms =1,ispin
                  indx  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx) = e1
               end do

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
      end do

      end if

      if (move) then
         call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
!$OMP DO
         do iii=1,nion_paw
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
            fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
            fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
         end do
!$OMP END DO
      end if

*     **** popping stack memory ****
      value = .true.
      if (move) then
         value = value.and.BA_pop_stack(ftmp(2))
         value = value.and.BA_pop_stack(glmvcmp(2))
      end if
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(glm_cmp(2))
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:popping stack',
     >               1,MA_ERR)


      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dElocal_Qlm(1)))

      call nwpw_timing_end(38)
      return
      end


*     *************************************************
*     *                                               *
*     *        nwpw_compcharge_add_dElocal_Qlm        *
*     *                                               *
*     *************************************************
*
      subroutine nwpw_compcharge_add_dElocal_Qlm(ispin)
      implicit none
      integer ispin

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      integer k
c      call daxpy(2*nion_paw*lm_size_max,1.0d0,
c     >           dbl_mb(dElocal_Qlm(1)),1,
c     >           dbl_mb(dE_Qlm(1)),1)
!$OMP DO
      do k=1,2*nion_paw*lm_size_max
         dbl_mb(dE_Qlm(1)+k-1) = dbl_mb(dE_Qlm(1)+k-1) 
     >                         + dbl_mb(dElocal_Qlm(1)+k-1) 
      end do
!$OMP END DO

      return
      end

*     *************************************************
*     *                                               *
*     *        nwpw_compcharge_gen_dEmult_Qlm         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_dEmult_Qlm(ispin)
      implicit none
      integer ispin

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer ms,k,indx1up,indx1dn,indx2up,indx2dn,indx1,indx2
      real*8  q1,q2
 
      call nwpw_timing_start(37)

      !call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dEmult_Qlm(1)),1)
      call Parallel_shared_vector_zero(.true.,2*nion_paw*lm_size_max,
     >                                 dbl_mb(dEmult_Qlm(1)))
!$OMP MASTER
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx1dn = indx1up + (ispin-1)*lm_size_max
         indx2up = int_mb(lm2_gauss(1)+k-1)
         indx2dn = indx2up + (ispin-1)*lm_size_max
         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)

         do ms=1,ispin
            indx1 = indx1up + (ms-1)*lm_size_max
            indx2 = indx2up + (ms-1)*lm_size_max
            dbl_mb(dEmult_Qlm(1)+indx1) = dbl_mb(dEmult_Qlm(1)+indx1)
     >                                  + q2*dbl_mb(e_gauss(1)+k-1)
            dbl_mb(dEmult_Qlm(1)+indx2) = dbl_mb(dEmult_Qlm(1)+indx2)
     >                                  + q1*dbl_mb(e_gauss(1)+k-1)
         end do

      end do
!$OMP END MASTER
!$OMP BARRIER
      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dEmult_Qlm(1)))

      call nwpw_timing_end(37)
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 q1,q2

c      call daxpy(2*nion_paw*lm_size_max,1.0d0,
c     >           dbl_mb(dEmult_Qlm(1)),1,
c     >           dbl_mb(dE_Qlm(1)),1)
!$OMP DO
      do k=1,2*nion_paw*lm_size_max
         dbl_mb(dE_Qlm(1)+k-1) = dbl_mb(dE_Qlm(1)+k-1)
     >                         + dbl_mb(dEmult_Qlm(1)+k-1)
      end do
!$OMP END DO

cc      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,dbl_mb(dEmult_Qlm(1)),1)
c      do k=1,ngauss
c         indx1up = int_mb(lm1_gauss(1)+k-1)
c         indx1dn = indx1up + (ispin-1)*lm_size_max
c         indx2up = int_mb(lm2_gauss(1)+k-1)
c         indx2dn = indx2up + (ispin-1)*lm_size_max
c         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
c         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)
c
c         dbl_mb(dE_Qlm(1)+indx1up) = dbl_mb(dE_Qlm(1)+indx1up) 
c     >                                 + q2*dbl_mb(e_gauss(1)+k-1)
c         dbl_mb(dE_Qlm(1)+indx1dn) = dbl_mb(dE_Qlm(1)+indx1dn) 
c     >                                 + q2*dbl_mb(e_gauss(1)+k-1)
c
c         dbl_mb(dE_Qlm(1)+indx2up) = dbl_mb(dE_Qlm(1)+indx2up) 
c     >                                 + q1*dbl_mb(e_gauss(1)+k-1)
c         dbl_mb(dE_Qlm(1)+indx2dn) = dbl_mb(dE_Qlm(1)+indx2dn) 
c     >                                 + q1*dbl_mb(e_gauss(1)+k-1)
c
c      end do
cc      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
cc     >                            dbl_mb(dE_Qlm(1)))
      
      return
      end

*     *************************************************
*     *                                               *
*     *        nwpw_compcharge_E_multipole_zv         *
*     *                                               *
*     *************************************************

      real*8 function nwpw_compcharge_E_multipole_zv(ispin,zv)
      implicit none
      integer ispin
      real*8 zv(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      real*8 E,q1,q2,q11,q22,fourpi
      integer iii,jjj,l1,m1,l2,m2,ia,iia,ja,jja

      fourpi = 16.0d0*datan(1.0d0)

      E = 0.0d0
!$OMP MASTER
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
            q22 = -zv(ja)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if


         q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)+q11
         q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)+q22
         E = E + q1*q2*dbl_mb(e_gauss(1)+k-1)

      end do
!$OMP END MASTER
!$OMP BARRIER
      call Parallel_SumAll(E)

      nwpw_compcharge_E_multipole_zv = E
      return
      end


*     *************************************************
*     *                                               *
*     *        nwpw_compcharge_E_multipole_zv_zv      *
*     *                                               *
*     *************************************************

      real*8 function nwpw_compcharge_E_multipole_zv_zv(ispin,zv)
      implicit none
      integer ispin
      real*8 zv(*)

#include "bafdecls.fh"
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
            q22 = -zv(ja)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if

         E = E + q11*q22*dbl_mb(e_gauss(1)+k-1)

      end do
      call Parallel_SumAll(E)

      nwpw_compcharge_E_multipole_zv_zv = E
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iia,shift

      call nwpw_timing_start(35)

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


      call nwpw_timing_end(35)
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
      real*8  coeff,w,w1,scal,eh0

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      scal = 1.0d0/lattice_omega()
      eh0 = 0.0d0

      do i=1,nindx
         iprj  = iprj_hartree(i)
         jprj  = jprj_hartree(i)
         iprj1 = iprj1_hartree(i)
         jprj1 = jprj1_hartree(i)
         coeff = coeff_hartree(i)
         w  = wmatrix(iprj,jprj,1)   + wmatrix(iprj,jprj,ispin)
         w1 = wmatrix(iprj1,jprj1,1) + wmatrix(iprj1,jprj1,ispin)
         eh0 = eh0 + coeff*w*w1*scal*scal*0.5d0
!$OMP DO
         do n=1,ne(1)+ne(2)
            sw2(n,iprj) = sw2(n,iprj) 
     >                  + scal*coeff*sw1(n,jprj)*w1*0.5d0
            sw2(n,iprj1) = sw2(n,iprj1) 
     >                   + scal*coeff*sw1(n,jprj1)*w*0.5d0
         end do
!$OMP END DO
      end do

!$OMP MASTER
      eh = eh0
!$OMP END MASTER
   
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

#include "bafdecls.fh"
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
      common /nwpw_compcharge_coulomb_e_sub_tmp/ eh

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      scal = 1.0d0/lattice_omega()
!$OMP MASTER
      eh = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:eh)
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
!$OMP END DO

      nwpw_compcharge_coulomb_e_sub = eh
      return
      end



      subroutine nwpw_compcharge_indxiiilm(indx,iii,l,m)
      implicit none
      integer indx
      integer iii,l,m

#include "bafdecls.fh"
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

      subroutine nwpw_compcharge_iiilmindx(iii,l,m,indx)
      implicit none
      integer iii,l,m
      integer indx

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      indx = (iii-1)*(2*lm_size_max) +  l*(l+1) + m
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c debug cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



*     *************************************************
*     *                                               *
*     *   nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_test   *
*     *                                               *
*     *************************************************
*
      subroutine nwpw_compcharge_gen_dE_ncmp_vloc_Qlm_test(ispin,zv,
     >                                           vcmp,vcmp_smooth,
     >                                           move,fion)
      implicit none
      integer ispin
      real*8     zv(*)
      complex*16 vcmp(*)
      complex*16 vcmp_smooth(*)
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     *** local variables ****
      logical value
      integer np,np_j,np_k,taskid,taskid_j,taskid_k,pcount
      integer iii,iia,ii,ia,lm,l,m,k,jjj,ja,jja,l1,m1,l2,m2,ms
      integer indx1up,indx2up,indx1,indx2,Gx,Gy,Gz
      integer glm_cmp(2),glm_cmp_smooth(2),gvtmp(2),ftmp(2)
      real*8  q11,q22,e1,e2,fourpi
      real*8   qq,e1x,e1y,e1z,e2x,e2y,e2z

      integer  Pack_G_indx,c_G_indx
      external Pack_G_indx,c_G_indx

      write(*,*) "I am here"
      call nwpw_timing_start(38)

      call dcopy(2*nion_paw*lm_size_max,0.0d0,0,
     >           dbl_mb(dElocal_Qlm(1)),1)

*     **** allocating stack memory ****
      value = BA_push_get(mt_dcpl,npack0,'glm_cmp',
     >                    glm_cmp(2),glm_cmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack0,'glm_cmp_smooth',
     >                    glm_cmp_smooth(2),glm_cmp_smooth(1))
      if (move) then
         value = value.and.
     >           BA_push_get(mt_dbl,npack0,'gvtmp',
     >                       gvtmp(2),gvtmp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,3*nion_paw,'ftmp',
     >                       ftmp(2),ftmp(1))
         call dcopy(3*nion_paw,0.0d0,0,dbl_mb(ftmp(1)),1)
      end if
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:out of stack',
     >               0,MA_ERR)

      if (isgamma) then
      if (move) then
         Gx = Pack_G_indx(0,1)
         Gy = Pack_G_indx(0,2)
         Gz = Pack_G_indx(0,3)
      end if
     
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)
         if (move) call dcopy(npack0,0.0d0,0,dbl_mb(gvtmp(1)),1)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np_j).eq.taskid_j) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           e1)
               call Pack_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           e2)
               !lm = l*(l+1)+m
               do ms=1,ispin
                  indx1  = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx1) = (e1+e2)
               end do
               if (move) then
                  indx1  = (iii-1)*2*lm_size_max+lm
                  indx2  = (iii-1)*2*lm_size_max
     >                   + (ispin-1)*lm_size_max+lm
                  qq = dbl_mb(Qlm(1)+indx1) + dbl_mb(Qlm(1)+indx2)
                  call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp(1)),
     >                                       vcmp,
     >                                       dbl_mb(gvtmp(1)))
                  call Pack_cct_iaconjgMulAdd(0,qq,
     >                                       dcpl_mb(glm_cmp_smooth(1)),
     >                                       vcmp_smooth,
     >                                       dbl_mb(gvtmp(1)))
               end if

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
         if (move) then
            call Pack_tt_idot(0,dbl_mb(Gx),dbl_mb(gvtmp(1)),e1x)
            call Pack_tt_idot(0,dbl_mb(Gy),dbl_mb(gvtmp(1)),e1y)
            call Pack_tt_idot(0,dbl_mb(Gz),dbl_mb(gvtmp(1)),e1z)
            dbl_mb(ftmp(1)+3*(iii-1))  = dbl_mb(ftmp(1)+3*(iii-1))  -e1x
            dbl_mb(ftmp(1)+3*(iii-1)+1)= dbl_mb(ftmp(1)+3*(iii-1)+1)-e1y
            dbl_mb(ftmp(1)+3*(iii-1)+2)= dbl_mb(ftmp(1)+3*(iii-1)+2)-e1z
         end if
      end do !** iii **


      else

      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)
      taskid = taskid_j + np_j*taskid_k
      np     = np_j+np_k
      pcount = 0
      do iii=1,nion_paw
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         iia = int_mb(katm_paw(1)+iii-1)
         ia  = int_mb(katm_pawtoion(1)+iia-1)

         lm = 0
         do l=0,int_mb(mult_l(1)+iia-1)
         do m=-l,l
            if (mod(pcount,np).eq.taskid) then
               call nwpw_compcharge_gen_glm2(ii,l,m,
     >                              dcpl_mb(glm_cmp(1)),
     >                              dcpl_mb(glm_cmp_smooth(1)))
               call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp(1)),
     >                           vcmp,
     >                           e1)
               call Cram_cc_idot(0,
     >                           dcpl_mb(glm_cmp_smooth(1)),
     >                           vcmp_smooth,
     >                           e2)
               !lm = l*(l+1)+m
               do ms=1,ispin
                  indx1 = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max+lm
                  dbl_mb(dElocal_Qlm(1)+indx1) = (e1+e2)
               end do

            end if
            lm     = lm + 1
            pcount = pcount + 1
         end do
         end do
      end do !** iii **

      end if

      if (move) then
         call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(ftmp(1)))
         do iii=1,nion_paw
            ii  = int_mb(ion_pawtoion(1)+iii-1)
            fion(1,ii) = fion(1,ii) + dbl_mb(ftmp(1)+3*(iii-1))
            fion(2,ii) = fion(2,ii) + dbl_mb(ftmp(1)+3*(iii-1)+1)
            fion(3,ii) = fion(3,ii) + dbl_mb(ftmp(1)+3*(iii-1)+2)
         end do
      end if

*     **** popping stack memory ****
      value = .true.
      if (move) then
         value = value.and.BA_pop_stack(ftmp(2))
         value = value.and.BA_pop_stack(gvtmp(2))
      end if
      value = value.and.BA_pop_stack(glm_cmp_smooth(2))
      value = value.and.BA_pop_stack(glm_cmp(2))
      if (.not.value)
     >  call errquit('nwpw_compcharge_gen_Elocal_Qlm:popping stack',
     >               1,MA_ERR)


      fourpi = 16.0d0*datan(1.0d0)
      do k=1,ngauss
         indx1up = int_mb(lm1_gauss(1)+k-1)
         indx2up = int_mb(lm2_gauss(1)+k-1)

         call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
         call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)

         if (l1.eq.0) then
            iia = int_mb(katm_paw(1)+iii-1)
            ia  = int_mb(katm_pawtoion(1)+iia-1)
            q11 = -1.0d0*zv(ia)/dsqrt(fourpi)
         else
            q11 = 0.0d0
         end if

         if (l2.eq.0) then
            jja = int_mb(katm_paw(1)+jjj-1)
            ja  = int_mb(katm_pawtoion(1)+jja-1)
            q22 = -1.0d0*zv(ja)/dsqrt(fourpi)
         else
            q22 = 0.0d0
         end if

         do ms=1,ispin
            indx1 = indx1up + (ms-1)*lm_size_max
            indx2 = indx2up + (ms-1)*lm_size_max
            dbl_mb(dElocal_Qlm(1)+indx1) = dbl_mb(dElocal_Qlm(1)+indx1)
     >                                   + q22*dbl_mb(e_gauss(1)+k-1)
            dbl_mb(dElocal_Qlm(1)+indx2) = dbl_mb(dElocal_Qlm(1)+indx2)
     >                                   + q11*dbl_mb(e_gauss(1)+k-1)
         end do
      end do
      call Parallel_Vector_SumAll(2*nion_paw*lm_size_max,
     >                            dbl_mb(dElocal_Qlm(1)))

      call nwpw_timing_end(38)
      return
      end

*     *************************************************
*     *                                               *
*     *           nwpw_compcharge_F_multipole_test    *
*     *                                               *
*     *************************************************
*
*  This routine computes the multiple only contribution to the force.
*  Note that the components of the force due to the variation wrt ncmp
*  are included in the standard non-local psp forces.
*
      subroutine nwpw_compcharge_F_multipole_test(ispin,zv,fion)
      implicit none
      integer ispin
      real*8 zv(*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer k,indx1up,indx1dn,indx2up,indx2dn
      integer iii1,iii2,ii1,ii2
      integer iii,l1,m1,iia,ia
      integer jjj,l2,m2,jja,ja
      integer fmult(2)
      real*8 q1,q1x,q1y,q1z,q2,q2x,q2y,q2z
      real*8 q11,q22,fourpi
      real*8 tmp,tmpx,tmpy,tmpz,R,R12(3)

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      fourpi = 16.0d0*datan(1.0d0)

      if (.not.BA_push_get(mt_dbl,3*nion_paw,"fmult",fmult(2),fmult(1)))
     >   call errquit("nwpw_compcharge_F_multipole:stack",1,MA_ERR)
      call dcopy(3*nion_paw,0.0d0,0,dbl_mb(fmult(1)),1)

      do k=1,ngauss
         iii1 = int_mb(iii1_gauss(1)+k-1)
         iii2 = int_mb(iii2_gauss(1)+k-1)
         if (iii1.ne.iii2) then
            indx1up = int_mb(lm1_gauss(1)+k-1)
            indx1dn = indx1up + (ispin-1)*lm_size_max
            indx2up = int_mb(lm2_gauss(1)+k-1)
            indx2dn = indx2up + (ispin-1)*lm_size_max

            call nwpw_compcharge_indxiiilm(indx1up,iii,l1,m1)
            call nwpw_compcharge_indxiiilm(indx2up,jjj,l2,m2)

            if (l1.eq.0) then
               iia = int_mb(katm_paw(1)+iii-1)
               ia  = int_mb(katm_pawtoion(1)+iia-1)
               q11 = -1.0d0*zv(ia)/dsqrt(fourpi)
            else
               q11 = 0.0d0
            end if
            if (l2.eq.0) then
               jja = int_mb(katm_paw(1)+jjj-1)
               ja  = int_mb(katm_pawtoion(1)+jja-1)
               q22 = -1.0d0*zv(ja)/dsqrt(fourpi)
            else
               q22 = 0.0d0
            end if

            q1 = dbl_mb(Qlm(1)+indx1up) + dbl_mb(Qlm(1)+indx1dn)
            q2 = dbl_mb(Qlm(1)+indx2up) + dbl_mb(Qlm(1)+indx2dn)

            !tmp  = (q1*q2 + q11*q2 + q1*q22)
            tmp  = (q11*q2 + q1*q22)
            !tmp  = q1*q2

            call daxpy(3,-tmp,dbl_mb(f_gauss(1)+3*(k-1)), 1,
     >                        dbl_mb(fmult(1)+3*(iii1-1)),1)
            call daxpy(3, tmp,dbl_mb(f_gauss(1)+3*(k-1)), 1,
     >                        dbl_mb(fmult(1)+3*(iii2-1)),1)

         end if
      end do
      call Parallel_Vector_SumAll(3*nion_paw,dbl_mb(fmult(1)))

      do iii1=1,nion_paw
         ii1 = int_mb(ion_pawtoion(1)+iii1-1)
         fion(1,ii1) = fion(1,ii1) + dbl_mb(fmult(1)+3*(iii1-1))
         fion(2,ii1) = fion(2,ii1) + dbl_mb(fmult(1)+3*(iii1-1)+1)
         fion(3,ii1) = fion(3,ii1) + dbl_mb(fmult(1)+3*(iii1-1)+2)
      end do

      if (.not.BA_pop_stack(fmult(2)))
     >  call errquit("nwpw_compcharge_F_multipole_test:stack",2,MA_ERR)

      return
      end


