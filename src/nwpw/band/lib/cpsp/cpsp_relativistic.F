*
* $Id$
*


*     *******************************************
*     *                                         *
*     *             Multiply_Kijl_SO            *
*     *                                         *
*     *******************************************
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Create the G(i,j,L)*<psi|(LS)|prj(j,L)> array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Multiply_Kijl_SO(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             zsw1,zsw2)
      implicit none
      integer nn
      integer nprj,nmax,lmax,nh
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax),xx,dl2,dm2,dlm
      complex*16 zsw1(nn,nprj)
      complex*16 zsw2(nn,nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb

      nh=nn/2
      call dcopy(2*nn*nprj,0.0d0,0,zsw2,1)
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)

         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)

            if ((la.eq.lb).and.(ma.eq.mb)) then
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dble(ma)  
              call daxpy(nn,xx,zsw1(1,a),1,zsw2(1,b),1)
              call daxpy(nn,(-xx),zsw1(nh+1,a),1,zsw2(nh+1,b),1)
            end if

         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb+1)).and.(mb.ne.lb)) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1(1,a),1,zsw2(nh+1,b),1)
            end if
         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb-1)).and.(mb.ne.(-lb))) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1(nh+1,a),1,zsw2(1,b),1)
            end if
         end do
      end do
      return
      end


*     *******************************************
*     *                                         *
*     *             Multiply_Kijl_SO_x          *
*     *                                         *
*     *******************************************
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Create the G(i,j,L)*<psi|(LS)|prj(j,L)> array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Multiply_Kijl_SO_x(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             zsw1u,zsw1d,zsw2u,zsw2d)
      implicit none
      integer nn
      integer nprj,nmax,lmax
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax),xx,dl2,dm2,dlm
      complex*16 zsw1u(nprj),zsw1d(nprj)
      complex*16 zsw2u(nprj),zsw2d(nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb

      call dcopy(2*nn*nprj,0.0d0,0,zsw2u,1)
      call dcopy(2*nn*nprj,0.0d0,0,zsw2d,1)
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)

         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)

            if ((la.eq.lb).and.(ma.eq.mb)) then
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dble(ma)  
              call daxpy(nn,xx,zsw1u(a),1,zsw2u(b),1)
              call daxpy(nn,(-xx),zsw1d(a),1,zsw2d(b),1)
            end if

         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb+1)).and.(mb.ne.lb)) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1u(a),1,zsw2d(b),1)
            end if
         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb-1)).and.(mb.ne.(-lb))) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1d(a),1,zsw2u(b),1)
            end if
         end do
      end do
      return
      end

*     *******************************************
*     *				  		*
*     *	 	 cpsp_v_nonlocal_orb_2com  	*
*     *						*
*     *******************************************

      subroutine cpsp_v_nonlocal_orb_2com(nb,orb1,orb2)
      implicit none
      integer    nb
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone
c      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,npack1,npack,nbrill
      integer ii,ia,l,prj_shift
      integer shift,l_prj,nproj,shifts,ne1
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1u(2),zsw2u(2),zsw1d(2),zsw2d(2)
      logical value,sd_function
 
*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,brillioun_nbrillioun
      integer  cpsi_ne,cpsp_projector_get_ptr
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,brillioun_nbrillioun
      external lattice_omega,cpsi_ne,cpsp_projector_get_ptr

      one=dcmplx(1.0d0,0.0d0)
      mone=dcmplx(-1.0d0,0.0d0)
      
      call nwpw_timing_start(6)
         
      prj_shift=vso_shift
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nbrill = brillioun_nbrillioun()

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw1u',zsw1u(2),zsw1u(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw2u',zsw2u(2),zsw2u(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw1d',zsw1d(2),zsw1d(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw2d',zsw2d(2),zsw2d(1))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb2com:pushing stack',0,MA_ERR)

      ne1=cpsi_ne(1)
      shifts = npack1*ne1 + 1
      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
         ia=ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
c            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
            call zscal(npack,cxr,dcpl_mb(exi(1)),1)

            do l=1,nproj

c              shift = vnl(1)+(l-1) *npack1
c     >                      +(nb-1)*npack1*vnl_stride
c     >                      +(ia-1)*npack1*vnl_stride*nbrill
              shift = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nb,l)

              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*jmmax_max)

              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if
*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****
*             *** current function is s or d ****
              if (sd_function) then
                 call Cram_rc_Mul(nb,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
              else
                 call Cram_irc_Mul(nb,dbl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
              end if

*             **** compute 1Xnproj matrix zsw1 = <orb1|prj> ****
              call Cram_cc_izdot(nb,
     >                      orb1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1u(1)+(l-1)))
              call Cram_cc_izdot(nb,
     >                      orb1(shifts),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1d(1)+(l-1)))
            end do !**l**
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1u(1)))
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1d(1)))

*           **** zsw2 = Gijl*zsw1 ******
            call Multiply_Gijl_zsw1(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb( Gijl(1)
     >                         + (ia-1)*gij_stride),
     >                         dcpl_mb(zsw1u(1)),
     >                         dcpl_mb(zsw2u(1)))
            call Multiply_Gijl_zsw1(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*gij_stride),
     >                         dcpl_mb(zsw1d(1)),
     >                         dcpl_mb(zsw2d(1)))

*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2u(1)),1)
            call dscal(2*nproj,scal,dcpl_mb(zsw2d(1)),1)
            call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2u(1)),   1,
     >                 one,
     >                 orb2, npack1)
            call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2d(1)),   1,
     >                 one,
     >                 orb2(shifts), npack1)

         end if !** nproj>0 **
      end do !** ii***

      value =           BA_pop_stack(zsw2d(2))
      value = value.and.BA_pop_stack(zsw1d(2))
      value = value.and.BA_pop_stack(zsw2u(2))
      value = value.and.BA_pop_stack(zsw1u(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb:popping stack',3,MA_ERR)

      call nwpw_timing_end(6)

      return 
      end



*     ***********************************
*     *					*
*     *	 	 cpsp_v_spin_orbit  	*
*     *					*
*     ***********************************

      subroutine cpsp_v_spin_orbit(ispin,ne,
     >                           psi1_tag,psi2_tag,move,fion)
      implicit none
      integer    ispin,ne(2)
      integer    psi1_tag
      integer    psi2_tag
c      complex*16 psi1(*),psi2(*)
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     *** local variables ***
      complex*16 one,mone
      integer nfft3d,G(3),npack1,npack
      integer ii,ia,l,n,nn,nb,nbq,nbrill
      integer shift,l_prj,nproj
      integer psi1_shift,psi2_shift,psi_shift,nshift
      real*8  omega,weight,scal
      complex*16 cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx,Pneb_nbrillq,Pneb_convert_nb
      integer  cpsp_projector_get_ptr,cpsi_data_get_chnk
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,c_G_indx,Pneb_nbrillq,Pneb_convert_nb
      external cpsp_projector_get_ptr,cpsi_data_get_chnk
      external lattice_omega,brillioun_weight


      if (.not.do_spin_orbit) return

      one  = dcmplx( 1.0d0,0.0d0)
      mone = dcmplx(-1.0d0,0.0d0)

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      nbrill = Pneb_nbrillq()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nshift = 2*npack1

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     &       MA_ERR)

      if (move) then
        value = BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        value = value.and.
     >        BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
        if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      end if

      omega = lattice_omega()
      scal  = 1.0d0/omega
      do 200 ii=1,ion_nion()
        ia=ion_katm(ii)
cccccccccc if this atom is not HGH PPOT skip it....
        if (int_mb(psp_type(1)+ia-1).ne.1) goto 200
cccccccccccccccccccccccccccccccccccccccccccc
        nproj = int_mb(nprj(1)+ia-1)
        if (nproj.gt.0) then
        do nbq=1,nbrill
           nb = Pneb_convert_nb(nbq)
           psi1_shift = cpsi_data_get_chnk(psi1_tag,nbq)
           psi2_shift = cpsi_data_get_chnk(psi2_tag,nbq)
           call Cram_npack(nb,npack)

*       **** structure factor pseudopotential ****
           call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
           call cstrfac_k(ii,nb,cxr)
c           call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
           call zscal(npack,cxr,dcpl_mb(exi(1)),1)


*          **** generate zsw1's and projectors ****
           do 105 l=1,nproj

c              shift = vnlso(1)+(l-1)*npack1
c     >                      +(nb-1)*npack1*vso_stride
c     >                      +(ia-1)*npack1*vso_stride*nbrill
              shift = cpsp_projector_get_ptr(
     >                     int_mb(vnlso(1)+ia-1),nb,l)

              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*jmmax_max)
              if (l_prj.eq.0) then
                call dcopy(npack1*2,0.0d0,0,
     >               dcpl_mb(prjtmp(1)+(l-1)*npack1),1)
                goto 105
              end if

              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****

c             *** current function is s or d ****
                if (sd_function) then
                  call Cram_cc_Mul(nb,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
                else
                  call Cram_icc_Mul(nb,dbl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
                end if

*             **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
                call Cram_cc_inzdot(nb,nn,
     >                      dbl_mb(psi1_shift),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1(1)+(l-1)*nn))
           
 105       continue 
           call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))
         

*          **** zsw2 = Kijl*zsw1 ******
           call Multiply_Kijl_SO(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb(Kijl(1)+(ia-1)*kij_stride),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))

*          **** do Kleinman-Bylander Multiplication ****
           call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
           call ZGEMM('N','C',npack,nn,nproj,
     >                mone,
     >                dcpl_mb(prjtmp(1)), npack1,
     >                dcpl_mb(zsw2(1)),   nn,
     >                one,
     >                dbl_mb(psi2_shift), npack1)

           if (move) then
              weight = brillioun_weight(nb)
              if (ispin.eq.1) 
     >           call dscal(2*nn*nproj,2.0d0,dcpl_mb(zsw2(1)),1)
          
              do l=1,nproj
          
                 psi_shift = psi1_shift
                 do n=1,nn

                    call Cram_zccr_Multiply2(nb,
     >                                 dcpl_mb(zsw2(1)+(l-1)*nn+n-1),
     >                                 dbl_mb(psi_shift),
     >                                 dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                                 dbl_mb(xtmp(1)))
                    psi_shift = psi_shift + nshift

c                    do i=1,npack
c                    ctmp = psi1(i+(n-1)*npack1
c     >                           +(nb-1)*neall*npack1)
c     >                   *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
c     >                   *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
c                    dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
c                    end do

*                   **** define Gx,Gy and Gz in packed space ****
                    call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
                    call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
                    call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
                    call Cram_r_pack(nb,dbl_mb(Gx(1)))
                    call Cram_r_pack(nb,dbl_mb(Gy(1)))
                    call Cram_r_pack(nb,dbl_mb(Gz(1)))
                    call Cram_rr_idot(nb,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do !**n**

                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
                 do n=1,nn
                    fion(1,ii) = fion(1,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+3*(n-1))
                    fion(2,ii) = fion(2,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+1+3*(n-1))
                    fion(3,ii) = fion(3,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+2+3*(n-1))
                 end do !** nn **

              end do !** l **
           end if !** move **
       
        end do !** nb **
        end if !** nproj>0**

200   continue !**ii**

      if (move) then
        value = BA_pop_stack(sum(2))
        value = value.and.BA_pop_stack(Gz(2))
        value = value.and.BA_pop_stack(Gy(2))
        value = value.and.BA_pop_stack(Gx(2))
        value = value.and.BA_pop_stack(xtmp(2))
        if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',2,
     &       MA_ERR)
      end if

      value =           BA_pop_stack(zsw2(2))
      value = value.and.BA_pop_stack(zsw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)

      return 
      end


*     *******************************************
*     *				  		*
*     *	 	 cpsp_v_spin_orbit_orb          *
*     *						*
*     *******************************************

      subroutine cpsp_v_spin_orbit_orb(nb,orb1,orb2)
      implicit none
      integer    nb
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"

*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,npack1,npack
      integer ii,ia,l
      integer shift,l_prj,nproj,shifts,ne1
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1u(2),zsw2u(2),zsw1d(2),zsw2d(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm
      integer  cpsi_ne,cpsp_projector_get_ptr
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm
      external lattice_omega,cpsi_ne,cpsp_projector_get_ptr

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw1u',zsw1u(2),zsw1u(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw2u',zsw2u(2),zsw2u(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw1d',zsw1d(2),zsw1d(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw2d',zsw2d(2),zsw2d(1))
      if (.not.value) 
     > call errquit('cpsp_v_spin_orbit_orb2com:pushing stack',0,MA_ERR)

      ne1=cpsi_ne(1)
      shifts = npack1*ne1 + 1
      omega = lattice_omega()
      scal  = 1.0d0/omega

      do 200 ii=1,ion_nion()
         ia=ion_katm(ii)

         if (int_mb(psp_type(1)+ia-1).ne.1) goto 200

         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
c            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
            call zscal(npack,cxr,dcpl_mb(exi(1)),1)

            do l=1,nproj

c              shift = vnlso(1)+(l-1) *npack1
c     >                      +(nb-1)*npack1*vso_stride
c     >                      +(ia-1)*npack1*vso_stride*nbrill
              shift = cpsp_projector_get_ptr(
     >                     int_mb(vnlso(1)+ia-1),nb,l)

              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*jmmax_max)

              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****
*             *** current function is s or d ****
              if (sd_function) then
                 call Cram_cc_Mul(nb,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
              else
                 call Cram_icc_Mul(nb,dbl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
              end if

*             **** compute 1Xnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_izdot(nb,
     >                      orb1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1u(1)+(l-1)))
              call Cram_cc_izdot(nb,
     >                      orb1(shifts),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1d(1)+(l-1)))
            end do !**l**
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1u(1)))
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1d(1)))

*           **** zsw2 = Gijl*zsw1 ******
            call Multiply_Kijl_SO_x(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb(Kijl(1)
     >                         + (ia-1)*kij_stride),
     >                         dcpl_mb(zsw1u(1)),dcpl_mb(zsw1d(1)),
     >                         dcpl_mb(zsw2u(1)),dcpl_mb(zsw2d(1)))
*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2u(1)),1)
            call dscal(2*nproj,scal,dcpl_mb(zsw2d(1)),1)
            call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2u(1)),   1,
     >                 one,
     >                 orb2, npack1)
            call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2d(1)),   1,
     >                 one,
     >                 orb2(shifts), npack1)

         end if !** nproj>0 **
200   continue !** ii***

      value =           BA_pop_stack(zsw2d(2))
      value = value.and.BA_pop_stack(zsw1d(2))
      value = value.and.BA_pop_stack(zsw2u(2))
      value = value.and.BA_pop_stack(zsw1u(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('cpsp_v_spin_orbit_orb:popping stack',3,MA_ERR)

      call nwpw_timing_end(6)

      return 
      end

*     ***********************************************
*     *                                             *
*     *             cpsp_v_nonlocal_rel             *
*     *                                             *
*     ***********************************************

c Dirac Two Component Relativistic 
c   Non-Local Pseudopotential
c    pjn --- use at your own peril

      subroutine cpsp_v_nonlocal_rel(nb,ii,ia,nproj,npack1,nbrill,nn,ne,
     >                               exi,prtmp,zsw1,psi1,psi2,
     >                               Gij,bz_vol)
      implicit none
      integer    nb,nn,npack1,nproj,ia,ii,nbrill,ne(2)  
      complex*16 prtmp(*)
      complex*16 zsw1(*)
      complex*16 exi(*)
      complex*16 psi1(*)
      complex*16 psi2(*) 
      real*8     bz_vol,Gij(*)

#include "cpsp_common.fh"
#include "bafdecls.fh"

*     **** local variables ****
      complex*16 one,mone,ione
      parameter(one=(1.0d0,0.0d0),mone=(-1.0d0,0.0d0))
      parameter(ione=(0.0d0,1.0d0))

      complex*16 cxr
      real*8 xx
      integer l,shifts,shiftxu,shiftxd,l_prj,npack,ne1
      integer zshft,vshiftu,vshiftd
      logical sd_function,val

*     **** external functions ****
      integer  cpsp_projector_get_ptr
      external cpsp_projector_get_ptr

cccccccccccccccccccccccccccccccccccccccccc
c These should be declared as parameters
c however there seems to be no "standard"
c way to declare a complex parameter that
c works with all compilers...gfortran,ifort,
c pf90,etc.
ccccccccccccccccccccccccccccccccccccccccc
c      one  = dcmplx( 1.0d0,0.0d0)
c      mone = dcmplx(-1.0d0,0.0d0)
c      ione = dcmplx( 0.0d0,1.0d0)

      ne1=ne(1)
      shifts=1+npack1*ne1
      call Cram_npack(nb,npack)

*       **** structure factor pseudopotential ****
      call cstrfac_pack(nb,ii,exi)
      call cstrfac_k(ii,nb,cxr)
      call zscal(npack,cxr,exi,1)


*     **** generate zsw1's and projectors ****
      do 105 l=1,nproj

              vshiftu = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nb,l)
              vshiftd = cpsp_projector_get_ptr(
     >                     int_mb(vnlso(1)+ia-1),nb,l)
              shiftxu=1 + (l-1)*npack1 
              shiftxd=shiftxu + nproj*npack1
              zshft=1+(l-1)*ne1
              l_prj = int_mb(l_projector(1)+(l-1)
     >                          +(ia-1)*jmmax_max)
              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if
*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****

c             *** current function is s or d ****
              if (sd_function) then
                call Cram_cc_Mul(nb,dbl_mb(vshiftu),
     >                              exi,
     >                              prtmp(shiftxu))
                call Cram_cc_Mul(nb,dbl_mb(vshiftd),
     >                              exi,
     >                              prtmp(shiftxd))
*             *** current function is p or f ****
              else
                call Cram_icc_Mul(nb,dbl_mb(vshiftu),
     >                               exi,
     >                               prtmp(shiftxu))

                call Cram_icc_Mul(nb,dbl_mb(vshiftd),
     >                              exi,
     >                              prtmp(shiftxd))
              end if
         
*             **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_inzdot(nb,ne1,
     >                            psi1,
     >                            prtmp(shiftxu),
     >                            zsw1(zshft))
              call Cram_cc_inzdotAdd(nb,ne1,
     >                            psi1(shifts),
     >                            prtmp(shiftxd),
     >                            zsw1(zshft))
              xx=Gij(l)*bz_vol
              call dscal(2*ne1,xx,zsw1(zshft),1)
 105  continue 


      call C3dB_Vector_SumAll((2*ne1*nproj),zsw1)

      call ZGEMM('N','C',npack,ne1,nproj,
     >                mone,
     >                prtmp, npack1,
     >                zsw1, ne1,
     >                one,
     >                psi2(1), npack1)
      call ZGEMM('N','C',npack,ne1,nproj,
     >                mone,
     >                prtmp(1+npack1*nproj), npack1,
     >                zsw1, ne1,
     >                one,
     >                psi2(shifts), npack1)

      return
      end


ccccccccccc
*     ***********************************
*     *					*
*     *	 	 cpsp_v_nonlocal_orb  	*
*     *					*
*     ***********************************
      subroutine cpsp_v_nonlocal_rel_orb(nb,orb1,orb2,
     >  zsw1,Gij,exi,nfft3d,ia,ii,ne1,npack1,nproj)
      implicit none
      integer    nb
      real*8 Gij(*)
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone,ione
      integer nfft3d,npack1,npack
      integer ii,ia,l,ne1
      integer l_prj,nproj
      integer shiftu,shiftd,shifts,shiftux,shiftdx
      real*8  omega,scal
      complex*16 cxr
      complex*16 exi(*)
      complex*16 zsw1(*)
      logical sd_function

*     **** external functions ****
      logical  is_sORd,cpsi_spin_orbit
      integer  ion_nion,ion_katm
      integer  cpsp_projector_get_ptr
      real*8   lattice_omega
      external is_sORd,cpsi_spin_orbit
      external ion_nion,ion_katm
      external cpsp_projector_get_ptr
      external lattice_omega

      one=dcmplx(1.0d0,0.0d0)
      mone=dcmplx(-1.0d0,0.d0)
      ione=dcmplx(0.0d0,1.d0)

      omega = lattice_omega()
      scal  = 1.0d0/omega
      shifts= ne1*npack1+1
*           **** structure factor ****
      
      call Cram_npack(nb,npack)
      call cstrfac_pack(nb,ii,exi)
      call cstrfac_k(ii,nb,cxr)
c      call Cram_c_ZMul(nb,cxr,exi,exi)
      call zscal(npack,cxr,exi,1)

      do l=1,nproj
              shiftu = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nb,l)
              shiftd = cpsp_projector_get_ptr(
     >                     int_mb(vnlso(1)+ia-1),nb,l)
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*jmmax_max)
              shiftux=(l-1)*npack1
              shiftdx=shiftux+npack1*nproj
              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****
*             *** current function is s or d ****
              if (sd_function) then
                 call Cram_cc_Mul(nb,dbl_mb(shiftu),exi,
     >                     dcpl_mb(prjtmp(1)+shiftux))
                 call Cram_cc_Mul(nb,dbl_mb(shiftd),exi,
     >                     dcpl_mb(prjtmp(1)+shiftdx))

*             *** current function is p or f ****
              else
                 call Cram_icc_Mul(nb,dbl_mb(shiftu),exi,
     >                     dcpl_mb(prjtmp(1)+shiftux))
                 call Cram_icc_Mul(nb,dbl_mb(shiftd),exi,
     >                     dcpl_mb(prjtmp(1)+shiftdx))
              end if

*             **** compute 1Xnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_izdot(nb,
     >                      orb1,
     >                      dcpl_mb(prjtmp(1)+shiftux),
     >                      zsw1(1+(l-1)))
              call Cram_cc_izdotAdd(nb,
     >                      orb1(shifts),
     >                      dcpl_mb(prjtmp(1)+shiftdx),
     >                      zsw1(1+(l-1)))
              zsw1(1+(l-1))=zsw1(1+(l-1))*Gij(l)*scal
      end do !**l**

      call C3dB_Vector_SumAll((2*nproj),zsw1)

*           **** do Kleinman-Bylander Multiplication ****
      call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 zsw1,   1,
     >                 one,
     >                 orb2, npack1)

      call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)+nproj*npack1), npack1,
     >                 zsw1,   1,
     >                 one,
     >                 orb2(shifts), npack1)

      return 
      end
ccccccccccccccccc
      subroutine cpsp_f_nonlocal_rel(nb,ii,ia,nproj,npack1,ne1,
     >  exi,zsw1,psi,prtmp,xtmp,sum1,Gij,Gx,Gy,Gz,
     >  fx,fy,fz,weight,scal)
      implicit none
      integer ia,ii,nproj,ne1,npack1,nb,npack
      complex*16 exi(*),zsw1(*),psi(*),prtmp(*),cxr
      real*8 xtmp(*),sum1(*),fx,fy,fz,xx
      real*8 weight,scal,Gij(*),Gx(*),Gy(*),Gz(*)
#include "cpsp_common.fh"
#include "bafdecls.fh"
ccccccccccc locals
      integer shiftxu,shiftxd,vshiftu,vshiftd,zshft
      integer l_prj,l,n,pshftu,pshftd,sshft,shifts
      logical sd_function
ccccccccccc external
      integer cpsp_projector_get_ptr
      external cpsp_projector_get_ptr
cccccccccccccccccccccccccccccccccccccccccccc
*             **** structure factor ****
      call Cram_npack(nb,npack)
      call cstrfac_pack(nb,ii,exi)
      call cstrfac_k(ii,nb,cxr)
c      call Cram_c_ZMul(nb,cxr,exi,exi)
      call zscal(npack,cxr,exi,1)

      shifts=1+npack1*ne1
      do l=1,nproj
         vshiftu = cpsp_projector_get_ptr(
     >       int_mb(vnl(1)+ia-1),nb,l)
         vshiftd = cpsp_projector_get_ptr(
     >       int_mb(vnlso(1)+ia-1),nb,l)
         shiftxu=1 + (l-1)*npack1 
         shiftxd=shiftxu + nproj*npack1
         
         l_prj = int_mb(l_projector(1)+(l-1)
     >      +(ia-1)*jmmax_max)

         sd_function=.true.
         if (mod(l_prj,2).ne.0) then
           sd_function=.false.
         end if

*       **** phase factor does not matter therefore ****
*       **** (-i)^l is the same as (i)^l in the     ****
*       **** Rayleigh scattering formula            ****
*       *** current function is s or d ****
         if (sd_function) then
           call Cram_cc_Mul(nb,dbl_mb(vshiftu),exi,prtmp(shiftxu))
           call Cram_cc_Mul(nb,dbl_mb(vshiftd),exi,prtmp(shiftxd))
         else
           call Cram_icc_Mul(nb,dbl_mb(vshiftu),exi,prtmp(shiftxu))
           call Cram_icc_Mul(nb,dbl_mb(vshiftd),exi,prtmp(shiftxd))
         end if

         zshft=1+(l-1)*ne1

         call Cram_cc_inzdot(nb,ne1,psi,prtmp(shiftxu),
     >      zsw1(zshft))
         call Cram_cc_inzdotAdd(nb,ne1,psi(shifts),prtmp(shiftxd),
     >      zsw1(zshft))

         xx=Gij(l)*scal
         call dscal(ne1*2,xx,zsw1(zshft),1)

      end do
      call C3db_Vector_SumAll(2*ne1*nproj,zsw1)
      do l=1,nproj
         shiftxu=1+(l-1)*npack1
         shiftxd=shiftxu+nproj*npack1
         do n=1,ne1
           zshft=(l-1)*ne1+n
           pshftu=1+(n-1)*npack1
           pshftd=pshftu+ne1*npack1
           call Cram_zccr_Multiply2(nb,
     >         zsw1(zshft),
     >         psi(pshftu),
     >         prtmp(shiftxu),
     >         xtmp)
           call Cram_zccr_Multiply2Add(nb,
     >         zsw1(zshft),
     >         psi(pshftd),
     >         prtmp(shiftxd),
     >         xtmp)
           
           sshft=1+3*(n-1)
           call Cram_rr_idot(nb,Gx,xtmp,sum1(sshft))
           call Cram_rr_idot(nb,Gy,xtmp,sum1(sshft+1))
           call Cram_rr_idot(nb,Gz,xtmp,sum1(sshft+2))
         end do
         call C3db_Vector_Sumall(3*ne1,sum1)
         do n=1,ne1
            sshft=1+3*(n-1)
            fx=fx+2.0d0*weight*sum1(sshft)
            fy=fy+2.0d0*weight*sum1(sshft+1)
            fz=fz+2.0d0*weight*sum1(sshft+2)
         end do
      end do
      return
      end
ccccccccccccccccccccccccc    
 


   


  


