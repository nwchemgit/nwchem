*
* $Id: paw_mult.F,v 1.8 2003-02-21 03:01:48 bylaska Exp $
*


!**************************************************
!
!       Name: paw_mult_init
!
!       Purpose: initializes  paw_mult
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_init()
      implicit none

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer i,j,k,nfft3d,npack0,nion,nkatm
      integer lm_size
      integer ia,ii,jj
      integer l,m
      integer Gx(2),Gy(2),Gz(2),Ylm(2)
      real*8  scal,gg,fourpi,omega
      complex*16 itol,cscal
    
      !**** external functions ****
      integer  control_ncut,G_indx
      real*8   control_rcut,double_factorial,lattice_omega
      real*8   lattice_unita
      external control_ncut,G_indx
      external control_rcut,double_factorial,lattice_omega
      external lattice_unita


      fourpi = 16.0d0*datan(1.0d0)
      omega = lattice_omega()

      !*** allocate paw mult memory from heap ***
      call D3dB_nfft3d(1,nfft3d)
      nion = ion_nion()
      nkatm = ion_nkatm()
      call Pack_npack(0,npack0)
      lm_size = (paw_basis_max_mult_l()+1)**2

      !*** allocate gk_smooth,gk,and glm ***
      ok =        my_alloc(mt_dbl,npack0,"gk_smooth",gk_smooth)
      ok = ok.and.my_alloc(mt_dbl,npack0*nkatm,"gk",gk)
      ok = ok.and.my_alloc(mt_dcpl,npack0*lm_size,"g_lm",g_lm)

      !*** allocate paw mult arrays ***
      ok = ok.and.my_alloc(mt_int,nion,"i_v_mult",i_v_mult) !** same as i_paw_qlm **
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",0)

      lm_size = 0
      do ii=1,nion
        ia = ion_katm(ii)

        int_mb(i_v_mult(1) + ii - 1) = lm_size
        lm_size = lm_size + (paw_basis_mult_l(ia) + 1)**2
      end do
      ok = my_alloc(mt_dcpl,lm_size,"v_mult",v_mult)
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",1)

      !*** allocate self_energy_coeff and mult_energy_coeff arrays ***
      lm_size  = paw_basis_max_mult_l()
      ok = my_alloc(mt_dcpl,lm_size,"self_energy_coeff",
     >              self_energy_coeff)

      lm_size  = (paw_basis_max_mult_l()+1)**2
      ok = ok.and.
     >     my_alloc(mt_dcpl,lm_size,"mult_energy_coeff",
     >              mult_energy_coeff)


      !**** initialize sigma_smooth and ncut***
      sigma_smooth = control_rcut()
      ncut         = control_ncut()


      !**** initialize gk_smooth, gk, and g_lm ****

      !**** allocate stack memory ****
      ok = MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'Ylm',Ylm(2),Ylm(1))
      if (.not.ok) 
     > call errquit("paw_mult_init: out of stack memory",2)

      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      !**** g_lm ****
      jj = 0
      itol = dcmplx(1.0d0,0.0d0)
      do l=0,paw_basis_max_mult_l()
        k = 2*l+1
        cscal = itol/double_factorial(k)
        itol = itol*dcmplx(0.0d0,-1.0d0)

        do m = -l,l

          !*** generate Ylm ***
          call spher_harmonics_generate(l,m,npack0,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(Ylm(1)))

          !**** (-i)**l*Ylm(k)*|k|**l/(2l+1)!! ****
          do k=1,npack0
            gg = dbl_mb(Gx(1)+k-1)**2
     >         + dbl_mb(Gy(1)+k-1)**2
     >         + dbl_mb(Gz(1)+k-1)**2
            dcpl_mb(g_lm(1)+k-1+(jj)*npack0)
     >            =cscal
     >            *dsqrt(gg)**l
     >            *dcpl_mb(Ylm(1)+k-1)
     >
          end do !*k*
          jj = jj + 1

        end do !*m*
      end do !*l*

      !**** gk_smooth and gk ****
      do k=1,npack0
         gg = dbl_mb(Gx(1)+k-1)**2
     >      + dbl_mb(Gy(1)+k-1)**2
     >      + dbl_mb(Gz(1)+k-1)**2
         scal = 0.25d0 * sigma_smooth**2
         dbl_mb(gk_smooth(1)+k-1) = fourpi*dexp(-scal*gg)
     >                             /omega
         do ia=1,nkatm
           scal = 0.25d0 * paw_basis_sigma(ia)**2
           dbl_mb(gk(1)+k-1+(ia-1)*npack0) = fourpi*dexp(-scal*gg)
     >                                      /omega
         end do !*ia*
      end do !*k*

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(Ylm(2))
      ok = ok.and.MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      if (.not.ok) 
     > call errquit('paw_mult_init:error popping stack',2)



*     **** allocate rcell memory ****
      nshl3d=(2*ncut+1)**3
      ok = my_alloc(mt_dbl,(3*nshl3d),'rcell',rcell)
      if (.not. ok) call errquit('out of heap memory',0)

*     **** get lattice vectors in real space ****
      l=0
      do k=-ncut,ncut
        do j=-ncut,ncut
          do i=-ncut,ncut
             l = l+1
             dbl_mb(rcell(1)+3*(l-1) )
     >                = i*lattice_unita(1,1)
     >                + j*lattice_unita(1,2)
     >                + k*lattice_unita(1,3)
             dbl_mb(rcell(1)+3*(l-1)+1)
     >                = i*lattice_unita(2,1)
     >                + j*lattice_unita(2,2)
     >                + k*lattice_unita(2,3)
             dbl_mb(rcell(1)+3*(l-1)+2)
     >                = i*lattice_unita(3,1)
     >                + j*lattice_unita(3,2)
     >                + k*lattice_unita(3,3)
          end do
        end do
      end do


      return
      end


*************************************************
!
!       Name: paw_mult_end
!
!       Purpose: deallocates heap memory
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_end()
      implicit none

#include "paw_mult_data.fh"
#include "paw_ma.fh"

      !**** local varables ****
      logical ok

      ok =        my_dealloc(i_v_mult)
      ok = ok.and.my_dealloc(v_mult)
      ok = ok.and.my_dealloc(g_lm)
      ok = ok.and.my_dealloc(gk)
      ok = ok.and.my_dealloc(gk_smooth)
      ok = ok.and.my_dealloc(self_energy_coeff)
      ok = ok.and.my_dealloc(mult_energy_coeff)
      ok = ok.and.my_dealloc(rcell)
      if (.not.ok)
     > call errquit("paw_mult_end: error freeing heap",0)

      return
      end 




!**************************************************
!
!       Name: paw_mult_dn_cmp_get
!
!       Purpose: returns dn_cmp and dn_cmp_smooth 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_dn_cmp_get(dn_cmp,
     >                               dn_cmp_smooth)
      implicit none
      complex*16 dn_cmp(*)
      complex*16 dn_cmp_smooth(*)

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer ia,ii,jj,kk,l,m,nion,npack0,mult_l
      integer exi(2),tmp(2),QYlm(2)
      real*8  sum

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'QYlm',QYlm(2),QYlm(1))
      tmp(1) = exi(1)
      tmp(2) = exi(2)
      if (.not.ok)
     > call errquit(
     >  'paw_mult_dn_cmp_get: out of stack memory',0)


      call dcopy(2*npack0,0.0d0,0,dn_cmp,1)
      call dcopy(2*npack0,0.0d0,0,dn_cmp_smooth,1)
      sum = 0.0d0
      nion = ion_nion()
      do ii=1,nion
         ia = ion_katm(ii)

         !**** Define QYlm ****
         mult_l = paw_basis_mult_l(ia)
         jj     = int_mb(i_paw_qlm(1)+ii-1)
         write(*,*) "Qlm(0):",ii,dcpl_mb(paw_qlm(1)+jj)
         sum = sum + dble(dcpl_mb(paw_qlm(1)+jj))
         kk     = 0
         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QYlm(1)),1)
         do l=0,mult_l
         do m=-l,l
           call Pack_cc_zaxpy(0,
     >               dcpl_mb(paw_qlm(1)+jj),
     >               dcpl_mb(g_lm(1)+(kk)*npack0),
     >               dcpl_mb(QYlm(1)))

           jj = jj + 1
           kk = kk + 1
         end do
         end do

         !**** Multiply by Structure Factor ****
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_cc_Mul(0,
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(QYlm(1)))


         !**** add up ncmp_smooth^ii  ****
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
         call Pack_cc_Sum(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp_smooth,
     >                    dn_cmp_smooth)

         !**** add up ncmp^ii ***
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
         call Pack_cc_Sum(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp,
     >                    dn_cmp)


      end do !*ii*

      write(*,*) "Qtotal:",sum*dsqrt(16.0d0*datan(1.0d0))
      sum = 2*dbl_mb(paw_q00(1)) + dbl_mb(paw_q00(1)+1)
      write(*,*) "Qcore :",sum*dsqrt(16.0d0*datan(1.0d0))

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(QYlm(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_dn_cmp_get: error popping stack',1)

      return
      end

!**************************************************
!
!       Name: paw_mult_coeff_set
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_coeff_set(vh,vcmp)
      implicit none
      complex*16 vh(*)
      complex*16 vcmp(*)

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"


      !**** local variables ***
      logical ok
      integer ia,ii,jj,kk,l,m,nion,npack0,mult_l,mabs
      integer  isgn
      integer exi(2),tmp1(2),gls(2),gl(2),t_mult(2)
      real*8  omega
      complex*16 csum1,csum2

      !**** external functions ****
      real*8   lattice_omega
      external lattice_omega

      omega = lattice_omega()

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'tmp1',tmp1(2),tmp1(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'gl',gl(2),gl(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'gls',gls(2),gls(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,v_mult(3),
     >                 't_mult',t_mult(2),t_mult(1))
      if (.not.ok)
     > call errquit(
     >  'paw_mult_coeff_set: out of stack memory',0)

      nion = ion_nion()
      do ii=1,nion
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         ia = ion_katm(ii)

         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(gl(1)))
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(gls(1)))

         mult_l   = paw_basis_mult_l(ia)
         jj       = int_mb(i_v_mult(1) + ii - 1) 
         kk       = 0
         do l=0,mult_l
           do m=-l,l
             call Pack_cc_Mul(0,
     >                    dcpl_mb(gl(1)),
     >                    dcpl_mb(g_lm(1)+(kk)*npack0),
     >                    dcpl_mb(tmp1(1)))
             call Pack_cc_izdot(0,
     >                        vh,
     >                        dcpl_mb(tmp1(1)),
     >                        csum1)
             call Pack_cc_Mul(0,
     >                    dcpl_mb(gls(1)),
     >                    dcpl_mb(g_lm(1)+(kk)*npack0),
     >                    dcpl_mb(tmp1(1)))
             call Pack_cc_izdot(0,
     >                        vcmp,
     >                        dcpl_mb(tmp1(1)),
     >                        csum2)

             !**** v_mult(l,m,ii) = <g_lm^a|vh> + <tg_lm^a|vcmp> ****
             dcpl_mb(t_mult(1)+jj) = (csum1 + csum2)*omega

             jj = jj + 1
             kk = kk + 1
           end do !*m*
         end do !*l*
      end do !*ii*


      !**** unscramble multipole coefficients ****
      write(*,*)
      write(*,*)
      do ii=1,nion
         ia = ion_katm(ii)
         mult_l   = paw_basis_mult_l(ia)
         jj       = int_mb(i_v_mult(1) + ii - 1) 
         do l=0,mult_l
         do m=-l,l
            mabs = abs(m)
            if (mod(mabs,2).eq.0) then
              isgn = 1
            else
              isgn = -1
            end if
           kk = jj - 2*m
           dcpl_mb(v_mult(1)+jj) = 0.5d0*(dcpl_mb(t_mult(1)+jj) 
     >                           + isgn*dconjg(dcpl_mb(t_mult(1)+kk)))

           write(*,*) "v_mult:",ii,l,m,dcpl_mb(v_mult(1)+jj)
           jj = jj + 1
         end do
         end do
      end do
      call D3dB_Vector_SumAll(2*v_mult(3),dcpl_mb(v_mult(1)))
      write(*,*)
      write(*,*)


      !**** deallocate stack memory ****
      ok =        MA_pop_stack(t_mult(2))
      ok = ok.and.MA_pop_stack(gls(2))
      ok = ok.and.MA_pop_stack(gl(2))
      ok = ok.and.MA_pop_stack(tmp1(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_coeff_set: error popping stack',1)

      return
      end

!**************************************************
!
!       Name: paw_mult_vzero
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_vzero(vzero)
      implicit none
      real*8 vzero

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ****
      integer ia,ii,jj,nion
      real*8  fourpi
 
      !**** external functions ****
      real*8   lattice_omega
      external lattice_omega

      fourpi = 16.0d0*datan(1.0d0)
      nion = ion_nion()
      vzero = 0.0d0
      do ii=1,nion
         ia = ion_katm(ii)
         jj = int_mb(i_paw_qlm(1)+ii-1)
         vzero = vzero 
     >         + dble(dcpl_mb(paw_qlm(1)+jj))
     >          *(sigma_smooth**2-paw_basis_sigma(ia)**2)
      end do
      vzero = vzero*fourpi*dsqrt(fourpi)/lattice_omega()/4.0d0

      return
      end


!**************************************************
!
!       Name: find_mult_energy_coeff
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine find_mult_energy_coeff(sigma_smooth,
     >                            nshl3d,rcell,
     >                            nion,mult_l_max,mult_energy_coeff)
      implicit none
      real*8 sigma_smooth
      integer nshl3d
      real*8 rcell(3,nshl3d)
      integer nion,mult_l_max
      complex*16 mult_energy_coeff(nion,mult_l_max,nion,mult_l_max)


#include "paw_gaunt.fh"
#include "paw_hypergeom.fh"
#include "paw_spher_func.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer il,ilp
      integer lshell
      integer ii,ia,jj,ja
      integer l,l1,l2,m,m1,m2
      integer mult_l,mult_l_prime 

      real*8 sigma_a, sigma_a_prime
      real*8 ra(3), ra_prime(3)
      real*8 R

      real*8 sigma4,sigma3,sigma2,sigma1
      real*8 cos_theta,phi
      real*8 gaunt_coeff
      real*8 r_over_s4, r_over_s3, r_over_s2, r_over_s1
      real*8 phase_factor4,phase_factor3
      real*8 phase_factor2,phase_factor1, phase_factor
      real*8 const4_l,const2_l,const1_l,const,twopi

      complex*16 f,Ylm

      !**** external functions ****
      real*8     paw_get_mult_coeff,double_factorial
      complex*16 find_Y_lm
      external   paw_get_mult_coeff,double_factorial
      external   find_Y_lm


      twopi = 8.0d0*datan(1.0d0)

      call dcopy(2*(nion*mult_l)**2,0.0d0,0,mult_energy_coeff,1)

      do ii=1,nion
       ia      = ion_katm(ii)
       mult_l  = paw_basis_mult_l(ia)
       sigma_a = paw_basis_sigma(ia)

       sigma4 = dsqrt((sigma_smooth**2+sigma_a**2)/2)
       sigma2 = dsqrt((sigma_smooth**2+sigma_smooth**2)/2)
       sigma1 = dsqrt((sigma_a**2+sigma_a**2)/2)

        il = 1
        do l1=0,mult_l
          const = 4.0d0*dsqrt(twopi)
     >           / ((2*l1+1)*double_factorial(2*l1+1))
        
          const4_l = const/(sigma4**(2*l1+1))
          const2_l = const/(sigma2**(2*l1+1))
          const1_l = const/(sigma1**(2*l1+1))
          do m1=-l1,l1
            f = (-1)**m1*(-2*const4_l+const2_l+const1_l)
        
            ilp = 1
            do l2=0,mult_l
            do m2=-l2,l2

              if ((l1.eq.l2).and.(m1.eq.(-m2)) ) then
                mult_energy_coeff(ii,il,ii,ilp) = f
              end if

              ilp = ilp+1
            end do !*m2*
            end do !*l2*

            il = il+1
          end do !*m1*
        end do !*l1*

      end do !*ii*



      do ii=1,nion
        ia   = ion_katm(ii)
        ra(1)   = ion_rion(1,ii)
        ra(2)   = ion_rion(2,ii)
        ra(3)   = ion_rion(3,ii)
        sigma_a = paw_basis_sigma(ia)
        mult_l  = paw_basis_mult_l(ia)

        do jj=1,ii-1
        do lshell=1,nshl3d
          ja = ion_katm(jj)
          ra_prime(1)   = ion_rion(1,jj) + rcell(1,lshell)
          ra_prime(2)   = ion_rion(2,jj) + rcell(2,lshell)
          ra_prime(3)   = ion_rion(3,jj) + rcell(3,lshell)
          sigma_a_prime = paw_basis_sigma(ja)
          mult_l_prime  = paw_basis_mult_l(ja)
        
          R = dsqrt( (ra(1)-ra_prime(1))**2
     >              +(ra(2)-ra_prime(2))**2
     >              +(ra(3)-ra_prime(3))**2)

          if (R.lt.(4*sigma_smooth)) then
            sigma4 = dsqrt((sigma_smooth**2+sigma_a_prime**2)/2)
            sigma3 = dsqrt((sigma_a**2+sigma_smooth**2)/2)
            sigma2 = dsqrt((sigma_smooth**2+sigma_smooth**2)/2)
            sigma1 = dsqrt((sigma_a**2+sigma_a_prime**2)/2)

            r_over_s4 = R/sigma4
            r_over_s3 = R/sigma3
            r_over_s2 = R/sigma2
            r_over_s1 = R/sigma1

            cos_theta = (ra(3)-ra_prime(3))/R

            if (      (ra(2).eq.ra_prime(2)) 
     >          .and. (ra(1).eq.ra_prime(1))) then
              phi = 0.0d0
            else
              phi = datan2(ra(2)-ra_prime(2),ra(1)-ra_prime(1))
            end if

            il = 1
            do l1=0,mult_l
            do m1=-l1,l1

              ilp = 1
              do l2=0,mult_l_prime
              do m2=-l2,l2
                      
                f = dcmplx(0.0d0,0.0d0)
                m = m1 + m2
                phase_factor4 = (-1)**(m1+l1)/sigma4**(l1+l2+1)
                phase_factor3 = (-1)**(m1+l1)/sigma3**(l1+l2+1)
                phase_factor2 = (-1)**(m1+l1)/sigma2**(l1+l2+1)
                phase_factor1 = (-1)**(m1+l1)/sigma1**(l1+l2+1)
                phase_factor  = 1
                do l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2

                gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

                f = f + paw_get_mult_coeff(l1,l2,l)*phase_factor
     >                  *gaunt_coeff*find_Y_lm(l,m,cos_theta,phi)
     >                  *(
     >  -r_over_s4**l
     >   *paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s4)**2/2)*phase_factor4
     > - r_over_s3**l
     >   *paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s3)**2/2)*phase_factor3
     > + r_over_s2**l
     >   *paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s2)**2/2)*phase_factor2
     > + r_over_s1**l
     >   *paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s1)**2/2)*phase_factor1
     >                   ) 
                phase_factor = -phase_factor

                end do !*l*
                
                mult_energy_coeff(ii,il,jj,ilp) 
     >           = mult_energy_coeff(ii,il,jj,ilp) + f 

                ilp = ilp + 1
              end do !*m2*
              end do !*l2*

              il = il + 1
            end do !*m1*
            end do !*l1*

          end if

        end do !*lshell*
        end do !*jj*
      end do !*ii*

      do ii=1,nion
      do jj=ii+1,nion
        do il=1,mult_l_max
        do ilp=1,mult_l_max
              mult_energy_coeff(ii,il,jj,ilp) 
     >         = mult_energy_coeff(jj,ilp,ii,il)
        end do
        end do
      end do
      end do

      return
      end






