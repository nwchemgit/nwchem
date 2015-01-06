

**** This file contains debug and old versions of code *****


*     ***********************************
*     *					*
*     *	 	   v_nonlocal_old       *
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 5-13-2002 to improve
*         parallel efficiency.
*
      subroutine v_nonlocal_old(ispin,ne,psi1,psi2,move,fion,
     >                          fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"


*     *** local variables ***
      integer G(3),npack1,nion,nu
      integer i,j,ii,ia,l,n,nn
      integer k,shift,l_prj,m_prj,nproj,Gijl_indx
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),xtmp(2),sw1(2),sw2(2),sw3(2),sum(2)
      logical value,sd_function
      real*8 vmm(50)
      integer ld_ptr

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('v_nonlocal: out of stack',0, MA_ERR)

      if (move) then
       value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
       value = value.and.BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
       if (.not. value) 
     >  call errquit('v_nonlocal:out of stack memory',1,MA_ERR)

       G(1)  = Pack_G_indx(1,1)
       G(2)  = Pack_G_indx(1,2)
       G(3)  = Pack_G_indx(1,3)
      end if

      omega = lattice_omega()
      scal = 1.0d0/(omega)

      do ii=1,nion
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
c        call strfac(ii,dcpl_mb(exi(1)))
c        call Pack_c_pack(1,dcpl_mb(exi(1)))
        call strfac_pack(1,ii,dcpl_mb(exi(1)))


*       **** generate sw1's and projectors ****
        do l=1,nproj

           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))
           m_prj = int_mb(m_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))


           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)

*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****


*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))

           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

*           ***** scale psp by factor - used for generating antiferromagnetic structures ****
*           **** nwchem input: pspspin up/down scale l ion_numbers                       ****
            if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999)))
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw1(1)+(l-1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw1(1)+(l-1)*nn+ne(1)),1)
            end if
        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Gijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
        call Multiply_Gijl_sw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl_indx),
     >                         dbl_mb(sw1(1)),
     >                         dbl_mb(sw2(1)))

*       **** do Kleinman-Bylander Multiplication ****
        !scal = 1.0d0/(omega)
        call dscal(nn*int_mb(nprj(1)+ia-1),
     >             scal,dbl_mb(sw2(1)),1)

*       **** add xc and coulomb paw parts to sw2 ***
        if ((int_mb(psp_type(1)+ia-1).eq.4)) then

*          **** sw2 = sw2 + Vxcijl*sw1 ******
           call nwpw_xc_solve(ii,ia,
     >        int_mb(n1dgrid(1)+ia-1),
     >        int_mb(n1dbasis(1)+ia-1),
     >        dbl_mb(psi_data_get_chnk(int_mb(phi_ae(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(phi_ps(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(dphi_ae(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(dphi_ps(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(core_ae(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(core_ps(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(core_ae_prime(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(core_ps_prime(1)+ia-1))),
     >        dbl_mb(psi_data_get_chnk(int_mb(rgrid(1)+ia-1))),
     >        dbl_mb(log_amesh(1)+ia-1),
     >        ispin,ne,int_mb(nprj(1)+ia-1),
     >        dbl_mb(sw1(1)),dbl_mb(sw2(1)))

        end if

        call DGEMM('N','T',2*npack1,nn,int_mb(nprj(1)+ia-1),
     >             (-1.0d0),
     >             dcpl_mb(prjtmp(1)), 2*npack1,
     >             dbl_mb(sw2(1)),     nn,
     >             (1.0d0),
     >             psi2,               2*npack1)


        if (move) then
        do l=1,nproj
             do n=1,nn
                if (ispin.eq.1) 
     >            dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >            =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
                do i=1,npack1
                   ctmp = psi1(i+(n-1)*npack1)
     >                  *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1 + i-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do
                call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+1+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+2+3*(n-1)))

             end do

             call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

             !**** fractional weighting ****
             if (fractional) then
              do n=1,nn
               call Dneall_qton(n,i)
               dbl_mb(sum(1)+3*(n-1))=dbl_mb(sum(1)+3*(n-1))*occ(i)
               dbl_mb(sum(1)+1+3*(n-1))=dbl_mb(sum(1)+1+3*(n-1))*occ(i)
               dbl_mb(sum(1)+2+3*(n-1))=dbl_mb(sum(1)+2+3*(n-1))*occ(i)
              end do
             end if
  
             ff(1) = 0.0d0
             ff(2) = 0.0d0
             ff(3) = 0.0d0
             do n=1,nn
                ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+  3*(n-1))
                ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+1+3*(n-1))
                ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(sum(1)+2+3*(n-1))
             end do
             call D1dB_Vector_SumAll(3,ff)
             fion(1,ii) = fion(1,ii)  + ff(1)
             fion(2,ii) = fion(2,ii)  + ff(2)
             fion(3,ii) = fion(3,ii)  + ff(3)
       
        end do !** l **
        end if !** move **


        end if !** nproj>0 **
      end do !** ii **

      value = .true.
      if (move) then
      value = value.and.BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      end if
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('v_nonlocal: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      return 
      end


*     ***********************************
*     *					*
*     *	 	   f_vnonlocal_old1	*
*     *					*
*     ***********************************

      subroutine f_vnonlocal_old1(ispin,ne,psi1,fion,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"
cccccccccccc#include "frac_occ.fh"

*     *** local variables ***
      integer G(3),npack1,shift,Gijl_indx
      integer i,ii,ia,k,l,n,nn,l_prj,m_prj,nproj
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),vtmp(2),xtmp(2),sw1(2),sw2(2),sum(2)
c      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx,psi_data_get_ptr
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx,psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, npack1,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      value = value.and.
     >      BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) 
     >  call errquit('f_vnonlocal: out of stack memory',0, MA_ERR)

c     **** define Gx,Gy and Gz in packed space ****
      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

        do l=1,nproj

           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  +(ia-1)*(nmax_max*lmmax_max))
           m_prj = int_mb(m_projector(1)+(l-1) 
     >                                  +(ia-1)*(nmax_max*lmmax_max))

           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))
           end if

           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

*           ***** scale psp by factor - used for generating antiferromagnetic structures ****
*           **** nwchem input: pspspin up/down scale l ion_numbers                       ****
            if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999)))
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw1(1)+(l-1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw1(1)+(l-1)*nn+ne(1)),1)
            end if

        end do
        call D3dB_Vector_Sumall((nn*nproj),dbl_mb(sw1(1)))
 

*       **** sw2 = Gijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
        call Multiply_Gijl_sw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl_indx),
     >                         dbl_mb(sw1(1)),
     >                         dbl_mb(sw2(1)))




*       **** do Kleinman-Bylander Multiplication ****
        do l=1,nproj


           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))

           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(vtmp(1)))
           end if


          scal = 1.0d0/(omega)
          call dscal(nn,scal,dbl_mb(sw2(1)+(l-1)*nn),1)


          do n=1,nn
             if (ispin.eq.1) 
     >         dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >         =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
             do i=1,npack1
                ctmp = psi1(i+(n-1)*npack1)
     >               *dconjg(dcpl_mb(vtmp(1)+i-1))
                dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
             end do
             call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+3*(n-1)))
             call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+1+3*(n-1)))
             call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                         dbl_mb(sum(1)+2+3*(n-1)))

          end do

          call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

          !**** fractional weighting ****
          if (fractional) then
           do n=1,nn
            call Dneall_qton(n,i)
            dbl_mb(sum(1)+3*(n-1))  =dbl_mb(sum(1)  +3*(n-1))*occ(i)
            dbl_mb(sum(1)+1+3*(n-1))=dbl_mb(sum(1)+1+3*(n-1))*occ(i)
            dbl_mb(sum(1)+2+3*(n-1))=dbl_mb(sum(1)+2+3*(n-1))*occ(i)
           end do
          end if
  
          ff(1) = 0.0d0
          ff(2) = 0.0d0
          ff(3) = 0.0d0
          do n=1,nn
             ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+3*(n-1))
             ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+1+3*(n-1))
             ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                            *dbl_mb(sum(1)+2+3*(n-1))
          end do
          call D1dB_Vector_SumAll(3,ff)
          fion(1,ii) = fion(1,ii) + ff(1)
          fion(2,ii) = fion(2,ii) + ff(2)
          fion(3,ii) = fion(3,ii) + ff(3)
       
        end do !** l **

        end if !** nproj>0 **

      end do !** ii **

      value = BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('f_vnonlocal: popping stack memory',0, MA_ERR)


      call nwpw_timing_end(6)
      return 
      end



*     ***********************************
*     *					*
*     *	       f_vnonlocal_old2         *
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 12-1-2013 to handle PAW operators.
*
*  To Do -  For very large numbers of atoms the code will need to be restructured
*           to distribute the sw1a and sw2a matrices over np_i.  Basically, if the orthogonalization matrices
*           need to be distributed then sw1a and sw2a will need to be distributed as well.  A simple algorithm
*           to do this will be to keep the loop structure the same but instead of the
*           call to D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1))), this will need to be, 
*           changed to D3dB_Vector_SumAll(nn*nproj,dbl_mb(sw1(1))) and then place sw1 --> sw1a(distributed)
*
      subroutine f_vnonlocal_old2(ispin,ne,psi1,fion,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer G(3),npack1,nion,nu
      integer i,j,ii,ia,l,n,nn,l1,ip,iip,iipmax,jp,swstart
      integer k,shift,l_prj,m_prj,nproj,Gijl_indx
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),xtmp(2),sum(2)
      logical value,sd_function
      real*8 vmm(50)
      integer ld_ptr

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega,ddot
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega,ddot

      call nwpw_timing_start(6) 


*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      call psp_prj_indx_alloc_sw1a_sw2a(nn)
      value = BA_push_get(mt_dcpl,npack1,'exi',exi(2), exi(1))
      value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
      value = value.and.BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not.value) 
     >  call errquit('f_vnonlocal:out of stack',0, MA_ERR)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()
      scal = 1.0d0/(omega)

      jp = 0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)

*        **** structure factor and local pseudopotential ****
         call strfac_pack(1,ii,dcpl_mb(exi(1)))

         do l=1,nproj
            shift       = int_mb(shift_prj_indx(1)+jp)
            sd_function = log_mb(sd_function_prj_indx(1)+jp)
            jp = jp + 1

*           **** phase factor does not matter therefore ****
*           **** (-i)^l is the same as (i)^l in the     ****
*           **** Rayleigh scattering formula            ****

*           *** current function is s or d ****
            if (sd_function) then
               call Pack_tc_Mul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

*           *** current function is p or f ****
            else
               call Pack_tc_iMul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

            end if
            call Pack_cc_indot(1,nn,
     >                       psi1,
     >                       dcpl_mb(prjtmp(1)),
     >                       dbl_mb(sw1a(1)+(jp-1)*nn))
         end do
      end do
      call D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1)))


*     **** Compute sw2  ****
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)
         swstart     = int_mb(swstart_prj_indx(1)+ip-1)

*        **** sw2 = Gijl*sw1 ******
         Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
         call Multiply_Gijl_sw1(nn,
     >                          nproj,
     >                          int_mb(nmax(1)+ia-1),
     >                          int_mb(lmax(1)+ia-1),
     >                          int_mb(n_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(l_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(m_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          dbl_mb(Gijl_indx),
     >                          dbl_mb(sw1a(1)+swstart*nn),
     >                          dbl_mb(sw2a(1)+swstart*nn))

      end do

*     **** do Kleinman-Bylander Multiplication ****
      call dscal(nn*n_prj_indx,scal,dbl_mb(sw2a(1)),1)

*     **** apply the sw2 to psi ****
      jp  = 0
      do iip=1,nion_prj_indx,nprj_mult 

         swstart = int_mb(swstart_prj_indx(1)+iip-1)
         l1      = 0
         iipmax  = (iip+nprj_mult-1)
         if (iipmax.gt.nion_prj_indx) iipmax = nion_prj_indx

         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

*           **** structure factor and local pseudopotential ****
            call strfac_pack(1,ii,dcpl_mb(exi(1)))

            do l=1,nproj
               shift       = int_mb(shift_prj_indx(1)+jp)
               l_prj       = int_mb(l_prj_prj_indx(1)+jp)
               m_prj       = int_mb(m_prj_prj_indx(1)+jp)
               sd_function = log_mb(sd_function_prj_indx(1)+jp)
               jp = jp + 1

*              **** phase factor does not matter therefore ****
*              **** (-i)^l is the same as (i)^l in the     ****
*              **** Rayleigh scattering formula            ****

*              *** current function is s or d ****
               if (sd_function) then
                  call Pack_tc_Mul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))

*              *** current function is p or f ****
               else
                  call Pack_tc_iMul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))
               end if

*              ***** scale (sw2a) psp by factor - used for generating antiferromagnetic structures ****
*              **** nwchem input: pspspin up/down scale l ion_numbers                              ****
               if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999)))
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(l-1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(l-1)*nn+ne(1)),1)
               end if

               l1 = l1 + 1
            end do
         end do

         l1 = 0
         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

            do l=1,nproj
               do n=1,nn
                  call Pack_cct_iconjgMul(1,
     >                               dcpl_mb(prjtmp(1)+l1*npack1),
     >                               psi1(1+(n-1)*npack1),
     >                               dbl_mb(xtmp(1)))
                 call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+1+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+2+3*(n-1)))

               end do
               call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

               !**** fractional weighting ****
               if (fractional) then
                do n=1,nn
                 call Dneall_qton(n,i)
                 dbl_mb(sum(1)+3*(n-1))  
     >              =dbl_mb(sum(1)  +3*(n-1))*occ(i)
                 dbl_mb(sum(1)+1+3*(n-1))
     >              =dbl_mb(sum(1)+1+3*(n-1))*occ(i)
                 dbl_mb(sum(1)+2+3*(n-1))
     >              =dbl_mb(sum(1)+2+3*(n-1))*occ(i)
                end do
               end if
 
               ff(1) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)),3)
               ff(2) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+1),3)
               ff(3) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+2),3)
               call D1dB_Vector_SumAll(3,ff)
               fion(1,ii) = fion(1,ii)  + ff(1)*(3-ispin)
               fion(2,ii) = fion(2,ii)  + ff(2)*(3-ispin)
               fion(3,ii) = fion(3,ii)  + ff(3)*(3-ispin)

               l1 = l1 + 1
            end do !** l **
         end do !** ip **

      end do !** iip **

*     *** add multipole force here ***
      if (pawexist)
     >   call nwpw_compcharge_F_multipole_all(ispin,dbl_mb(zv(1)),fion)

      value =           BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('f_vnonlocal: popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)
      return 
      end



*     ***********************************
*     *					*
*     *	 	  psp_check_print 	*
*     *					*
*     ***********************************
      subroutine psp_check_print(ia)      
      implicit none
      integer ia

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "psp.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer i,j,r,s,l,n,nbasis,indx,indx0

*     **** external functions ****
      logical  control_print
      external control_print
      integer  psp_psp_type,psp_lmax,psp_n1dbasis,psp_nmax
      external psp_psp_type,psp_lmax,psp_n1dbasis,psp_nmax
      integer  psi_data_get_ptr
      external psi_data_get_ptr

      call Parallel_taskid(taskid)
      oprint= (taskid.eq.MASTER)

      nbasis = psp_n1dbasis(ia)
      n   = psp_nmax(ia)
      if (oprint) then

*        **** print vcore matrix ****
         if (psp_psp_type(ia).eq.4) then
            do l=0,2*psp_lmax(ia)

              indx0 = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
              write(*,*) "overlap matrix: ia=",ia,"  l=",l
              do i=1,n
                 indx = indx0 + l*n*n + (i-1)
                 write(*,'(10E11.3)') (dbl_mb(indx+(j-1)*n),j=1,n)
              end do
              write(*,*) 

c              indx0 = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),3)
c              write(*,*) "kinetic matrix: ia=",ia,"  l=",l
c              do i=1,n
c                 indx = indx0 + l*n*n + (i-1)
c                 write(*,'(10E11.3)') (dbl_mb(indx+(j-1)*n),j=1,n)
c              end do
c              write(*,*) 
c
c              indx0 = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),4)
c              write(*,*) "vpseudo matrix: ia=",ia,"  l=",l
c              do i=1,n
c                 indx = indx0 + l*n*n + (i-1)
c                 write(*,'(10E11.3)') (dbl_mb(indx+(j-1)*n),j=1,n)
c              end do
c              write(*,*) 

c               indx0 = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),5)
c               write(*,*) "vcore matrix: ia=",ia,"  l=",l
c               do i=1,n
c                  indx = indx0 + l*n*n + (i-1)
c                  write(*,'(10E11.3)') (dbl_mb(indx+(j-1)*n),j=1,n)
c               end do
c               write(*,*) 
               
c               write(*,*) "comp_charge matrix: ia=",ia,"  l=",l
c               indx0 = 
c     >          psi_data_get_ptr(int_mb(comp_charge_matrix(1)+ia-1),l+1)
c               do i=1,nbasis
c                  indx = indx0 + (i-1)
c                  write(*,'(10E13.3)') (dbl_mb(indx+(j-1)*nbasis),
c     >                                   j=1,nbasis)
c               end do
c               write(*,*) 

c               write(*,*) "comp_pot matrix: ia=",ia,"  l=",l
c               indx0 = 
c     >          psi_data_get_ptr(int_mb(comp_pot_matrix(1)+ia-1),l+1)
c               do i=1,nbasis
c                  indx = indx0 + (i-1)
c                  write(*,'(10E13.3)') (dbl_mb(indx+(j-1)*nbasis),
c     >                                   j=1,nbasis)
c               end do
c               write(*,*) 

c               write(*,*) "hartree matrix: ia=",ia,"  l=",l
c               indx0 = 
c     >          psi_data_get_ptr(int_mb(hartree_matrix(1)+ia-1),l+1)
c               do i=1,nbasis
c                  do j=1,nbasis
c                     do r=1,nbasis
c                        indx = indx0 + (i-1)*nbasis*nbasis*nbasis 
c     >                               + (j-1)*nbasis*nbasis 
c     >                               + (r-1)*nbasis
c                        write(*,'(10E13.3)') (dbl_mb(indx+(s-1)),
c     >                                        s=1,nbasis)
c                     end do
c                  end do
c               end do
c               write(*,*) 

            end do
         end if

      end if

      return
      end 



*     ***********************************
*     *					*
*     *	       f_vnonlocal_check        *
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 12-1-2013 to handle PAW operators.
*
*  To Do -  For very large numbers of atoms the code will need to be restructured
*           to distribute the sw1a and sw2a matrices over np_i.  Basically, if the orthogonalization matrices
*           need to be distributed then sw1a and sw2a will need to be distributed as well.  A simple algorithm
*           to do this will be to keep the loop structure the same but instead of the
*           call to D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1))), this will need to be, 
*           changed to D3dB_Vector_SumAll(nn*nproj,dbl_mb(sw1(1))) and then place sw1 --> sw1a(distributed)
*
      subroutine f_vnonlocal_check(ispin,ne,psi1,fractional,occ)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer pcount,taskid_j,np_j
      integer G(3),npack0,npack1,nion,nu,mult_l,m,ms,lm
      integer i,j,ii,ia,l,n,nn,l1,ip,iip,iipmax,jp,swstart
      integer k,shift,l_prj,m_prj,nproj,Gijl_indx
      integer nx,ny,nz
      real*8  omega,scal,ff(3),scal1,dv
      complex*16 ctmp
      integer exi(2),xtmp(2),sum(2)
      integer dng_cmp(2),dng_cmp_smooth(2)
      integer vcmp(2),vcmp_smooth(2)
      logical value,sd_function,periodic
      real*8 vmm(50),eh_atom
      real*8 fion(3,50)
      integer ld_ptr

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega,ddot
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega,ddot
      integer  nwpw_compcharge_mult_l,control_version
      external nwpw_compcharge_mult_l,control_version
      real*8   nwpw_compcharge_Qlm
      external nwpw_compcharge_Qlm


      call dcopy(3*50,0.0d0,0,fion,1)

      call nwpw_timing_start(6) 

      periodic = (control_version().eq.3)
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      call psp_prj_indx_alloc_sw1a_sw2a(nn)
      value = BA_push_get(mt_dcpl,npack1,'exi',exi(2), exi(1))
      if (.not.value) 
     >  call errquit('v_nonlocal:out of stack',0, MA_ERR)

      value = value.and.BA_push_get(mt_dbl,npack1,
     >                              'xtmp',xtmp(2),xtmp(1))
      value = value.and.BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) 
     > call errquit('v_nonlocal:out of stack',1,MA_ERR)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()
      scal = 1.0d0/(omega)

      jp = 0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)

*        **** structure factor and local pseudopotential ****
         call strfac_pack(1,ii,dcpl_mb(exi(1)))

         do l=1,nproj
            shift       = int_mb(shift_prj_indx(1)+jp)
            sd_function = log_mb(sd_function_prj_indx(1)+jp)
            jp = jp + 1

*           **** phase factor does not matter therefore ****
*           **** (-i)^l is the same as (i)^l in the     ****
*           **** Rayleigh scattering formula            ****

*           **** phase fact DOES matter for compensation charge!!!!     ****
*           **** assume that sign factor for proj is in kbpp formatting ****

*           *** current function is s or d ****
            if (sd_function) then
               call Pack_tc_Mul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

*           *** current function is p or f ****
            else
               call Pack_tc_iMul(1,dbl_mb(shift),
     >                            dcpl_mb(exi(1)),
     >                            dcpl_mb(prjtmp(1)))

            end if
            call Pack_cc_indot(1,nn,
     >                       psi1,
     >                       dcpl_mb(prjtmp(1)),
     >                       dbl_mb(sw1a(1)+(jp-1)*nn))
         end do
      end do
      call D3dB_Vector_SumAll(nn*n_prj_indx,dbl_mb(sw1a(1)))


*     **** Compute sw2  ****
      eh_atom = 0.0d0
      do ip=1,nion_prj_indx
         ii          = int_mb(ii_prj_indx(1)+ip-1)
         ia          = int_mb(ia_prj_indx(1)+ip-1)
         nproj       = int_mb(nproj_prj_indx(1)+ip-1)
         swstart     = int_mb(swstart_prj_indx(1)+ip-1)


*        **** sw2 = Gijl*sw1 ******
         Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
cebug         call Multiply_Gijl_sw1(nn,
cebug     >                          nproj,
cebug     >                          int_mb(nmax(1)+ia-1),
cebug     >                          int_mb(lmax(1)+ia-1),
cebug     >                          int_mb(n_projector(1)
cebug     >                                 + (ia-1)*(nmax_max*lmmax_max)),
cebug     >                          int_mb(l_projector(1)
cebug     >                                 + (ia-1)*(nmax_max*lmmax_max)),
cebug     >                          int_mb(m_projector(1)
cebug     >                                 + (ia-1)*(nmax_max*lmmax_max)),
cebug     >                          dbl_mb(Gijl_indx),
cebug     >                          dbl_mb(sw1a(1)+swstart*nn),
cebug     >                          dbl_mb(sw2a(1)+swstart*nn))

     
*        **** paw operations #1 - generate it's compcharge, add atomic coulomb, and add atomic xc potential ****
         if ((int_mb(psp_type(1)+ia-1).eq.4)) then

            call nwpw_timing_start(39)

*           **** paw atom - generate it's atomic density matrix ****
            call psp_gen_density_matrix(ispin,ne,nproj,
     >                                  dbl_mb(sw1a(1)+swstart*nn),
     >                                  dbl_mb(wtmp(1)))

*           **** paw atom - generate it's compcharge ***
            call nwpw_compcharge_gen_Qlm(ii,ia,ispin,nproj,
     >                                   dbl_mb(wtmp(1)))
            call nwpw_timing_end(39)

*           **** atomic coulomb matrix - sw2 = sw2 + Vhatomijl*sw1  ****
cebug            call nwpw_compcharge_coulomb_atom(ii,ia,ispin,ne,nproj,
cebug     >                                   dbl_mb(wtmp(1)),
cebug     >                                   dbl_mb(sw1a(1)+swstart*nn),
cebug     >                                   dbl_mb(sw2a(1)+swstart*nn),
cebug     >                                   eh_atom)

*           **** xc matrix - sw2 = sw2 + Vxcijl*sw1 ******
cebug            call nwpw_xc_solve(ii,ia,
cebug     >        int_mb(n1dgrid(1)+ia-1),
cebug     >        int_mb(n1dbasis(1)+ia-1),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(phi_ae(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(phi_ps(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(dphi_ae(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(dphi_ps(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(core_ae(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(core_ps(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(core_ae_prime(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(core_ps_prime(1)+ia-1))),
cebug     >        dbl_mb(psi_data_get_chnk(int_mb(rgrid(1)+ia-1))),
cebug     >        dbl_mb(log_amesh(1)+ia-1),
cebug     >        ispin,ne,nproj,
cebug     >        dbl_mb(sw1a(1)+swstart*nn),dbl_mb(sw2a(1)+swstart*nn))

        end if
      end do



*     **** paw operations #2 - generate vcmp,.... and Gaussian Multipole  ****
      if (pawexist) then
         call nwpw_timing_start(39)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         dv = lattice_omega()*scal1

         if (periodic) then
            call Pack_npack(0,npack0)
         else
            call D3dB_n2ft3d(1,npack0)
            npack0 = npack0/2
         end if


         if (use_grid_cmp) then


            value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                          dng_cmp(2),dng_cmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp',
     >                          vcmp(2),vcmp(1))
            if (.not.value) 
     >      call errquit('v_nonlocal:out of stack',4, MA_ERR)

*           **** zero out dE/dQlm array ****
            call nwpw_compcharge_zero_dE_Qlm()

            call nwpw_compcharge_gen_dn_cmp(ispin,dcpl_mb(dng_cmp(1)))

            call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp_tmp(1)))
            call Pack_cc_Sum(0,
     >                       dcpl_mb(vc_tmp(1)),
     >                       dcpl_mb(vcmp_tmp(1)),
     >                       dcpl_mb(vcmp(1)))

            pcount = 0
            do ip=1,nion_prj_indx
               ii          = int_mb(ii_prj_indx(1)+ip-1)
               ia          = int_mb(ia_prj_indx(1)+ip-1)
               if ((int_mb(psp_type(1)+ia-1).eq.4)) then
                  mult_l = nwpw_compcharge_mult_l(ia)
                  do l=0,mult_l
                  do m=-l,l
                     if (mod(pcount,np_j).eq.taskid_j) then
                        call nwpw_compcharge_gen_glm(ii,l,m,
     >                                      dcpl_mb(dng_cmp(1)))
                        call Pack_cc_idot(0,dcpl_mb(dng_cmp(1)),
     >                                      dcpl_mb(vcmp(1)),ff(1))
                        call nwpw_compcharge_add_dE_Qlm(ispin,ii,l,m,
     >                                                  ff(1)*omega)
                     end if
                     pcount = pcount + 1
                  end do
                  end do
               end if
            end do
            value =           BA_pop_stack(vcmp(2))
            value = value.and.BA_pop_stack(dng_cmp(2))
            if (.not.value) 
     >      call errquit('v_nonlocal:popping stack',4, MA_ERR)

         else


            value = BA_push_get(mt_dcpl,npack0,'dng_cmp',
     >                          dng_cmp(2),dng_cmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'dng_cmp_smooth',
     >                          dng_cmp_smooth(2),dng_cmp_smooth(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp',
     >                          vcmp(2),vcmp(1))
            value = value.and.
     >              BA_push_get(mt_dcpl,npack0,'vcmp_smooth',
     >                          vcmp_smooth(2),vcmp_smooth(1))
            if (.not.value) 
     >        call errquit('v_nonlocal:out of stack',4, MA_ERR)

            call nwpw_compcharge_gen_dn_cmp2(ispin,
     >                                       dcpl_mb(dng_cmp(1)),
     >                                       dcpl_mb(dng_cmp_smooth(1)))

            !*** compute hartree potential of ntilde + ncmp_tilde ***
            !*** compute hartree potential ncmp   - ncmp_tilde ***
            if (periodic) then
               call coulomb_v(dcpl_mb(dng_cmp(1)),dcpl_mb(vcmp(1)))
               call Pack_c_Copy(0,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))
               call coulomb_v(dcpl_mb(dng_cmp_smooth(1)),
     >                        dcpl_mb(vcmp_smooth(1)))
               call Pack_cc_Sub2(0,
     >                           dcpl_mb(vcmp_smooth(1)),
     >                           dcpl_mb(vcmp(1)))
               call Pack_cc_Sum2(0,dcpl_mb(vc_tmp(1)),
     >                             dcpl_mb(vcmp_smooth(1)))
               call Pack_c_SMul1(0,omega,dcpl_mb(vcmp(1)))
               call Pack_c_SMul1(0,omega,dcpl_mb(vcmp_smooth(1)))
            else
               !*** dng_cmp(G),dng_cmp_smooth(G) --> dng_cmp(r),dng_cmp_smooth(r) ***
               call Pack_c_unpack(0,dcpl_mb(dng_cmp(1)))
               call Pack_c_unpack(0,dcpl_mb(dng_cmp_smooth(1)))
               call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp(1)))
               call D3dB_cr_pfft3b(1,0,dcpl_mb(dng_cmp_smooth(1)))

               !*** generate vcmp and vcmp_smooth ***
               call coulomb2_v(dcpl_mb(dng_cmp(1)),
     >                         dcpl_mb(vcmp(1)))
               call coulomb2_v(dcpl_mb(dng_cmp_smooth(1)),
     >                         dcpl_mb(vcmp_smooth(1)))
               call D3dB_r_Copy(1,dcpl_mb(vcmp(1)),dcpl_mb(vcmp_tmp(1)))

               !*** vcmp        = vcmp-vcmp_smooth ***
               !*** vcmp_smooth = vcmp_smooth + vc ***
               call D3dB_rr_Sub2(1,
     >                           dcpl_mb(vcmp_smooth(1)),
     >                           dcpl_mb(vcmp(1)))
cdebug               call D3dB_rr_Sum2(1,
cdebug     >                           dcpl_mb(vc_tmp(1)),
cdebug     >                           dcpl_mb(vcmp_smooth(1)))

               !*** vcmp(r),vcmp_smooth(r) --> vcmp(G),vcmp_smooth(G) ***
               !*** May want to change nwpw_compcharge_gen_dE_Qlm to remove these FFTs ***
               call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp(1)))
               call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp_smooth(1)))
               call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp(1)))
               call D3dB_rc_pfft3f(1,0,dcpl_mb(vcmp_smooth(1)))
               call Pack_c_pack(0,dcpl_mb(vcmp(1)))
               call Pack_c_pack(0,dcpl_mb(vcmp_smooth(1)))

            end if

c            call dcopy(2*npack0,dcpl_mb(vc_tmp(1)),1,
c     >                          dcpl_mb(vcmp_smooth(1)),1)
c            call dcopy(2*npack0,0.0d0,0,dcpl_mb(vcmp(1)),1)
c            call D3dB_r_SMul1(1,dv,dcpl_mb(vcmp_smooth(1)))
c            call D3dB_rc_fft3f(1,dcpl_mb(vcmp_smooth(1)))
c            call Pack_c_pack(0,dcpl_mb(vcmp_smooth(1)))

cdebug            call nwpw_compcharge_gen_dE_Qlm(ispin,
cdebug     >                                   dcpl_mb(vcmp_smooth(1)),
cdebug     >                                   dcpl_mb(vcmp(1)),.true.,fion)
cdebug            call write_force(nion,fion,"f_vnonlocal_check a")
         
cdebug            call nwpw_compcharge_gen_dEmult_Qlm(ispin)
cdebug            call nwpw_compcharge_add_dEmult_Qlm(ispin)

            value =           BA_pop_stack(vcmp_smooth(2))
            value = value.and.BA_pop_stack(vcmp(2))
            value = value.and.BA_pop_stack(dng_cmp_smooth(2))
            value = value.and.BA_pop_stack(dng_cmp(2))
            if (.not.value) 
     >        call errquit('v_nonlocal:popping stack',5,MA_ERR)
         end if

*        **** update ncmp*Vloc ****
        call nwpw_compcharge_zero_dE_Qlm()
        call psp_dE_ncmp_vloc_Qlm_test(ispin,.true.,fion)
        call nwpw_compcharge_add_dElocal_Qlm(ispin)
        call write_force(nion,fion,"f_vnonlocal_check b")

         do ip=1,nion_prj_indx
            ii          = int_mb(ii_prj_indx(1)+ip-1)
            ia          = int_mb(ia_prj_indx(1)+ip-1)
            nproj       = int_mb(nproj_prj_indx(1)+ip-1)
            swstart     = int_mb(swstart_prj_indx(1)+ip-1)
            call nwpw_compcharge_gen_sw2(ii,ia,ispin,ne,nproj,
     >                                   dbl_mb(sw1a(1)+swstart*nn),
     >                                   dbl_mb(sw2a(1)+swstart*nn))
         end do
     

         call nwpw_timing_end(39)
      end if


*     **** do Kleinman-Bylander Multiplication ****
      call dscal(nn*n_prj_indx,scal,dbl_mb(sw2a(1)),1)

*     **** apply the sw2 to psi ****
      jp  = 0
      do iip=1,nion_prj_indx,nprj_mult 

         swstart = int_mb(swstart_prj_indx(1)+iip-1)
         l1      = 0
         iipmax  = (iip+nprj_mult-1)
         if (iipmax.gt.nion_prj_indx) iipmax = nion_prj_indx

         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

*           **** structure factor and local pseudopotential ****
            call strfac_pack(1,ii,dcpl_mb(exi(1)))

            do l=1,nproj
               shift       = int_mb(shift_prj_indx(1)+jp)
               l_prj       = int_mb(l_prj_prj_indx(1)+jp)
               m_prj       = int_mb(m_prj_prj_indx(1)+jp)
               sd_function = log_mb(sd_function_prj_indx(1)+jp)
               jp = jp + 1

*              **** phase factor does not matter therefore ****
*              **** (-i)^l is the same as (i)^l in the     ****
*              **** Rayleigh scattering formula            ****

*              *** current function is s or d ****
               if (sd_function) then
                  call Pack_tc_Mul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))

*              *** current function is p or f ****
               else
                  call Pack_tc_iMul(1,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+l1*npack1))
               end if

*              ***** scale (sw2a) psp by factor - used for generating antiferromagnetic structures ****
*              **** nwchem input: pspspin up/down scale l ion_numbers                              ****
               if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_upl(1)+ii-1)) .and.
     >           ((m_prj.ne.int_mb(pspspin_upm(1)+ii-1)).or.
     >            (int_mb(pspspin_upm(1)+ii-1).gt.999))) 
     >            call dscal(ne(1),dbl_mb(pspspin_upscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(l-1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.int_mb(pspspin_downl(1)+ii-1)).and.
     >           ((m_prj.ne.int_mb(pspspin_downm(1)+ii-1)).or.
     >            (int_mb(pspspin_downm(1)+ii-1).gt.999)))
     >            call dscal(ne(2),dbl_mb(pspspin_downscale(1)+ii-1),
     >                       dbl_mb(sw2a(1)+(l-1)*nn+ne(1)),1)
               end if

               l1 = l1 + 1
            end do
         end do

         l1 = 0
         do ip=iip,iipmax
            ii    = int_mb(ii_prj_indx(1)+ip-1)
            ia    = int_mb(ia_prj_indx(1)+ip-1)
            nproj = int_mb(nproj_prj_indx(1)+ip-1)

            do l=1,nproj
               do n=1,nn
                  call Pack_cct_iconjgMul(1,
     >                               dcpl_mb(prjtmp(1)+l1*npack1),
     >                               psi1(1+(n-1)*npack1),
     >                               dbl_mb(xtmp(1)))
                 call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+1+3*(n-1)))
                 call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                             dbl_mb(sum(1)+2+3*(n-1)))

               end do
               call D3dB_Vector_SumAll(3*(nn),dbl_mb(sum(1)))

               !**** fractional weighting ****
               if (fractional) then
                do n=1,nn
                 call Dneall_qton(n,i)
                 dbl_mb(sum(1)+3*(n-1))  
     >              =dbl_mb(sum(1)  +3*(n-1))*occ(i)
                 dbl_mb(sum(1)+1+3*(n-1))
     >              =dbl_mb(sum(1)+1+3*(n-1))*occ(i)
                 dbl_mb(sum(1)+2+3*(n-1))
     >              =dbl_mb(sum(1)+2+3*(n-1))*occ(i)
                end do
               end if
 
               ff(1) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)),3)
               ff(2) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+1),3)
               ff(3) =2.0d0*ddot(nn,dbl_mb(sw2a(1)+(swstart+l1)*nn),1,
     >                              dbl_mb(sum(1)+2),3)
               call D1dB_Vector_SumAll(3,ff)
               fion(1,ii) = fion(1,ii)  + ff(1)*(3-ispin)
               fion(2,ii) = fion(2,ii)  + ff(2)*(3-ispin)
               fion(3,ii) = fion(3,ii)  + ff(3)*(3-ispin)

               l1 = l1 + 1
            end do !** l **
         end do !** ip **

      end do !** ii **

      call write_force(nion,fion,"f_vnonlocal_check c")

*     *** add multipole force here ***
      if (pawexist)
     >  call nwpw_compcharge_F_multipole_test(ispin,dbl_mb(zv(1)),fion)
      call write_force(nion,fion,"f_vnonlocal_check d")

      value =           BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('f_vnonlocal_check: popping stack',3,
     >       MA_ERR)


      call write_force(nion,fion,"f_vnonlocal_check")

      call nwpw_timing_end(6)
      return 
      end



