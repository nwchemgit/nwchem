c
c $Id$
c

*     **********************************
*     *	                               *
*     *    Model_Potential_init        *
*     *                                *
*     **********************************

*   This reads in Gaussian Basis and initializes
*  variables for the SR contribution 
*
      subroutine Model_Potential_init(rtdb)
      implicit none
      integer rtdb

#include "errquit.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "bafdecls.fh"

***** Mode Potential common block *****
#include "model_potential.fh"

*     **** local variables ***
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,ia,element_length,basis_length
      real*8 sqrtpi
      character*4 element
      character*30 rtdb_name
      integer count,tmp(2)

*     **** external functions ****
      character*2 pspw_qmmm_atom
      integer     pspw_qmmm_nkatm
      real*8      pspw_Q_zv_MM
      external    pspw_qmmm_atom
      external    pspw_qmmm_nkatm
      external    pspw_Q_zv_MM




      call Parallel_taskid(taskid)
      value = btdb_parallel(.true.)


*     ***************************
*     **** read in lr_cutoff ****
*     ***************************

*     **** Allocate Memory from heap ****
      count = pspw_qmmm_nkatm()
      value = BA_alloc_get(mt_dbl,count,'lr_cutoff',
     >                       lr_cutoff(2),lr_cutoff(1))

*     **** read lr_cutoff from rtdb ****
      do ia=1,pspw_qmmm_nkatm()
        element = '    '
        element = pspw_qmmm_atom(ia)

        element_length = index(element,' ') - 1
        rtdb_name = 'pspw_lr_cutoff:'//element(1:element_length)
        if (.not.btdb_get(rtdb,rtdb_name,
     >                   mt_dbl,1,
     >                   dbl_mb(lr_cutoff(1)+ia-1))) then
           dbl_mb(lr_cutoff(1)+ia-1) = 1.0d0
        end if
      end do

      if (taskid.eq.MASTER) then
        write(luout,*)
        write(luout,*) "Model Potential - Long-Range Potential   "
        write(luout,*) "-----------------------------------------"
        write(luout,*)
        write(luout,'  (A15,2x,A15,2x,A4)') "q","rcut","atom"
        do ia=1,pspw_qmmm_nkatm()
           write(luout,'  (E15.6,2x,E15.6,2x,A4)') 
     >         pspw_Q_zv_MM(ia),
     >         dbl_mb(lr_cutoff(1)+ia-1),
     >         pspw_qmmm_atom(ia)
        end do
        write(luout,*)
      end if


*     **************************
*     **** read in SR Basis ****
*     **************************

*     **** determine Nbasis = number of total basis functions ****
      Nbasis = 0
      do ia=1,pspw_qmmm_nkatm()
        element = '    '
        element = pspw_qmmm_atom(ia)

        element_length = index(element,' ') - 1
        rtdb_name = 'pspw_SR_basis_length:'//element(1:element_length)
        if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,basis_length)) then
          basis_length = 0
        end if
        Nbasis = Nbasis + basis_length
      end do

      if (Nbasis.gt.0) then

*       **** Allocate Memory from heap ****
        value = BA_alloc_get(mt_int,Nbasis,'iatom',
     >                       iatom(2),iatom(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'c',c(2),c(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'g',g(2),g(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'coeff',coeff(2),coeff(1))
        if (.not. value) 
     >   call errquit('Model_Potentail_init:out of heap memory',0,
     &       MA_ERR)


*       **** Read basis from rtdb ****
        count = 0
        do ia=1,pspw_qmmm_nkatm()
          element = '    '
          element = pspw_qmmm_atom(ia)

          element_length = index(element,' ') - 1
          rtdb_name = 'pspw_SR_basis_length:'//element(1:element_length)
          if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,basis_length)) then
            basis_length = 0
          end if

          if (basis_length.gt.1) then
            value = value.and.
     >         BA_push_get(mt_dbl,(2*basis_length),'tmp',tmp(2),tmp(1))
              if (.not. value)
     >         call errquit(
     >        'Model_Potential_init:failed pushing stack:',0, MA_ERR)

            rtdb_name = 'pspw_SR_basis:'//element(1:element_length)
            value = btdb_get(rtdb,rtdb_name,
     >                 mt_dbl,(2*basis_length),dbl_mb(tmp(1)))


            do i=1,basis_length
              int_mb(iatom(1)+count+i-1) = ia
              dbl_mb(c(1)    +count+i-1) = dbl_mb(tmp(1)+2*(i-1)  )
              dbl_mb(g(1)    +count+i-1) = dbl_mb(tmp(1)+2*(i-1)+1) 
            end do
            value = BA_pop_stack(tmp(2))
            if (.not. value)
     >       call errquit(
     >      'Model_Potential_init:failed popping stack:',0, MA_ERR)
          end if

          count = count + basis_length
        end do

        if (taskid.eq.MASTER) then
           write(luout,*)
           write(luout,*) "Model Potential - Gaussian Basis   "
           write(luout,*) "-----------------------------------"
           write(luout,*)
           write(luout,'  (A15,2x,A15,2x,A4)') "coeff","exponent","atom"
           do i=1,Nbasis
              write(luout,'  (E15.6,2x,E15.6,2x,A4)') 
     >            dbl_mb(c(1)+i-1),
     >            dbl_mb(g(1)+i-1),
     >            pspw_qmmm_atom(int_mb(iatom(1)+i-1))
           end do
           write(luout,*)
        end if

         sqrtpi = dsqrt(4.0d0*datan(1.0d0))
         do i=1,Nbasis
            dbl_mb(coeff(1)+i-1) = dbl_mb(c(1)+i-1)
     >                           * (1.0d0/(sqrtpi*dbl_mb(g(1)+i-1)))**3
            dbl_mb(g(1)+i-1) = 1.0d0/dbl_mb(g(1)+i-1)**2
         end do

      end if


      return
      end

*     **********************************
*     *	                               *
*     *     Model_Potential_end	       *
*     *                                *
*     **********************************

      subroutine Model_Potential_end()

#include "bafdecls.fh"

***** Mode Potential common block *****
#include "model_potential.fh"

*     **** local variables ***
      logical value

      value = BA_free_heap(lr_cutoff(2))
      if (Nbasis.gt.0) then
         value = value.and.BA_free_heap(iatom(2))
         value = value.and.BA_free_heap(c(2))
         value = value.and.BA_free_heap(g(2))
         value = value.and.BA_free_heap(coeff(2))
      end if
      if (.not.value) call errquit('error BA_free_heap',0, MA_ERR)

      return
      end




*     **********************************
*     *	                               *
*     *   Model_Potential_Generate_V   *
*     *                                *
*     **********************************

*   This routine calculate the Model Potential
*
      subroutine Model_Potential_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,n2ft3d)
      real*8 Vqm(n2ft3d)

#include "bafdecls.fh"

***** Mode Potential common block *****
#include "model_potential.fh"


*     **** local variables ****
      integer ii,i,k,ia
      real*8 r1(3)
      real*8 q,d,d2,rcut
      real*8 sqrtpi

*     **** external functions ****
      integer  pspw_qmmm_katm,pspw_qmmm_nion
      real*8   pspw_qmmm_rion,pspw_Q_zv_MM,util_erf
      external pspw_qmmm_katm,pspw_qmmm_nion
      external pspw_qmmm_rion,pspw_Q_zv_MM,util_erf


      sqrtpi = dsqrt(4.0d0*datan(1.0d0))

      do ii=1,pspw_qmmm_nion()
         r1(1) = pspw_qmmm_rion(1,ii)
         r1(2) = pspw_qmmm_rion(2,ii)
         r1(3) = pspw_qmmm_rion(3,ii)
         ia =pspw_qmmm_katm(ii)

*        **** Add the SR contribution to the model potential ****
         do i=1,Nbasis
           if (int_mb(iatom(1)+i-1).eq.ia) then
             do k=1,n2ft3d
                d2 = (rgrid(1,k)-r1(1))**2 
     >             + (rgrid(2,k)-r1(2))**2
     >             + (rgrid(3,k)-r1(3))**2 

                Vqm(k) = Vqm(k) 
     >                 + dbl_mb(coeff(1)+i-1)*dexp(-d2*dbl_mb(g(1)+i-1))
             end do
           end if
         end do


*        **** Add the LR contribution to the Model Potential ****
*        **** get charges from pspw_Q data structure         ****
         q    = pspw_Q_zv_MM(ia)
         rcut = dbl_mb(lr_cutoff(1)+ia-1)
         if (abs(q).gt.0.0d0) then
           do k=1,n2ft3d
             d = dsqrt( (rgrid(1,k)-r1(1))**2
     >                + (rgrid(2,k)-r1(2))**2
     >                + (rgrid(3,k)-r1(3))**2)

             if (d.gt.small) then
               Vqm(k) = Vqm(k) - q*util_erf(d/rcut)/d
             else
               Vqm(k) = Vqm(k) - 2.0d0*q/(sqrtpi*rcut)
             end if
           end do
        end if


      end do

      return 
      end


*     **********************************
*     *	                               *
*     *      Model_Potential_Fmm       *
*     *                                *
*     **********************************

*     This routine calculates the force on the MM ions from the
*   Model potential interaction with the QM density and adds it
*   to the MM ion force, f_mm.
*
*     Exit - f_mm: force on the MM ions
*
*     Uses - pspw_qmmm_nion,pspw_qmmm_rion
*
*     Author - Eric Bylaska
*
      subroutine Model_Potential_Fmm(n2ft3d,rgrid,rho,dv,f_mm)
      implicit none
      integer n2ft3d
      real*8  rgrid(3,*),rho(*)
      real*8  dv

      real*8 f_mm(3,*)

#include "bafdecls.fh"

***** Mode Potential common block *****
#include "model_potential.fh"

*     **** local variables ****
      integer ii,i,k,ia
      real*8  r1(3),ftmp(3),desr
      real*8  x,y,z,d2,d
      real*8  q,rcut,sqrtpi

*     **** external functions ****
      integer  pspw_qmmm_katm,pspw_qmmm_nion
      real*8   pspw_qmmm_rion,pspw_Q_zv_MM,util_erf
      external pspw_qmmm_katm,pspw_qmmm_nion
      external pspw_qmmm_rion,pspw_Q_zv_MM,util_erf


      sqrtpi = dsqrt(4.0d0*datan(1.0d0))
      do ii=1,pspw_qmmm_nion()
         r1(1) = pspw_qmmm_rion(1,ii)
         r1(2) = pspw_qmmm_rion(2,ii)
         r1(3) = pspw_qmmm_rion(3,ii)
         ia =pspw_qmmm_katm(ii)

         call dcopy(3,0.0d0,0,ftmp,1)

*        **** Add the SR contribution ****
         do i=1,Nbasis
           if (int_mb(iatom(1)+i-1).eq.ia) then
             do k=1,n2ft3d
                x = (rgrid(1,k)-r1(1))
                y = (rgrid(2,k)-r1(2))
                z = (rgrid(3,k)-r1(3))
                d2 = (x*x+y*y+z*z)

                desr = dv*rho(k)
     >                 * dbl_mb(coeff(1)+i-1)
     >                 * dexp(-d2*dbl_mb(g(1)+i-1))
     >                 * (-2.0d0*dbl_mb(g(1)+i-1))
                ftmp(1) = ftmp(1) + x*desr
                ftmp(2) = ftmp(2) + y*desr
                ftmp(3) = ftmp(3) + z*desr
             end do
           end if
         end do

*        **** Add the LR contribution ****
         q    = pspw_Q_zv_MM(ia)
         rcut = dbl_mb(lr_cutoff(1)+ia-1)
         if (abs(q).gt.0.0d0) then
           do k=1,n2ft3d
             x = (rgrid(1,k)-r1(1))
             y = (rgrid(2,k)-r1(2))
             z = (rgrid(3,k)-r1(3))
             d = dsqrt(x*x+y*y+z*z)
             if (d.gt.small) then
                desr = dv*rho(k)*q
     >                 * ((-2.0d0/sqrtpi)
     >                     *(d/rcut)
     >                     *dexp(-(d/rcut)**2)
     >                   + util_erf((d/rcut)))
     >                 /(d**3)
                ftmp(1) = ftmp(1) + x*desr
                ftmp(2) = ftmp(2) + y*desr
                ftmp(3) = ftmp(3) + z*desr
             end if
           end do
        end if


*        **** add all force contributions and add to f_mm ****
         call D3dB_Vector_SumAll(3,ftmp)
         f_mm(1,ii) = f_mm(1,ii) + ftmp(1)
         f_mm(2,ii) = f_mm(2,ii) + ftmp(2)
         f_mm(3,ii) = f_mm(3,ii) + ftmp(3)


      end do

      return
      end
