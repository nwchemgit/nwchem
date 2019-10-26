
*     ******************************
*     *                            *
*     *     rho_symmetrizer_init   *
*     *                            *
*     ******************************
      subroutine rho_symmetrizer_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer p,indx,nops
      integer i,j,k,nxh,nyh,nzh,nfft3d
      real*8  f,msym(3,4)

*     **** rhosymr_block common blocks ****
      integer nx,ny,nz,n2ft3d_serial,n2ft3d
      integer i_grid(2),f_grid(2),sf_grid(2),dn_serial(2),rho0(2)
      logical use_sym,use_translation
      common /rhosymr_block0/ i_grid,f_grid,sf_grid,dn_serial,rho0,
     >                        n2ft3d,n2ft3d_serial,nx,ny,nz,use_sym,
     >                        use_translation

*     **** external functions ****
      logical  control_rho_use_symmetry
      external control_rho_use_symmetry
      integer  ion_sym_number_ops
      external ion_sym_number_ops

      use_sym = control_rho_use_symmetry()
     >         .and.(ion_sym_number_ops().gt.0)

      if (use_sym) then
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         nxh = nx/2
         nyh = ny/2
         nzh = nz/2
         n2ft3d_serial = (nx+2)*ny*nz

         call D3dB_n2ft3d(1,n2ft3d)
         call D3dB_nfft3d(1,nfft3d)

         nops = ion_sym_number_ops()
         use_translation = .false.
         do i=1,nops
            call ion_sym_get_op(i,msym)
            use_translation = use_translation.and.
     >                    (dabs(msym(1,4)).gt.1.0d-6).or.
     >                    (dabs(msym(2,4)).gt.1.0d-6).or.
     >                    (dabs(msym(3,4)).gt.1.0d-6)
         end do

*        **** allocate from heap for rhosym_block common block ****
         value = BA_alloc_get(mt_int,(4*n2ft3d),'i_grid',
     >                        i_grid(2),i_grid(1))
         value = BA_alloc_get(mt_dbl,(3*n2ft3d),'f_grid',
     >                        f_grid(2),f_grid(1))
         value = BA_alloc_get(mt_dbl,(3*n2ft3d),'sf_grid',
     >                        sf_grid(2),sf_grid(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,n2ft3d_serial,'dn_serial',
     >                        dn_serial(2),dn_serial(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,2*n2ft3d,'rho0',
     >                        rho0(2),rho0(1))
         if (.not.value) 
     >      call errquit('rho_symmetrizer_init: out of heap',0, MA_ERR)

         call lattice_r_grid(dbl_mb(sf_grid(1)))
         call lattice_r1_to_frac(n2ft3d,
     >                           dbl_mb(sf_grid(1)),
     >                           dbl_mb(f_grid(1)))
         do i=1,3
            do k=1,n2ft3d
               f = dbl_mb(f_grid(1)+(k-1)*3+i-1) + 0.5d0
               !f = dbl_mb(f_grid(1)+(k-1)*3+i-1) + 1.0d0
               do while (f.ge.1.0d0) 
                  f = f - 1.0d0
               end do
               do while (f.lt.0.0d0) 
                  f = f + 1.0d0
               end do
               dbl_mb(f_grid(1)+(k-1)*3+i-1) = f
            end do
         end do

         call lattice_i_grid(1,int_mb(i_grid(1)))

      end if

      return
      end


*     ******************************
*     *                            *
*     *     rho_symmetrizer_end    *
*     *                            *
*     ******************************
      subroutine rho_symmetrizer_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

*     **** rhosymr_block common blocks ****
      integer nx,ny,nz,n2ft3d_serial,n2ft3d
      integer i_grid(2),f_grid(2),sf_grid(2),dn_serial(2),rho0(2)
      logical use_sym,use_translation
      common /rhosymr_block0/ i_grid,f_grid,sf_grid,dn_serial,rho0,
     >                        n2ft3d,n2ft3d_serial,nx,ny,nz,use_sym,
     >                        use_translation

      if (use_sym) then
         value = BA_free_heap(dn_serial(2))
         value = value.and.BA_free_heap(rho0(2))
         value = value.and.BA_free_heap(i_grid(2))
         value = value.and.BA_free_heap(f_grid(2))
         value = value.and.BA_free_heap(sf_grid(2))
         if (.not.value) 
     >      call errquit('rho_symmetrizer_end: freeing heap',0, MA_ERR)
      end if
      return
      end


*     ******************************************
*     *                                        *
*     *         rho_symmetrizer_gen_serial0    *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_gen_serial0(newdn,n2ft3d,
     >                                  i_grid,sf_grid,
     >                                  rho,
     >                                  nx,ny,nz,n2ft3d_serial,
     >                                  dn_serial)
      implicit none
      logical newdn
      integer n2ft3d
      integer i_grid(4,n2ft3d)
      real*8  sf_grid(3,n2ft3d),rho(n2ft3d)
      integer nx,ny,nz,n2ft3d_serial
      real*8  dn_serial(nx+2,ny,nz)

*     **** local variables ****
      integer i,k,ic1,ic2,ic3
      real*8  c1,c2,c3

      if (newdn) then
        call Parallel_shared_vector_zero(.true.,n2ft3d_serial,dn_serial)
      end if

!$OMP DO
      do k=1,n2ft3d
         if (i_grid(4,k).gt.0) then
            do i=1,3
               do while (sf_grid(i,k).ge.1.0d0)
                  sf_grid(i,k) = sf_grid(i,k) - 1.0d0
               end do
               do while (sf_grid(i,k).lt.0.0d0)
                  sf_grid(i,k) = sf_grid(i,k) + 1.0d0
               end do
               if ((dabs(sf_grid(i,k)-1.0d0).lt.1.0d-6) .or.
     >             (dabs(sf_grid(i,k)).lt.1.0d-6)) then
                  sf_grid(i,k) = 0.0d0
               end if
            end do

            c1 = sf_grid(1,k)*nx
            c2 = sf_grid(2,k)*ny
            c3 = sf_grid(3,k)*nz

            ic1 = NINT(c1) + 1
            ic2 = NINT(c2) + 1
            ic3 = NINT(c3) + 1

            if (newdn) then
               dn_serial(ic1,ic2,ic3) = rho(k)
            else
               dn_serial(ic1,ic2,ic3) = dn_serial(ic1,ic2,ic3) + rho(k)
            end if
         end if
      end do
!$OMP END DO
      return
      end

*     ******************************************
*     *                                        *
*     *         rho_symmetrizer_to_rho         *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_to_rho(n2ft3d,n2ft3d_serial,nx,ny,nz,
     >                                  i_grid,
     >                                  dn_serial,rho)
      implicit none
      integer n2ft3d,n2ft3d_serial,nx,ny,nz
      integer i_grid(4,n2ft3d)
      real*8 dn_serial(nx+2,ny,nz)
      real*8 rho(n2ft3d)

*     **** local variables ****
      integer k,ic1,ic2,ic3,nxh,nyh,nzh

      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      !call Parallel_Vector_SumAll(n2ft3d_serial,dn_serial)
      call D3dB_Vector_SumAll(n2ft3d_serial,dn_serial)

!$OMP DO
      do k=1,n2ft3d
         if (i_grid(4,k).gt.0) then
            ic1 = i_grid(1,k) + nxh
            ic2 = i_grid(2,k) + nyh
            ic3 = i_grid(3,k) + nzh
            rho(k) = dn_serial(ic1+1,ic2+1,ic3+1)
         else
            rho(k) = 0.0d0
         end if
      end do
!$OMP END DO

      return
      end

*     ******************************************
*     *                                        *
*     *        rho_symmetrizer_translate       *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_translate(trans,scal1,rho)
      implicit none
      real*8     trans(3),scal1
      real*8     rho(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer exi(2),nfft3d

      call D3dB_nfft3d(1,nfft3d)
      value =  BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) 
     >   call errquit('rho_symmetrizer_translate:push stack',0,0)

      call D3dB_r_SMul1(1,scal1,rho)
      call D3dB_rc_fft3f(1,rho)
    
      call strfac_translate(trans,dcpl_mb(exi(1)))
      call D3dB_cc_Mul2(1,dcpl_mb(exi(1)),rho)

      call D3dB_cr_fft3b(1,rho)

      value = BA_pop_stack(exi(2))
      if (.not. value) 
     >   call errquit('rho_symmetrizer_translate:pop stack',1,0)
      return
      end


*     ******************************************
*     *                                        *
*     *         rho_symmetrizer_inv33          *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_inv33(m,minv)
      implicit none
      real*8 m(3,3),minv(3,3)

      real*8 det,invdet

      !*** computes the inverse of a matrix m ***
      det = m(1,1) * (m(2,2) * m(3,3) - m(3,2) * m(2,3)) -
     >      m(1,2) * (m(2,1) * m(3,3) - m(2,3) * m(3,1)) +
     >      m(1,3) * (m(2,1) * m(3,2) - m(2,2) * m(3,1))
      invdet = 1.0d0 / det

      minv(1,1) = (m(2,2) * m(3,3) - m(3,2) * m(2,3)) * invdet
      minv(1,2) = (m(1,3) * m(3,2) - m(1,2) * m(3,3)) * invdet
      minv(1,3) = (m(1,2) * m(2,3) - m(1,3) * m(2,2)) * invdet
      minv(2,1) = (m(2,3) * m(3,1) - m(2,1) * m(3,3)) * invdet
      minv(2,2) = (m(1,1) * m(3,3) - m(1,3) * m(3,1)) * invdet
      minv(2,3) = (m(2,1) * m(1,3) - m(1,1) * m(2,3)) * invdet
      minv(3,1) = (m(2,1) * m(3,2) - m(3,1) * m(2,2)) * invdet
      minv(3,2) = (m(3,1) * m(1,2) - m(1,1) * m(3,2)) * invdet
      minv(3,3) = (m(1,1) * m(2,2) - m(2,1) * m(1,2)) * invdet

      return
      end




*     ******************************************
*     *                                        *
*     *         rho_symmetrizer_rho0           *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_rho0(use_translation,
     >                               n2ft3d,i_grid,f_grid,sf_grid,rho,
     >                               nx,ny,nz,n2ft3d_serial,
     >                               dn_serial,rho0,rho1)
      implicit none
      logical use_translation
      integer n2ft3d
      integer i_grid(4,n2ft3d)
      real*8  f_grid(3,n2ft3d)
      real*8  sf_grid(3,n2ft3d)
      real*8  rho(*)

      integer nx,ny,nz
      integer n2ft3d_serial
      real*8  dn_serial(nx+2,ny,nz)
      real*8  rho0(*),rho1(*)

*     **** local variables ****
      logical translation
      integer op,nops
      real*8  msym(3,4),xyz(3),scal1,minv(3,3),xyz2(3),xyz2m(3)


*     **** external functions ****
      integer  ion_sym_number_ops
      external ion_sym_number_ops
      integer  control_rtdb
      external control_rtdb
      real*8   lattice_unita
      external lattice_unita

c      character*72 cube_comment
c      character*50 filename
c
c      character*7 c_index_name
c      external    c_index_name


      xyz2(1) = 0.5d0*(lattice_unita(1,1) 
     >                +lattice_unita(1,2) 
     >                +lattice_unita(1,3))
      xyz2(2) = 0.5d0*(lattice_unita(2,1) 
     >                +lattice_unita(2,2) 
     >                +lattice_unita(2,3))
      xyz2(3) = 0.5d0*(lattice_unita(3,1) 
     >                +lattice_unita(3,2) 
     >                +lattice_unita(3,3))
      xyz2m(1) = -xyz2(1)
      xyz2m(2) = -xyz2(2)
      xyz2m(3) = -xyz2(3)

c      cube_comment = "rhosym "
c      filename = "rho.cube             "
c      call dplot_gcube_write(control_rtdb(),filename,-1,
c     >           cube_comment,rho)


      scal1 = 1.0d0/dble(nx*ny*nz)

      nops = ion_sym_number_ops()
      call Parallel_shared_vector_zero(.true.,n2ft3d_serial,dn_serial)
      if (use_translation) call D3dB_r_Copy(1,rho,rho0)
      call rho_symmetrizer_translate(xyz2,scal1,rho0)

      op   = 1
      do while (op.le.nops)
         call ion_sym_get_op(op,msym)


         call lattice_frac_to_r1(1,msym(1,4),xyz)
         xyz(1) = xyz(1) + xyz2m(1)
         xyz(2) = xyz(2) + xyz2m(2)
         xyz(3) = xyz(3) + xyz2m(3)
         translation = (dabs(xyz(1)).gt.1.0d-6).or.
     >                 (dabs(xyz(2)).gt.1.0d-6).or.
     >                 (dabs(xyz(3)).gt.1.0d-6)

         call dgemm('N','N',3,n2ft3d,3,1.0d0,
     >              msym,3,
     >              f_grid,3,
     >              0.0d0,sf_grid,3)
c         write(*,*) "translation=",translation
c         write(*,*) "xyz=",xyz
c         write(*,*) "msym="
c         write(*,*) msym(1,1),msym(1,2),msym(1,3)
c         write(*,*) msym(2,1),msym(2,2),msym(2,3)
c         write(*,*) msym(3,1),msym(3,2),msym(3,3)
        
                  
         if (use_translation) then
            call rho_symmetrizer_gen_serial0(.true.,n2ft3d,
     >                                    i_grid,sf_grid,rho0,
     >                                    nx,ny,nz,n2ft3d_serial,
     >                                    dn_serial)
            call rho_symmetrizer_to_rho(n2ft3d,n2ft3d_serial,nx,ny,nz,
     >                                  i_grid,
     >                                  dn_serial,rho1)
            if (translation) then
               call rho_symmetrizer_translate(xyz,scal1,rho1)
            end if

c             filename = "rhosym"//c_index_name(op)//".cube          "
c             call dplot_gcube_write(control_rtdb(),filename,-1,
c     >           cube_comment,rho1)

            call D3dB_rr_Sum2(1,rho1,rho)
         else
            call rho_symmetrizer_gen_serial0(.false.,n2ft3d,
     >                                    i_grid,sf_grid,rho,
     >                                    nx,ny,nz,n2ft3d_serial,
     >                                    dn_serial)
         end if

         !end if
         op = op + 1
      end do


      if (.not.use_translation) then
         call rho_symmetrizer_to_rho(n2ft3d,n2ft3d_serial,nx,ny,nz,
     >                                  i_grid,
     >                                  dn_serial,rho0)
         call D3dB_rr_Sum2(1,rho0,rho)
      end if

      call D3dB_r_SMul1(1,1.0d0/dble(nops+1),rho)

c      filename = "rhosym.cube          "
c      call dplot_gcube_write(control_rtdb(),filename,-1,
c     >           cube_comment,rho)


      return
      end 


*     ******************************************
*     *                                        *
*     *          rho_symmetrizer_dn            *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_dn(ispin,n2ft3d_in,dn)
      implicit none
      integer ispin,n2ft3d_in
      real*8  dn(n2ft3d_in,ispin)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** rhosymr_block common blocks ****
      integer nx,ny,nz,n2ft3d_serial,n2ft3d
      integer i_grid(2),f_grid(2),sf_grid(2),dn_serial(2),rho0(2)
      logical use_sym,use_translation
      common /rhosymr_block0/ i_grid,f_grid,sf_grid,dn_serial,rho0,
     >                        n2ft3d,n2ft3d_serial,nx,ny,nz,use_sym,
     >                        use_translation

      if (use_sym) then
         call rho_symmetrizer_rho0(use_translation,
     >                            n2ft3d,int_mb(i_grid(1)),
     >                            dbl_mb(f_grid(1)),
     >                            dbl_mb(sf_grid(1)),
     >                            dn,
     >                            nx,ny,nz,n2ft3d_serial,
     >                            dbl_mb(dn_serial(1)),
     >                            dbl_mb(rho0(1)),
     >                            dbl_mb(rho0(1)+n2ft3d))

         
         if (ispin.eq.2) then
            call rho_symmetrizer_rho0(use_translation,
     >                               n2ft3d,int_mb(i_grid(1)),
     >                               dbl_mb(f_grid(1)),
     >                               dbl_mb(sf_grid(1)),
     >                               dn(1,2),
     >                               nx,ny,nz,n2ft3d_serial,
     >                               dbl_mb(dn_serial(1)),
     >                               dbl_mb(rho0(1)),
     >                               dbl_mb(rho0(1)+n2ft3d))
         end if
      end if

      return
      end

*     ******************************************
*     *                                        *
*     *          rho_symmetrizer_print         *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrizer_print(unit)
      implicit none
      integer unit

*     **** rhosymr_block common blocks ****
      integer nx,ny,nz,n2ft3d_serial,n2ft3d
      integer i_grid(2),f_grid(2),sf_grid(2),dn_serial(2),rho0(2)
      logical use_sym,use_translation
      common /rhosymr_block0/ i_grid,f_grid,sf_grid,dn_serial,rho0,
     >                        n2ft3d,n2ft3d_serial,nx,ny,nz,use_sym,
     >                        use_translation

*     **** external functions ****
      character*10 ion_sym_group_name
      external     ion_sym_group_name
      integer      ion_sym_number_ops
      external     ion_sym_number_ops

      if (use_sym) then
         write(unit,*)
         write(unit,*) "symmetrizing density:"
         write(unit,'(6x,A,A10)') "Group name  ",ion_sym_group_name()
         write(unit,'(6x,A,I10)') "Group order ",ion_sym_number_ops()
      end if
      return
      end
