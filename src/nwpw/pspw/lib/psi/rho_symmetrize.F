
*     ******************************
*     *                            *
*     *     rho_symmetrize_init    *
*     *                            *
*     ******************************
      subroutine rho_symmetrize_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer p,indx
      integer k1,k2,k3,i,j,k,nyh,nzh
      integer G(3),iG(3),g_itmp(2)

*     **** rhosym_block common blocks ****
      integer nx,ny,nz,nxh,nfft3d_serial,nfft3d
      integer g_indexes(2),dng_serial(2),dng_serial0(2)
      real*8 scal1
      logical use_sym
      common /rhosym_block0/ g_indexes,dng_serial,dng_serial0,
     >                     nfft3d,nfft3d_serial,nx,ny,nz,nxh
      common /rhosym_block1/ scal1,use_sym

*     **** external functions ****
      logical  control_rho_use_symmetry
      external control_rho_use_symmetry
      integer  ion_sym_number_ops
      external ion_sym_number_ops

      use_sym = control_rho_use_symmetry()
     >         .and.(ion_sym_number_ops().gt.0)

      write(*,*) "ion_sym_number_ops=",ion_sym_number_ops(),use_sym
      if (use_sym) then
         call Parallel2d_taskid_i(taskid)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         nxh = nx/2
         nyh = ny/2
         nzh = nz/2
         nfft3d_serial = (nxh+1)*ny*nz
         scal1 = 1.0d0/dble(nx*ny*nz)

         call D3dB_nfft3d(1,nfft3d)
         !call Pack_npack(0,npack0)

         !write(*,*) "NPACK0=",npack0
         write(*,*) "NFFT3D=",nfft3d
*        **** allocate from heap for rhosym_block common block ****
         value = BA_alloc_get(mt_dbl,(3*nfft3d),'f_indx',
     >                        f_indx(2),f_indx(1))
         value = value.and.
     >           BA_alloc_get(mt_dcpl,nfft3d_serial,'dng_serial',
     >                        dng_serial(2),dng_serial(1))
         value = value.and.
     >           BA_alloc_get(mt_dcpl,nfft3d_serial,'dng_serial0',
     >                        dng_serial0(2),dng_serial0(1))
         if (.not.value) 
     >      call errquit('rho_symmetrize_init: out of heap',0, MA_ERR)

*        **** allocate stack ****
         value = BA_push_get(mt_int,(3*nfft3d),'g_itmp',
     >                       g_itmp(2),g_itmp(1))
         if (.not.value) 
     >      call errquit('rho_symmetrize_init: out of stack',1, MA_ERR)

         iG(1) = g_indexes(1)
         iG(2) = iG(1) + nfft3d
         iG(3) = iG(2) + nfft3d

         G(1) = g_itmp(1)
         G(2) = G(1) + nfft3d
         G(3) = G(2) + nfft3d

         call Parallel2d_taskid_i(taskid)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         nxh = nx/2
         nyh = ny/2
         nzh = nz/2

         call icopy(3*nfft3d,0,0,int_mb(g_indexes(1)),1)
         do k3 = -nzh+1, nzh
            do k2 = -nyh+1, nyh
               do k1 = 0,nxh
                   i=k1
                   j=k2
                   k=k3
                   if (i .lt. 0) i = i + nx
                   if (j .lt. 0) j = j + ny
                   if (k .lt. 0) k = k + nz

                   call D3dB_ijktoindexp(1,i+1,j+1,k+1,indx,p)
                   if (p .eq. taskid) then
                      !int_mb(iG(1)+indx-1) = k1
                      !int_mb(iG(2)+indx-1) = k2
                      !int_mb(iG(3)+indx-1) = k3
                      int_mb(f_indx(1)+3*(indx-1))   = i/dble(nx)
                      int_mb(f_indx(2)+3*(indx-1)+1) = j/dble(ny)
                      int_mb(f_indx(3)+3*(indx-1)+2) = k/dble(nz)
                   end if
               end do
            end do
         end do
         !call Pack_i_pack(0,int_mb(G(1)))
         !call Pack_i_pack(0,int_mb(G(2)))
         !call Pack_i_pack(0,int_mb(G(3)))
         !call icopy(npack0,int_mb(G(1)),1,int_mb(iG(1)),1)
         !call icopy(npack0,int_mb(G(2)),1,int_mb(iG(2)),1)
         !call icopy(npack0,int_mb(G(3)),1,int_mb(iG(3)),1)

         value = BA_pop_stack(g_itmp(2))
         if (.not.value) 
     >      call errquit('rho_symmetrize_init:popping stack',2,MA_ERR)

      end if

      return
      end


*     ******************************
*     *                            *
*     *     rho_symmetrize_end     *
*     *                            *
*     ******************************
      subroutine rho_symmetrize_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

*     **** rhosym_block common blocks ****
      integer nx,ny,nz,nxh,nfft3d_serial,nfft3d
      integer g_indexes(2),dng_serial(2),dng_serial0(2)
      real*8 scal1
      logical use_sym
      common /rhosym_block0/ g_indexes,dng_serial,dng_serial0,
     >                     nfft3d,nfft3d_serial,nx,ny,nz,nxh
      common /rhosym_block1/ scal1,use_sym


      if (use_sym) then
         value = BA_free_heap(dng_serial(2))
         value = value.and.BA_free_heap(dng_serial0(2))
         value = value.and.BA_free_heap(g_indexes(2))
         if (.not.value) 
     >      call errquit('rho_symmetrize_end: freeing heap',0, MA_ERR)
      end if
      return
      end

*     ******************************************
*     *                                        *
*     *         rho_symmetrize_inv33           *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrize_inv33(m,minv)
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
*     *         rho_symmetrize_dng0            *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrize_dng0(nfft3d,ig,g,dng,
     >                               nxh1,nx,ny,nz,nfft3d_serial,
     >                               dng_serial,dng_serial0)
      implicit none
      integer    nfft3d
      integer    ig(nfft3d,3)
      real*8     g(nfft3d,3)
      complex*16 dng(*)

      integer    nxh1,nx,ny,nz
      integer    nfft3d_serial
      complex*16 dng_serial(nxh1,ny,nz)
      complex*16 dng_serial0(nxh1,ny,nz)

*     **** local variables ****
      logical translation,toobig
      integer op,nops,i,j,k
      integer ig1,ig2,ig3,ih1,ih2,ih3,nxh,nyh,nzh,t1,t2,t3,s1,s2,s3
      real*8  msym(3,4),xyz(3),gr,g1,g2,g3,h1,h2,h3,minv(3,3)
      complex*16 sdng,tdng,terr
      integer igg1,igg2,igg3,ihh1,ihh2,ihh3
      integer igg(3)

*     **** external functions ****
      integer  ion_sym_number_ops
      external ion_sym_number_ops


      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nops = ion_sym_number_ops()
      !nops = 1
      !write(*,*) "pi=",pi
      !write(*,*) "GERA,nops=",nops,nfft3d_serial

!$OMP DO
      do k=1,nfft3d
         ig1 = ig(k,1)
         ig2 = ig(k,2)
         ig3 = ig(k,3)
         dng_serial0(ig1+1,ig2+1,ig3+1) = dng(k)
      end do
!$OMP END DO
      call Parallel_Vector_SumAll(2*nfft3d_serial,dng_serial0)


      call Parallel_shared_vector_zero(.true.,2*nfft3d_serial,
     >                                 dng_serial)

      op = 1
      nops = 1
      do while (op.le.nops)
         call ion_sym_get_op(op,msym)
         translation = (dabs(msym(1,4)).gt.1.0d-6).or.
     >                 (dabs(msym(2,4)).gt.1.0d-6).or.
     >                 (dabs(msym(3,4)).gt.1.0d-6)
         call lattice_frac_to_r1(1,msym(1,4),xyz)

        write(*,*) "msym, op=",op
        write(*,*) msym(1,1),msym(1,2),msym(1,3)
        write(*,*) msym(2,1),msym(2,2),msym(2,3)
        write(*,*) msym(3,1),msym(3,2),msym(3,3)
        write(*,*) "xyz=",xyz

        op = op + 1
                  
!$OMP DO
         do k=1,nfft3d
            ig1 = ig(k,1)
            ig2 = ig(k,2)
            ig3 = ig(k,3)
            g1 = ig1*msym(1,1)+ig2*msym(2,1)+ig3*msym(3,1)
            g2 = ig1*msym(1,2)+ig2*msym(2,2)+ig3*msym(3,2)
            g3 = ig1*msym(1,3)+ig2*msym(2,3)+ig3*msym(3,3)
            ig1 = NINT(g1)
            ig2 = NINT(g2)
            ig3 = NINT(g3)

            do while (ig1.ge.nx)
               ig1 = ig1 - nx
            end do
            do while (ig1.lt.0)
               ig1 = ig1 + nx
            end do

            do while (ig2.ge.ny)
               ig2 = ig2 - ny
            end do
            do while (ig2.lt.0)
               ig2 = ig2 + ny
            end do

            do while (ig3.ge.nz)
               ig3 = ig3 - nz
            end do
            do while (ig3.lt.0)
               ig3 = ig3 + nz
            end do

            if (ig1.gt.nxh) ig1 = ig1 - nx
            if (ig2.gt.nyh) ig2 = ig2 - ny
            if (ig3.gt.nzh) ig3 = ig3 - nz

            sdng = dng(k)
            if (ig1.lt.0) then
               sdng = dconjg(sdng)
               ig1 = -ig1
               ig2 = -ig2
               ig3 = -ig3
            end if

            if (translation) then
               gr = g(k,1)*xyz(1)+g(k,2)*xyz(2)+g(k,3)*xyz(3)
               sdng = sdng*dcmplx(dcos(gr),dsin(gr))
            end if

            if (ig1.lt.0) ig1 = ig1 + nx
            if (ig2.lt.0) ig2 = ig2 + ny
            if (ig3.lt.0) ig3 = ig3 + nz

            write(*,*) "k,dngs=",k,ig(k,1),ig(k,2),ig(k,3), " ==>",
     >                 ig1,ig2,ig3
            write(*,*) dng_serial0(ig1+1,ig2+1,ig3+1),sdng, 
     >                 dng_serial0(ig1+1,ig2+1,ig3+1)-sdng
            write(*,*)
            dng_serial(ig1+1,ig2+1,ig3+1) 
     >         = dng_serial(ig1+1,ig2+1,ig3+1) + sdng

         end do
!$OMP END DO

      end do

      call Parallel_Vector_SumAll(2*nfft3d_serial,dng_serial)

!$OMP DO
      do k=1,nfft3d
         ig1 = ig(k,1)
         ig2 = ig(k,2)
         ig3 = ig(k,3)
         if (ig1.lt.0) ig1 = ig1 + nx
         if (ig2.lt.0) ig2 = ig2 + ny
         if (ig3.lt.0) ig3 = ig3 + nz
         dng(k) = dng(k) + dng_serial(ig1+1,ig2+1,ig3+1)
         !dng(k) = dng_serial(ig1+1,ig2+1,ig3+1)
      end do
!$OMP END DO
      !call Pack_c_SMul1(0,1.0d0/dble(nops+1),dng)
      call D3dB_c_SMul1(1,1.0d0/dble(nops+1),dng)

      return
      end 


*     ******************************************
*     *                                        *
*     *          rho_symmetrize_dng            *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrize_dng(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"
      
*     **** rhosym_block common blocks ****
      integer nx,ny,nz,nxh,nfft3d_serial,nfft3d
      integer g_indexes(2),dng_serial(2),dng_serial0(2)
      real*8 scal1
      logical use_sym
      common /rhosym_block0/ g_indexes,dng_serial,dng_serial0,
     >                     nfft3d,nfft3d_serial,nx,ny,nz,nxh
      common /rhosym_block1/ scal1,use_sym

*     **** external functions ****
      integer  G_indx
      external G_indx

      if (use_sym) then
         call rho_symmetrize_dng0(nfft3d,
     >                            int_mb(g_indexes(1)),
     >                            dbl_mb(G_indx(1)),
     >                            dng,
     >                            nxh+1,nx,ny,nz,nfft3d_serial,
     >                            dcpl_mb(dng_serial(1)),
     >                            dcpl_mb(dng_serial0(1)))
      end if

      return
      end


*     ******************************************
*     *                                        *
*     *           rho_symmetrize_dn            *
*     *                                        *
*     ******************************************
      subroutine rho_symmetrize_dn(ispin,n2ft3d,dn,tmp)
      implicit none
      integer ispin,n2ft3d
      real*8  dn(n2ft3d,ispin)
      real*8  tmp(n2ft3d)

*     **** rhosym_block1 common block ****
      integer nx,ny,nz,nxh,nfft3d_serial,nfft3d
      integer g_indexes(2),dng_serial(2),dng_serial0(2)
      real*8 scal1,trans(3)
      logical use_sym
      common /rhosym_block0/ g_indexes,dng_serial,dng_serial0,
     >                     nfft3d,nfft3d_serial,nx,ny,nz,nxh
      common /rhosym_block1/ scal1,use_sym

      integer i
      real*8 sumall,sumall0
      real*8  lattice_omega
      external lattice_omega


      if (use_sym) then
         if (ispin.eq.1) then
            call D3dB_r_Copy(1,dn,tmp)
            call D3dB_r_dsum(1,dn,sumall0)
            !call D3dB_rr_Sqrt1(1,dn)
            call D3dB_r_SMul1(1,scal1,dn)
            call D3dB_rc_fft3f(1,dn)

            !call dng_shift_center(dn)
            call rho_symmetrize_dng(dn)
            !call dng_unshift_center(dn)

            !call Pack_c_unpack(0,dn)
            call D3dB_cr_fft3b(1,dn)

            do i=1,n2ft3d
               write(*,*) "new dn=",i,dn(i,1),tmp(i),dn(i,1)-tmp(i)
            end do

            call D3dB_r_dsum(1,dn,sumall)
            write(*,*) "sumall=",sumall0*lattice_omega()*scal1,
     >                           sumall *lattice_omega()*scal1
c            call D3dB_r_SMul1(1,sumall0/sumall,dn)
         else
            call D3dB_rr_Minus(1,dn(1,1),dn(1,ispin),tmp)
            call D3dB_rr_Sum2(1,dn(1,ispin),dn)
            call D3dB_r_SMul1(1,scal1,dn)
            call D3dB_r_SMul1(1,scal1,tmp)
            call D3dB_rc_fft3f(1,dn)
            call D3dB_rc_fft3f(1,tmp)
            !call Pack_c_pack(0,dn)
            !call Pack_c_pack(0,tmp)

            call rho_symmetrize_dng(dn)
            call rho_symmetrize_dng(tmp)

            !call Pack_c_unpack(0,dn)
            !call Pack_c_unpack(0,tmp)
            call D3dB_cr_fft3b(1,dn)
            call D3dB_cr_fft3b(1,tmp)
            call D3dB_rr_Minus(1,dn,tmp,dn(1,2))
            call D3dB_rr_Sum2(1,tmp,dn)
            call D3dB_r_SMul1(1,0.5d0,dn(1,1))
            call D3dB_r_SMul1(1,0.5d0,dn(1,2))
         end if
      end if

      return
      end


      subroutine dng_shift_center(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer exi(2),nfft3d
      real*8 trans(3)

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) call errquit('dng_shift_center:push stack',0,0)

      trans(1) = (   lattice_unita(1,1)
     >             + lattice_unita(1,2)
     >             + lattice_unita(1,3) )/2.0d0
      trans(2) = (   lattice_unita(2,1)
     >             + lattice_unita(2,2)
     >             + lattice_unita(2,3) )/2.0d0
      trans(3) = (   lattice_unita(3,1)
     >             + lattice_unita(3,2)
     >             + lattice_unita(3,3) )/2.0d0
      call strfac_translate(trans,dcpl_mb(exi(1)))


      call D3dB_cc_Mul2(1,dcpl_mb(exi(1)),dng)

      value = BA_pop_stack(exi(2))
      if (.not. value) call errquit('dng_shift_center:pop stack',1,0)
      return
      end


      subroutine dng_unshift_center(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer exi(2),nfft3d
      real*8 trans(3)
      
*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) call errquit('dng_unshift_center:push stack',0,0)
      
      trans(1) = - (   lattice_unita(1,1)
     >               + lattice_unita(1,2)
     >               + lattice_unita(1,3) )/2.0d0
      trans(2) = - (   lattice_unita(2,1)
     >               + lattice_unita(2,2)
     >               + lattice_unita(2,3) )/2.0d0
      trans(3) = - (   lattice_unita(3,1)
     >               + lattice_unita(3,2)
     >               + lattice_unita(3,3) )/2.0d0
      call strfac_translate(trans,dcpl_mb(exi(1)))
     
     
      call D3dB_cc_Mul2(1,dcpl_mb(exi(1)),dng)
     
      value = BA_pop_stack(exi(2))
      if (.not. value) call errquit('dng_unshift_center:pop stack',1,0)
      return
      end
