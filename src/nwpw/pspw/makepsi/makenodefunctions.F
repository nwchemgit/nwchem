
      subroutine psi_makebasisfunctions(nnodes,
     >                                 ispin,ne,neq,
     >                                 Rw,psiin,psiout)
      implicit none
      integer nnodes,ispin,ne(2),neq(2)
      real*8 Rw(3,*)
      complex*16 psiin(*),psiout(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer k1,k2,k3,nbasis,n2ft3d,ftype,itype
      integer tmp1(2),tmp2(2),rgrid(2)
      character*255 filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name
      integer  Pack_G_indx
      external Pack_G_indx

      call D3dB_n2ft3d(1,n2ft3d)

      value = BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1))
      if (.not. value)
     >   call errquit('pspw_makebasisfunctions:out of stack',0,0)

      call lattice_r_grid(dbl_mb(rgrid(1)))

      nbasis = 0
      do k3=0,nnodes
      do k2=0,nnodes-k3
      do k1=0,nnodes-k3-k2
      !do k1=0,8

         itype = k1 + k2*100 + k3*10000
         call psi_makenodefunctions(0,itype,
     >                              ispin,neq,
     >                              Rw,psiin,psiout,
     >                              dbl_mb(rgrid(1)),
     >                              dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),
     >                              dbl_mb(Pack_G_indx(1,1)),
     >                              dbl_mb(Pack_G_indx(1,2)),
     >                              dbl_mb(Pack_G_indx(1,3)))
         call util_file_prefix(c_index_name(nbasis)//'.mbbasis',
     >                         filename)
         call psi_write_filename(filename,ispin,ne,psiout)
         nbasis = nbasis + 1

      end do
      end do
      end do

      value =           MA_pop_stack(rgrid(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value)
     >   call errquit('pspw_makebasisfunctions:popping stack',0,0)

      return
      end

*     *******************************************
*     *                                         *
*     *            psi_makenodefunctions        *
*     *                                         *
*     *******************************************
      subroutine psi_makenodefunctions(ftype,itype,
     >                                 ispin,neq,
     >                                 Rw,psiin,psiout,
     >                                 rgrid,psir,nodefunction,
     >                                 gx,gy,gz)
      implicit none
      integer ftype,itype
      integer ispin,neq(2)
      real*8  Rw(3,*)
      complex*16 psiin(*)
      complex*16 psiout(*)
      real*8     rgrid(3,*)
      real*8     psir(*)
      real*8     nodefunction(*)
      real*8     gx(*),gy(*),gz(*)

*     **** local variables ****
      integer q,n,nx,ny,nz,shift1,npack1,n2ft3d
      real*8  GR,scal1


      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)


      do q=1,(neq(1)+neq(2))
         call Dneall_qton(q,n)
         shift1 = 1+(q-1)*npack1


         if (ftype.eq.0) then
            call Pack_c_Copy(1,psiin(shift1),psir)
            call psi_derivative_f(itype,gx,gy,gz,psir)
         else
*           *** generate node function ***
            call psi_makenode_f(ftype,itype,n2ft3d,Rw(1,n),rgrid,
     >                          nodefunction)

*           *** generate psiout(r) = nodefunction(r)*psiin(r) ***
            call Pack_c_Copy(1,psiin(shift1),psir)
            call Pack_c_unpack(1,psir)
            call D3dB_cr_fft3b(1,psir)

            call D3dB_rr_Mul2(1,nodefunction,psir)

            call D3dB_r_SMul1(1,scal1,psir)
            call D3dB_rc_fft3f(1,psir)
            call Pack_c_pack(1,psir)
         end if
         call Pack_cc_dot(1,psir,psir,GR)
         GR = 1.0d0/dsqrt(GR)
         call Pack_c_SMul1(1,Gr,psir)
         call Pack_c_copy(1,psir,psiout(shift1))
      end do

      return
      end 

      subroutine psi_derivative_f(i,gx,gy,gz,f)
      implicit none
      integer    i
      real*8     gx(*),gy(*),gz(*)
      complex*16 f(*)
      integer k1,k2,k3,j

      k3 = i/10000
      k2 = (i - 10000*k3)/100
      k1 = (i - 10000*k3 - 100*k2)

      do j=1,k1
         call Pack_tc_iMul2(1,gx,f)
      end do
      do j=1,k2
         call Pack_tc_iMul2(1,gy,f)
      end do
      do j=1,k3
         call Pack_tc_iMul2(1,gz,f)
      end do

      return
      end 
 
      subroutine psi_makenode_f(ftype,i,n2ft3d,Rw,rgrid,f)
      implicit none
      integer ftype,i,n2ft3d
      real*8  Rw(3),rgrid(3,*)
      real*8 f(*)

*     **** local variables ****
      integer k1,k2,k3,k
      real*8  gx,gy,gz,x,y,z,fx,fy,fz,pi,r
      real*8  rxyz(3),frac(3)

      pi = 4.0d0*datan(1.0d0)

*     **** sine functions ****
      if (ftype.eq.1) then
         k3 = i/10000
         k2 = (i - 10000*k3)/100
         k1 = (i - 10000*k3 - 100*k2)
         k1= k1 + 1
         k2= k2 + 1
         k3= k3 + 1
         do k=1,n2ft3d
            rxyz(1) = rgrid(1,k)-Rw(1)
            rxyz(2) = rgrid(2,k)-Rw(2)
            rxyz(3) = rgrid(3,k)-Rw(3)
            call lattice_r1_to_frac(1,rxyz,frac)
            frac(1) = frac(1) - DNINT(frac(1))
            frac(2) = frac(2) - DNINT(frac(2))
            frac(3) = frac(3) - DNINT(frac(3))
            x = k1*pi*frac(1)
            y = k2*pi*frac(2)
            z = k3*pi*frac(3)
            fx = dsin(x)
            fy = dsin(y)
            fz = dsin(z)
            if (k1.le.1) fx = 1.0d0
            if (k2.le.1) fy = 1.0d0
            if (k3.le.1) fz = 1.0d0
            f(k) = fx*fy*fz
         end do

      else if (ftype.eq.2) then
         do k=1,n2ft3d
            rxyz(1) = rgrid(1,k)-Rw(1)
            rxyz(2) = rgrid(2,k)-Rw(2)
            rxyz(3) = rgrid(3,k)-Rw(3)
            call lattice_r1_to_frac(1,rxyz,frac)
            frac(1) = frac(1) - DNINT(frac(1))
            frac(2) = frac(2) - DNINT(frac(2))
            frac(3) = frac(3) - DNINT(frac(3))
            x = frac(1) 
            y = frac(2) 
            z = frac(3) 
            !x = rgrid(1,k)-Rw(1)
            !y = rgrid(2,k)-Rw(2)
            !z = rgrid(3,k)-Rw(3)
            r = dsqrt(x*x + y*y + z*z)
            if (r.lt.1.0d-9) r = 1.0d-9
            x = x/r
            y = y/r
            z = z/r
            if (i.eq.0) then 
               f(k) = 1.0d0
            else if (i.eq.1) then 
               f(k) = x
            else if (i.eq.2) then 
               f(k) = y
            else if (i.eq.3) then 
               f(k) = z
            else if (i.eq.4) then
               f(k) = x*y
            else if (i.eq.5) then
               f(k) = y*z
            else if (i.eq.6) then
               f(k) = 3.0d0*z**2 - 1.0d0
            else if (i.eq.7) then
               f(k) = z*x
            else if (i.eq.8) then
               f(k) = x**2 - y**2
            else if (i.eq.9) then
               f(k) = y*(3.0d0*(1.0d0-z**2)-4.0d0*y**2)
            else if (i.eq.10) then
               f(k) = x*y*z
            else if (i.eq.11) then
               f(k) = y*(5.0d0*z**2-1.0d0)
            else if (i.eq.12) then
               f(k) = z*(5.0d0*z**2-3.0d0)
            else if (i.eq.13) then
               f(k) = x*(5.0d0*z**2-1.0d0)
            else if (i.eq.14) then
               f(k) = z*(x**2 - y**2)
            else if (i.eq.15) then
               f(k) = x*(4.0d0*x**2 - 3.0d0*(1.0d0-z**2))
            else 
               f(k) = 0.0d0
            end if
         end do
      end if
   
      return
      end
       
