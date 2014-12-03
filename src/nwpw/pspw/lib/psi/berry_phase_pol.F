*
* $Id$
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Polarization calculated through Berry's Phase Approach
c  This is really simple since we have the Silvestrelli Overlap
c  Matrix routine from pspw_wannier. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine berry_phase_pol(psi,psir)
      implicit none
      integer psi(2),psir(2)
#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
      integer i,npack1,nfft3d,n2fft,nesq,idir,ms,ispin,ne(2)
      integer pivs(2),smat(2),nemax,info,n,ik,j
      real*8 pol(3,3),tx,ty,tz,ex,ey,ez,dx,dy,dz,tmp
      real*8 gs,g(3),dv
      complex*16 det
      logical val
ccccc external ccccccccccccccccccccccccccc
      integer ion_katm,ion_nion
      real*8 ion_rion,psp_zv,lattice_unitg,lattice_omega
      external ion_katm,ion_nion
      external ion_rion,psp_zv,lattice_unitg,lattice_omega
cccccccccccccccccccccccccccccccccccccccccc
c      VAL=BA_set_auto_verify(.TRUE.)
ccccccccccccccccccccccccccccccccccccccccccc
      call Pack_npack(1,npack1)
      call D3db_nfft3d(1,nfft3d)
      call psi_get_ne(ispin,ne)
      n2fft=nfft3d*2
      nesq=ne(1)*ne(1)*2
      val=BA_push_get(mt_int,ne(1),'pivs',pivs(2),pivs(1))
      val=val.and.BA_push_get(mt_dcpl,nesq,'smat',smat(2),smat(1))
      if (.not.val) then
        call errquit("berry_phase_pol: stack depleted",0,0)
      end if
      nemax = ne(1)+ne(2)
      do n=1,nemax
         call Pack_c_Copy(1,dcpl_mb(psi(1)+(n-1)*npack1),
     >                       dbl_mb(psir(1)+(n-1)*n2fft))
      end do
      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psir(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psir(1)))
      do idir=1,3
         det=dcmplx(1.0d0,0.0d0)
         g(1)=lattice_unitg(1,idir)
         g(2)=lattice_unitg(2,idir)
         g(3)=lattice_unitg(3,idir) 
         gs=(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
         if (dabs(gs).lt.1.d-12) then
           call errquit("g too small!",0,0)
         end if
         gs=1.d0/gs
         do ms=1,ispin
            call silvestrelli_overlap(g,ne(ms),
     &          dbl_mb( psir(1)+(ms-1)*n2fft*ne(1) ),
     &          dcpl_mb(smat(1)))
            call zgetrf(ne(ms),ne(ms),dcpl_mb(smat(1)),
     &         ne(ms),int_mb(pivs(1)),info)
            if (info.ne.0) then
              call errquit("berry_phase_pol:failed",0,0)
            end if
            do i=1,ne(ms)
              det=det*dcpl_mb(smat(1)+(i-1)*ne(ms)+i-1)
            end do
         end do
         tmp= -datan2(dimag(det),dble(det))
         if (ispin.eq.1) then 
           tmp=tmp*2.0d0
         end if 
         pol(1,idir)= tmp*g(1)*gs
         pol(2,idir)= tmp*g(2)*gs
         pol(3,idir)= tmp*g(3)*gs 
      end do
      dx=0.0d0
      dy=0.0d0
      dz=0.0d0
      do i=1,ion_nion()
        ik=ion_katm(i)
        dx=dx+psp_zv(ik)*ion_rion(1,i)
        dy=dy+psp_zv(ik)*ion_rion(2,i)
        dz=dz+psp_zv(ik)*ion_rion(3,i)
      end do
      ex=0.0d0
      ey=0.0d0
      ez=0.0d0
      do i=1,3
         ex=ex+pol(1,i)
         ey=ey+pol(2,i)
         ez=ez+pol(3,i)
      end do
      tx=dx-ex
      ty=dy-ey
      tz=dz-ez
c      if (nodeid.eq.MASTER) then
      write(6,*)"--------------------------------"
      write(6,*)"Berry Phase Results for the"
      write(6,*)"Components of the Polarization Tensor"
      write(6,1000)pol(1,1),pol(2,1),pol(3,1) 
      write(6,1000)pol(1,2),pol(2,2),pol(3,2) 
      write(6,1000)pol(1,3),pol(2,3),pol(3,3) 
      write(6,*)"--------------------------------"
      write(6,*)"Electronic DelP = P Omega in"
      write(6,*)"Cartesian Cooridinates"
      write(6,*)"--------------------------------"
      write(6,*)"Nuclear"
      write(6,1000)dx,dy,dz
      tmp=dsqrt(dx*dx+dy*dy+dz*dz)
      write(6,1010)tmp
      write(6,*)"Electronic"
      write(6,1000)ex,ey,ez
      tmp=dsqrt(ex*ex+ey*ey+ez*ez)
      write(6,1010)tmp
      write(6,*)"Total"
      write(6,1000)tx,ty,tz
      tmp=dsqrt(tx*tx+ty*ty+tz*tz)
      write(6,1010)tmp
      call flush(6)
c      end if
      val=BA_pop_stack(smat(2))
      val=val.and.BA_pop_stack(pivs(2))
      if (.not.val) then
        call errquit("berry_phase_pol: pop stack failed!",0,0)
      end if
1000  format("x= ",d15.6," y= ",d15.6," z= ",d15.6)
1010  format("|mu|= ",d15.6)
      return
      end
ccccccccccccc
 
             
             
