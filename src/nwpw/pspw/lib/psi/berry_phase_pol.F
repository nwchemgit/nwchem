cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Polarization calculated through Berry's Phase Approach
c  This is really simple since we have the Silvestrelli Overlap
c  Matrix routine from pspw_wannier. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine berry_phase_pol(psi,psir)
      implicit none
      integer psi(2),psir(2)
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
      integer i,npack1,nfft3d,n2fft,nesq,idir,ms,ispin,ne(2)
      integer pivs(2),smat(2),nemax,info,n
      real*8 pol(3),tx,ty,tz,ex,ey,ez,dx,dy,dz
      real*8 gs,g(3)
      complex*16 det
      logical val
ccccc external ccccccccccccccccccccccccccc
      integer ion_katm,ion_nion
      real*8 ion_rion,psp_zv,lattice_unitg
      external ion_katm,ion_nion
      external ion_rion,psp_zv,lattice_unitg
cccccccccccccccccccccccccccccccccccccccccc

      VAL=MA_SET_AUTO_VERIFY(.TRUE.)
ccccccccccccccccccccccccccccccccccccccccccc
      write(*,*)"t0o"
      call flush(6)
      call Pack_npack(1,npack1)
      write(*,*)"t0a"
      call flush(6)
      call D3db_nfft3d(1,nfft3d)
      call psi_get_ne(ispin,ne)
      write(*,*)"t0b"
      call flush(6)
      n2fft=nfft3d*2
      nesq=ne(1)*ne(1)*2
      val=MA_push_get(mt_int,ne(1),'pivs',pivs(2),pivs(1))
      val=val.and.MA_push_get(mt_dcpl,nesq,'smat',smat(2),smat(1))
      if (.not.val) then
        call errquit("berry_phase_pol: stack depleted",0,0)
      end if
      write(*,*)"t0c",ne(1),ne(2)
      call flush(6)

      nemax = ne(1)+ne(2)
      do n=1,nemax
         write(*,*)"t0cc",n
         call flush(6)
         call Pack_c_Copy(1,dcpl_mb(psi(1)+(n-1)*npack1),
     >                       dbl_mb(psir(1)+(n-1)*n2fft))
      end do
      write(*,*)"t0d"
      call flush(6)
      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psir(1)))
      write(*,*)"t0e"
      call flush(6)
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psir(1)))
      write(*,*)"t0f"
      call flush(6)
CCCCCCCCCCCCC
      VAL=MA_VERIFY_ALLOCATOR_STUFF() 
      WRITE(6,*)"T1"
      CALL FLUSH(6)
CCCCCCCCCCCCC
      do idir=1,3
         det=dcmplx(1.0d0,0.0d0)
         g(1)=lattice_unitg(1,idir)
         g(2)=lattice_unitg(2,idir)
         g(3)=lattice_unitg(3,idir)      
         gs=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
         if (dabs(gs).lt.1.d-12) then
           call errquit("g too small!",0,0)
         end if
         gs=1.d0/gs
CCCCCCCCCCCCC
      VAL=MA_VERIFY_ALLOCATOR_STUFF() 
      WRITE(6,*)"T2"
      CALL FLUSH(6)
CCCCCCCCCCCCC
         do ms=1,ispin
            call silvestrelli_overlap(g,ne(ms),
     &          dbl_mb( psir(1)+(ms-1)*n2fft*ne(1) ),
     &          dcpl_mb(smat(1)))
CCCCCCCCCCCCC
      WRITE(6,*)"T3"
      VAL=MA_VERIFY_ALLOCATOR_STUFF() 
      CALL FLUSH(6)
CCCCCCCCCCCCC
            call zgetrf(ne(ms),ne(ms),dcpl_mb(smat(1)),
     &         ne(ms),pivs,info)
CCCCCCCCCCCCC
      VAL=MA_VERIFY_ALLOCATOR_STUFF() 
      WRITE(6,*)"T4"
      CALL FLUSH(6)
CCCCCCCCCCCCC
            if (info.ne.0) then
              call errquit("berry_phase_pol:failed");
            end if
            do i=1,ne(ms)
              det=det*dcpl_mb(smat(1)+(i-1)*ne(ms)+i-1)
            end do
         end do
         pol(idir)= -datan2(dimag(det),dble(det))*gs  
      end do
CCCCCCCCCCCCC
      VAL=MA_VERIFY_ALLOCATOR_STUFF() 
      WRITE(6,*)"T5"
      CALL FLUSH(6)
CCCCCCCCCCCCC
      dx=0.0d0
      dy=0.0d0
      dz=0.0d0
      do i=1,ion_nion()
        dx=dx+psp_zv(ion_katm(i))*ion_rion(i,1)
        dy=dy+psp_zv(ion_katm(i))*ion_rion(i,2)
        dz=dz+psp_zv(ion_katm(i))*ion_rion(i,3)
      end do
      ex=0.0d0
      ey=0.0d0
      ez=0.0d0
      do i=1,3
         g(1)=lattice_unitg(1,i)
         g(2)=lattice_unitg(2,i)
         g(3)=lattice_unitg(3,i)      
         gs=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
         if (dabs(gs).lt.1.d-12) then
           call errquit("g too small!",0,0)
         end if
         gs=1.d0/gs
         g(1)=g(1)*gs
         g(2)=g(2)*gs
         g(3)=g(3)*gs
        ex=ex-pol(i)*g(1)
        ey=ey-pol(i)*g(2)
        ez=ez-pol(i)*g(3)
      end do
      tx=dx+ex
      ty=dy+ey
      tz=tz+ez
      write(6,*)"Berry Phase Results for dipole (P*Omega)"
      write(6,*)"Electronic Along Unit Cell vectors"
      write(6,1000)ex,ey,ez
      write(6,*)"Cartesian Cooridinates"
      write(6,*)"Nuclear"
      write(6,1000)dx,dy,dz
      write(6,*)"Electronic"
      write(6,1000)ex,ey,ez
      write(6,*)"Total"
      write(6,1000)tx,ty,tz
      call flush(6)
      val=MA_pop_stack(smat(2))
      val=val.and.MA_pop_stack(pivs(2))
      if (.not.val) then
        call errquit("berry_phase_pol: pop stack failed!",0,0)
      end if
1000  format("x= ",d15.6," y= ",d15.6," z= ",d15.6)
      return
      end
ccccccccccccc
 
             