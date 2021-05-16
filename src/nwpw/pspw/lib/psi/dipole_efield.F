
*     ******************************
*     *                            *
*     *     dipole_Efield_init     *
*     *                            *
*     ******************************
      subroutine dipole_Efield_init(ispin,ne,efield0,efield0_center)
      implicit none
      integer ispin,ne(2)
      real*8  efield0(3),efield0_center(3)

#include "bafdecls.fh"
#include "errquit.fh"

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

c     *** local variables ***
      logical value
      integer n2ft3d

      ispin0 = ispin
      ne0(1) = ne(1)
      ne0(2) = ne(2)

      efield(1) = efield0(1)
      efield(2) = efield0(2)
      efield(3) = efield0(3)

      efield_center(1) = efield0_center(1)
      efield_center(2) = efield0_center(2)
      efield_center(3) = efield0_center(3)

      call dipole_Efield_set_ion()
      call dipole_Efield_set_bv()

c     **** allocate from heap ****
      call D3dB_n2ft3d(1,n2ft3d)
      value =           BA_alloc_get(mt_dbl,n2ft3d,'Tc',Tc(2),Tc(1))
      value = value.and.BA_alloc_get(mt_dbl,n2ft3d,'Ts',Ts(2),Ts(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*n2ft3d),'rgrid',rgrid(2),rgrid(1))
      if (.not.value) 
     >   call errquit('dipole_efield_init:out of memory',0,MA_ERR)

c     **** generate rgrid ****
      call lattice_r_grid(dbl_mb(rgrid(1)))

      return
      end 

*     ******************************
*     *                            *
*     *     dipole_Efield_end      *
*     *                            *
*     ******************************
      subroutine dipole_Efield_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

c     *** local variables ***
      logical value

c     **** deallocate from heap ****
      value =           BA_free_heap(Tc(2))
      value = value.and.BA_free_heap(Ts(2))
      value = value.and.BA_free_heap(rgrid(2))
      if (.not.value) 
     >   call errquit('dipole_efield_end:freeing memory',0,MA_ERR)

      return
      end 


*     ******************************
*     *                            *
*     *   dipole_Efield_set_ion    *
*     *                            *
*     ******************************
      subroutine dipole_Efield_set_ion()
      implicit none

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     *** local variables ***
      integer i,ia

*     *** external functions ***
      integer  ion_katm,ion_nion
      external ion_katm,ion_nion
      real*8   ion_rion,psp_zv
      external ion_rion,psp_zv

      nx=0.0d0
      ny=0.0d0
      nz=0.0d0
      pcharge = 0.0d0
      do i=1,ion_nion()
        ia=ion_katm(i)
        nx=nx+psp_zv(ia)*ion_rion(1,i)
        ny=ny+psp_zv(ia)*ion_rion(2,i)
        nz=nz+psp_zv(ia)*ion_rion(3,i)
        pcharge = pcharge + psp_zv(ia)
      end do   

      return
      end 


*     ******************************
*     *                            *
*     *   dipole_Efield_set_bv     *
*     *                            *
*     ******************************
      subroutine dipole_Efield_set_bv()
      implicit none

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     *** local variables ***
      integer tmp_len
      parameter (tmp_len=140)

      integer i,j,info
      real*8 wrk(6,6),bmat(3,3),ixmat(3,6)

      real*8 xx,yy,zz,tmp1(tmp_len)
      real*8 tx,ty,tz,alpha,scal

*     *** external functions ***
      real*8   lattice_unitg,lattice_omega
      external lattice_unitg,lattice_omega

c     *** Silvestrelli G1 ***
      ixmat(1,1)=1.0d0
      ixmat(2,1)=0.0d0
      ixmat(3,1)=0.0d0

c     *** Silvestrelli G4 ***
      ixmat(1,2)=1.0d0
      ixmat(2,2)=1.0d0
      ixmat(3,2)=0.0d0

c     *** Silvestrelli G5 ***
      ixmat(1,3)=1.0d0
      ixmat(2,3)=0.0d0
      ixmat(3,3)=1.0d0

c     *** Silvestrelli G2 ***
      ixmat(1,4)=0.0d0
      ixmat(2,4)=1.0d0
      ixmat(3,4)=0.0d0

c     *** Silvestrelli G6 ***
      ixmat(1,5)=0.0d0
      ixmat(2,5)=1.0d0
      ixmat(3,5)=1.0d0

c     *** Silvestrelli G3 ***
      ixmat(1,6)=0.0d0
      ixmat(2,6)=0.0d0
      ixmat(3,6)=1.0d0

      do i=1,3
         bmat(i,1)=lattice_unitg(1,i)
         bmat(i,2)=lattice_unitg(2,i)
         bmat(i,3)=lattice_unitg(3,i)
      end do

      do i=1,6
         xx=0.0d0
         yy=0.0d0
         zz=0.0d0
         do j=1,3
           xx=xx+bmat(j,1)*ixmat(j,i)
           yy=yy+bmat(j,2)*ixmat(j,i)
           zz=zz+bmat(j,3)*ixmat(j,i)
         end do
         bv(1,i)=xx
         bv(2,i)=yy
         bv(3,i)=zz
      end do

      do i=1,6
         wrk(1,i)=bv(1,i)*bv(1,i)
         wrk(2,i)=bv(1,i)*bv(2,i)
         wrk(3,i)=bv(1,i)*bv(3,i)
         wrk(4,i)=bv(2,i)*bv(2,i)
         wrk(5,i)=bv(2,i)*bv(3,i)
         wrk(6,i)=bv(3,i)*bv(3,i)
         wts(i)=0.0d0
      end do

*     *** scal=(2*pi/L)**2 ***
      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0d0*datan(1.0d0)/scal
      scal = scal*scal
      wts(1)=scal
      wts(4)=scal
      wts(6)=scal
      call DGELS('N',6,6,1,wrk,6,wts,6,tmp1,tmp_len,info)
      if (info.ne.0) then
        write(*,*)"Illegal argument in call to dgels"
        call flush(6)
      end if
      rank=0
      do i=1,6
         if (dabs(wts(i)).gt.1.e-6) then
           rank=rank+1
           wrk(1,rank)=bv(1,i)
           wrk(2,rank)=bv(2,i)
           wrk(3,rank)=bv(3,i)
           wrk(4,rank)=wts(i)
         end if
      end do
      do i=1,rank
         bv(1,i)=wrk(1,i)
         bv(2,i)=wrk(2,i)
         bv(3,i)=wrk(3,i)
         wts(i)=wrk(4,i)
      end do

      return
      end



*     ******************************
*     *                            *
*     *      dipole_Efield_print   *
*     *                            *
*     ******************************
      subroutine dipole_Efield_print(iunit)
      implicit none
      integer iunit

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     *** local variables ***
      real*8 autoDebye
      parameter (autoDebye=2.5416d0)

      real*8 xx
      xx = dsqrt(dipole(1)**2 + dipole(2)**2 + dipole(3)**2)

      write(iunit,1771)
      write(iunit,1772) 'spin up   ',
     >                  dx(1)/dble(ne0(1)),
     >                  dy(1)/dble(ne0(1)),
     >                  dz(1)/dble(ne0(1))
      if (ne0(ispin0).ne.0)
     >   write(iunit,1772) 'spin down ',
     >                     dx(ispin0)/dble(ne0(ispin0)),
     >                     dy(ispin0)/dble(ne0(ispin0)),
     >                     dz(ispin0)/dble(ne0(ispin0))
      write(iunit,1772) 'electronic',
     >                   (dx(1)+dx(ispin0))/dble(ne0(1)+ne0(ispin0)),
     >                   (dy(1)+dy(ispin0))/dble(ne0(1)+ne0(ispin0)),
     >                   (dz(1)+dz(ispin0))/dble(ne0(1)+ne0(ispin0))
      write(iunit,1772) 'ionic     ',
     >                   nx/pcharge,
     >                   ny/pcharge,
     >                   nz/pcharge
      write(iunit,1778)
      write(iunit,1774) dipole
      write(iunit,1775) xx,xx*autoDebye

      return
*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1771 FORMAT(//'== Center of Charge =='/)
 1772 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')
 1773 FORMAT(//'== Wannier Crystal Dipole =='/)
 1774 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1775 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1776 FORMAT(/"ELECTRONIC DIPOLES")
 1777 FORMAT("DX =",F11.5," DY= ",F11.5," DZ= ",F11.5)
 1778 FORMAT(//'== Crystal Dipole (Resta) =='/)
 1780 FORMAT("NUCLEAR DIPOLES")
 1785 FORMAT("TOTAL DIPOLES")
      end 

*     ******************************
*     *                            *
*     *   dipole_Efield_gen_TcTs   *
*     *                            *
*     ******************************
      subroutine dipole_Efield_set_TcTs(i,ms,neq,n2ft3d,psi_r,psi_r2,W)
      implicit none
      integer i,ms,neq(2),n2ft3d
      real*8  psi_r(*),psi_r2(*)
      complex*16 W(*)

#include "bafdecls.fh"
#include "errquit.fh"

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     **** local variables ****
      logical value
      integer j,k,n1,n2,n3
      real*8  br,scal1
      integer Wc(2),Ws(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      
      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      scal1 = 1.0d0/dble(n1*n2*n3)

c     **** allocate memory from stack ****
      value =           Dneall_m_push_get(ms,Wc)
      value = value.and.Dneall_m_push_get(ms,Ws)
      if (.not. value)
     >  call errquit('dipole_Efield_set_TcTs:out of stack memory',0,0)


      !*** generate Tc and Ts ***
!$OMP DO
      do k=1,n2ft3d
         br = bv(1,i)*dbl_mb(rgrid(1)+(k-1)*3)
     >      + bv(2,i)*dbl_mb(rgrid(1)+(k-1)*3 + 1)
     >      + bv(3,i)*dbl_mb(rgrid(1)+(k-1)*3 + 2)
         dbl_mb(Tc(1)+k-1) =  cos(br)
         dbl_mb(Ts(1)+k-1) = -sin(br)
      end do
!$OMP END DO

*     **** generate W = <psi_r(i)|exp(-i b*r)|psi_r(j)> ****
      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Tc(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_ggm_sym_Multiply(ms,psi_r,psi_r2,n2ft3d,
     >                             dbl_mb(Wc(1)))

      call Dneall_m_scal(ms,scal1,dbl_mb(Wc(1)))

      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Ts(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_ggm_sym_Multiply(ms,psi_r,psi_r2,n2ft3d,
     >                             dbl_mb(Ws(1)))
      
      call Dneall_m_scal(ms,scal1,dbl_mb(Ws(1)))
      call Dneall_mmtow_Cmplx(ms,dbl_mb(Wc(1)),dbl_mb(Ws(1)),W)


*     **** pop memory ***
      value =           Dneall_m_pop_stack(Ws)
      value = value.and.Dneall_m_pop_stack(Wc)
      if (.not. value)
     >  call errquit('dipole_Efield_set_TcTs:popping stack memory',1,0)

      return
      end


c      subroutine dipole_Efield_gen_TcTs(Z,n2ft3d,ZTcTs)
c      implicit none
c      complex*16 Z
c      integer    n2ft3d
c      real*8     Tc(*),Ts(*)
c      real*8     ZTcTs(*)
c      integer i
c      do i=1,n2ft3d
c         ZTcTs(i) = dble(Z)*Tc(i) + dimag(Z)*Ts(i)
c      end do
c      return
c      end 

*     ******************************
*     *                            *
*     *   dipole_Efield_add_dipole *
*     *                            *
*     ******************************
      subroutine  dipole_Efield_add_dipole(ms,r,Z)
      implicit none
      integer ms,r
      complex*16 Z(*)
       
      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     **** local variables ***
      integer    i
      real*8     scal
      complex*16 arg

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0d0*datan(1.0d0)/scal
      scal = scal*scal

!$OMP MASTER
      !*** really just want complex eigenvalues of X here ***
      do i=1,ne0(ms)
         arg = -wts(r)*datan2(dimag(Z(i)),dble(Z(i)))

         dx(ms) = dx(ms) + bv(1,r)*arg/scal
         dy(ms) = dy(ms) + bv(2,r)*arg/scal
         dz(ms) = dz(ms) + bv(3,r)*arg/scal
      end do !* i *
!$OMP END MASTER

      return
      end

*     ******************************
*     *                            *
*     *    dipole_Efield_gen_Cij   *
*     *                            *
*     ******************************
      subroutine  dipole_Efield_gen_Cij(N,Z,U,V,ZENOM,Creal,Cimg)
      implicit none
      integer N
      complex*16 Z(N),U(N,N),V(N,N),ZENOM(N)
      real*8 Creal(N,N),Cimg(N,N)

*     **** local variables ****
      integer i,j,k
      complex*16 zzz

!$OMP DO
      do k=1,N
         ZENOM(k) = dcmplx(0.0d0,0.0d0)
         do i=1,N
            ZENOM(k) = ZENOM(k) + dconjg(U(i,k))*V(i,k)
         end do
         !write(*,*) "k,ZENOM=",k,ZENOM(k)
      end do
!$OMP END DO

!$OMP DO
      do i=1,N
      do j=1,N
         zzz = dcmplx(0.0d0,0.0d0)
         do k=1,N
c            zzz = zzz
c     >          + (dcmplx(-dimag(Z(k)),dble(Z(k)))
c     >            /(dble(Z(k))**2 + dimag(Z(k))**2))
c     >           * (dconjg(U(i,k))*V(j,k))/ZENOM(k)

            zzz = zzz
     >          + (dcmplx(dimag(Z(k)),dble(Z(k)))
     >            /(dble(Z(k))**2 + dimag(Z(k))**2))
     >           * (dconjg(U(i,k))*V(j,k))/ZENOM(k)

         end do
         Creal(i,j) = dble(zzz)
         Cimg(i,j)  = dimag(zzz)
      end do
      end do
!$OMP END DO 

      return
      end 


*     ******************************
*     *                            *
*     *  dipole_Efield_gen_hpsi_r  *
*     *                            *
*     ******************************
      subroutine dipole_Efield_gen_hpsi_r(ms,r,neq,n2ft3d,
     >                                    psi_r,psi_r2,
     >                                    Creal,Cimg,hpsi_r)
      implicit none
      integer ms,r,neq(2),n2ft3d
      real*8  psi_r(*),psi_r2(*)
      real*8  Creal(*),Cimg(*)
      real*8  hpsi_r(*) 

#include "bafdecls.fh"
#include "errquit.fh"

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     **** local variables ****
      integer j,k
      real*8 alpha,scal,br

*     *** external variables ****
      real*8   lattice_omega
      external lattice_omega


      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0d0*datan(1.0d0)/scal
      scal = scal*scal

      alpha = wts(r)*( efield(1)*bv(1,r)
     >               + efield(2)*bv(2,r)
     >               + efield(3)*bv(3,r))
     >              / scal

      !write(*,*) "alpha=",alpha,Creal(1),Cimg(1)


      !*** generate Tc and Ts ***
!$OMP DO
      do k=1,n2ft3d
         br = bv(1,r)*dbl_mb(rgrid(1)+(k-1)*3)
     >      + bv(2,r)*dbl_mb(rgrid(1)+(k-1)*3 + 1)
     >      + bv(3,r)*dbl_mb(rgrid(1)+(k-1)*3 + 2)
         dbl_mb(Tc(1)+k-1) =  cos(br)
         dbl_mb(Ts(1)+k-1) = -sin(br)
      end do
!$OMP END DO


*     **** generate W = <psi_r(i)|exp(-i b*r)|psi_r(j)> ****
      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Tc(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_gmg_Multiply(ms,psi_r2,n2ft3d,
     >                         Creal,-alpha,
     >                         hpsi_r,1.0d0)

      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Ts(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_gmg_Multiply(ms,psi_r2,n2ft3d,
     >                         Cimg,alpha,
     >                         hpsi_r,1.0d0)


      return
      end


*     ******************************
*     *                            *
*     *      dipole_Efield_Vnl     *
*     *                            *
*     ******************************
      subroutine dipole_Efield_Vnl(ispin,neq,n2ft3d,psi_r,hpsi_r,edpol)
      implicit none
      integer ispin,neq(2),n2ft3d
      real*8     psi_r(n2ft3d,*)
      real*8     hpsi_r(n2ft3d,*)
      real*8     edpol

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

*     **** local variables ****
      logical value
      integer i,j,k,r,ms
      integer n1,n2,n3

      real*8 b(3)
      real*8 xx,yy,zz,dv
      real*8 tx,ty,tz,scal,scal1
      complex*16 arg,wx

      integer X(2,6),Xeig(2,6)
      integer U(2,6),V(2,6),UV(2)
      integer Creal(2),Cimg(2)
      integer psi_r2(2)

*     **** external functions ****
      logical  Dneall_w_push_get,Dneall_w_pop_stack,control_print
      external Dneall_w_push_get,Dneall_w_pop_stack,control_print
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      real*8   lattice_omega,Dneall_m_trace
      external lattice_omega,Dneall_m_trace


      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      scal1 = 1.0d0/dble(n1*n2*n3)

*     *** scal=(2*pi/L)**2 ***
      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0d0*datan(1.0d0)/scal
      scal = scal*scal

*     ***** allocate X,Y,Z  ****
      value = .true.
      do j=1,6
         value = value.and.Dneall_w_push_get(1,X(1,j))
         value = value.and.
     >           BA_push_get(mt_dcpl,ne0(1),'Xeig',Xeig(2,j),Xeig(1,j))
         value = value.and.Dneall_w_push_get(1,U(1,j))
         value = value.and.Dneall_w_push_get(1,V(1,j))
      end do
      value = value.and.BA_push_get(mt_dcpl,ne0(1),'UV',UV(2),UV(1))
      value = value.and.Dneall_m_push_get(0,Creal)
      value = value.and.Dneall_m_push_get(1,Cimg)
      value = value.and.
     >        BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,
     >                    'Hpsi_r',psi_r2(2),psi_r2(1))
      if (.not.value) 
     >   call errquit('dipole_efield_Vnl:out of stack',0,MA_ERR)


!$OMP MASTER
      dx(1)=0.0d0
      dy(1)=0.0d0
      dz(1)=0.0d0
      dx(2)=0.0d0
      dy(2)=0.0d0
      dz(2)=0.0d0
!$OMP END MASTER

      
      do ms=1,ispin
         do r=1,rank
            call dipole_Efield_set_TcTs(r,ms,neq,n2ft3d,
     >                 psi_r,dbl_mb(psi_r2(1)),dcpl_mb(X(1,r)))

            call Dneall_w_eigenvaluesvectors(ms,dcpl_mb(X(1,r)),
     >                              dcpl_mb(Xeig(1,r)),
     >                              dcpl_mb(U(1,r)),
     >                              dcpl_mb(V(1,r)))

            call dipole_Efield_add_dipole(ms,r,dcpl_mb(Xeig(1,r)))
            call dipole_Efield_gen_Cij(ne0(ms),
     >                                 dcpl_mb(Xeig(1,r)),
     >                                 dcpl_mb(U(1,r)), 
     >                                 dcpl_mb(V(1,r)),
     >                                 dcpl_mb(UV(1)),
     >                                 dbl_mb(Creal(1)),dbl_mb(Cimg(1)))
             call dipole_Efield_gen_hpsi_r(ms,r,neq,n2ft3d,
     >                                 psi_r,
     >                                 dbl_mb(psi_r2(1)),
     >                                 dbl_mb(Creal(1)),dbl_mb(Cimg(1)),
     >                                 hpsi_r)
         end do !* r *
      end do !* ms *

      call Dneall_ggm_sym_Multiply(0,psi_r,hpsi_r,n2ft3d,
     >                             dbl_mb(Creal(1)))
!$OMP MASTER
      Pdipole = scal1*Dneall_m_trace(0,dbl_mb(Creal(1)))


cccccccccccccccccccccccccccccccccccccccccccccc
c  Molecular dipoles from Resta's theory!
ccccccccccccccccccccccccccccccccccccccccccccc
      tx=nx-dx(1)-dx(ispin)
      ty=ny-dy(1)-dy(ispin)
      tz=nz-dz(1)-dz(ispin)
      xx = dsqrt(tx*tx + ty*ty + tz*tz)
      dipole(1) = tx
      dipole(2) = ty
      dipole(3) = tz
!$OMP END MASTER
!$OMP BARRIER

      !write(*,*) "DipOlE=",dipole
      edpol = dipole(1) * efield(1) 
     >      + dipole(2) * efield(2) 
     >      + dipole(3) * efield(3) 

      Edipole = edpol


c     **** pop stack ****
      value = BA_pop_stack(psi_r2(2))
      value = value.and.Dneall_m_pop_stack(Cimg)
      value = value.and.Dneall_m_pop_stack(Creal)
      value = value.and.BA_pop_stack(UV(2))
      do j=6,1,-1
         value = value.and.Dneall_w_pop_stack(V(1,j))
         value = value.and.Dneall_w_pop_stack(U(1,j))
         value = value.and.BA_pop_stack(Xeig(2,j))
         value = value.and.Dneall_w_pop_stack(X(1,j))
      end do
      if (.not.value) 
     >   call errquit('dipole_efield_vnl:pop stack',2,MA_ERR)


      !write(*,*) "ispin,r,Pdipole,Edipole = ", ispin,r,Pdipole,Edipole
      !write(*,*) 

      return 
      end


*     ******************************
*     *                            *
*     *       dipole_Efield_e      *
*     *                            *
*     ******************************
      real*8 function dipole_Efield_e()
      implicit none

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

      dipole_efield_e = Edipole
      return
      end 


*     ******************************
*     *                            *
*     *       dipole_Efield_p      *
*     *                            *
*     ******************************
      real*8 function dipole_Efield_p()
      implicit none

      !**** dipole_efield common block ****
      real*8 efield(3),efield_center(3)
      real*8 dx(2),dy(2),dz(2)
      real*8 nx,ny,nz,pcharge
      real*8 dipole(3),Edipole,Pdipole
      real*8 bv(3,6),wts(6)
      integer rank
      integer Tc(2),Ts(2),rgrid(2)
      integer ispin0,ne0(2)
      common /dipole_efield/ efield,efield_center,
     >                       dx,dy,dz,
     >                       nx,ny,nz,pcharge,dipole,
     >                       Edipole,Pdipole,
     >                       bv,wts,rank,
     >                       Tc,Ts,rgrid,
     >                       ispin0,ne0

      dipole_efield_p = Pdipole
      return
      end




