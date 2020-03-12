*     *****************************************
*     *                                       *
*     *             vdw_DF_init               *
*     *                                       *
*     *****************************************
*
*
      subroutine vdw_DF_init()
      implicit none

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"

#include "vdw-DF.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical mprint,hprint,debug,does_it_exist,oprint,value
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j,k,l,nx
      integer Gindx(3)
      real*8  dk,gg,gx,gy,gz

      character*255 datafile
      logical found_datafile
      integer ifound

*     **** external functions ***
      logical  control_has_vdw,control_is_vdw2
      external control_has_vdw,control_is_vdw2
      integer  nwpw_splint_nx,Pack_G_indx
      external nwpw_splint_nx,Pack_G_indx


*     **** return if not has_vdw ****
      has_vdw = control_has_vdw()
      is_vdw2 = control_is_vdw2()
      if (.not.has_vdw) return

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call util_file_name_noprefix('vdw_kernels.dat',
     >                                .false.,
     >                                .false.,
     >                                datafile)
      if (taskid.eq.MASTER) then
         ifound = 0
         inquire(file=datafile,exist=found_datafile)
         if (found_datafile) ifound = 1
      end if
      call Parallel_Brdcst_ivalue(MASTER,ifound)
      found_datafile = .false.
      if (ifound.eq.1) found_datafile=.true.

      !**** generate vdw_kernels.dat file ****
      if (.not.found_datafile) then
         if (oprint) then
            l = index(datafile,' ') - 1
            write(luout,*) "Generating VDW kernel filename:",
     >                      datafile(1:l)
         end if
         call vdw_DF_kernel_gen_data(datafile)
      end if


      if (oprint) then
         l = index(datafile,' ') - 1
         write(luout,*) "Reading VDW kernel filename:",datafile(1:l)

      end if
*     **** create vdw data file ****
      l = index(datafile,' ') -1

      if (taskid.eq.MASTER) then
         call openfile(5,datafile,l,'r',l)
         call iread(5,Nqs,1)
         call iread(5,nk,1)
         call dread(5,kmax,1)
      end if
      call Parallel_Brdcst_ivalue(MASTER,Nqs)
      call Parallel_Brdcst_value(MASTER,kmax)
      call Parallel_Brdcst_ivalue(MASTER,nk)
      nk1 = nk + 1

      value= BA_alloc_get(mt_dbl,Nqs,'vdw_qmesh',qmesh(2),qmesh(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,Nqs*Nqs,'vdw_y',  ya(2), ya(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,Nqs*Nqs,'vdw_y2',y2a(2),y2a(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,nk1,'vdw_g',gphi(2),gphi(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,nk1*Nqs*(Nqs+1),
     >                    'vdw_phi',phi(2),phi(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,Nqs*n2ft3d,
     >                    'vdw_theta',theta(2),theta(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,Nqs*n2ft3d,
     >                    'vdw_ufunc',ufunc(2),ufunc(1))
      value= value.and.
     >       BA_alloc_get(mt_dbl,10*n2ft3d,
     >                    'vdw_xcp',xcp(2),xcp(1))
      xce(1)  = xcp(1) + 2*n2ft3d
      xxp(1)  = xcp(1) + 4*n2ft3d
      xxe(1)  = xcp(1) + 6*n2ft3d
      rho(1)  = xcp(1) + 8*n2ft3d
      value= value.and.
     >       BA_alloc_get(mt_dbl,npack0,
     >                    'vdw_Gpack',Gpack(2),Gpack(1))
      value= value.and.
     >       BA_alloc_get(mt_int,npack0,
     >                    'vdw_nxpack',nxpack(2),nxpack(1))
      if (.not.value) call errquit('vdw_DF_init:out of heap',2,MA_ERR) 

      if (taskid.eq.MASTER) then
         call dread(5,dbl_mb(qmesh(1)),Nqs)
         call dread(5,dbl_mb(phi(1)),nk1*Nqs*(Nqs+1))
         call closefile(5)

         !*** set qmin = 0.0d0 ***
         dbl_mb(qmesh(1)) = 0.0d0

         !***  set phi(:,1,1) = 0 ***
         call dcopy(2*nk1,0.0d0,0,dbl_mb(phi(1)),1)

      end if
      call Parallel_Brdcst_values(MASTER,Nqs,dbl_mb(qmesh(1)))
      call Parallel_Brdcst_values(MASTER,nk1*Nqs*(Nqs+1),
     >                            dbl_mb(phi(1)))

      call vdw_DF_init_poly(Nqs,dbl_mb(qmesh(1)),
     >                      dbl_mb(ya(1)),
     >                      dbl_mb(y2a(1)),
     >                      dbl_mb(gphi(1)))

      dk = kmax/dble(nk)
      do k=0,nk
         dbl_mb(gphi(1)+k) = k*dk
      end do

      Gindx(1)=Pack_G_indx(0,1)
      Gindx(2)=Pack_G_indx(0,2)
      Gindx(3)=Pack_G_indx(0,3)
      do k=1,npack0
         gx = dbl_mb(Gindx(1)+k-1)
         gy = dbl_mb(Gindx(2)+k-1)
         gz = dbl_mb(Gindx(3)+k-1)
         gg = dsqrt(gx*gx + gy*gy + gz*gz)
         dbl_mb(Gpack(1)+k-1) = gg

         nx = gg/dk
         int_mb(nxpack(1)+k-1) = nwpw_splint_nx(dbl_mb(gphi(1)),nx,gg)
      end do

      if (is_vdw2) then
         Zab = -1.887d0
      else
         Zab = -0.8491d0
      end if
      qmax = dbl_mb(qmesh(1)+Nqs-1)
      qmin = 0.0d0
    
      return

 999  continue
      call errquit('vdw_DF_init:error reading qmesh,Nqs=',Nqs,DISK_ERR)
      return

      end


*     **********************************************
*     *                                            *
*     *                vdw_DF_end                  *
*     *                                            *
*     **********************************************

      subroutine vdw_DF_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "vdw-DF.fh"

*     **** local variables ****
      logical value

      if (has_vdw) then
         value = BA_free_heap(nxpack(2))
         value = value.and.BA_free_heap(Gpack(2))
         value = value.and.BA_free_heap(xcp(2))
         value = value.and.BA_free_heap(ufunc(2))
         value = value.and.BA_free_heap(theta(2))
         value = value.and.BA_free_heap(phi(2))
         value = value.and.BA_free_heap(gphi(2))
         value = value.and.BA_free_heap(qmesh(2))
         value = value.and.BA_free_heap(ya(2))
         value = value.and.BA_free_heap(y2a(2))
         if (.not.value) 
     >      call errquit('vdw_DF_end:free heap failed',0,MA_ERR)
      end if

      return
      end

*     **********************************************
*     *                                            *
*     *                vdw_DF_init_poly            *
*     *                                            *
*     **********************************************
      subroutine vdw_DF_init_poly(n,x,ya,y2a,utmp)
      implicit none
      integer n
      real*8 x(*),ya(n,*),y2a(n,*),utmp(*)

*     **** local variables ****
      integer i,j
      real*8 yp1,ypn

      yp1 = 0.0d0
      ypn = 0.0d0

      do i=1,n
         do j=1,n
            if (i.eq.j) then
               ya(i,j) = 1.0d0
            else
               ya(i,j) = 0.0d0
            end if
         end do
      end do
      do j=1,n
         call nwpw_spline(x,ya(1,j),n,yp1,ypn,y2a(1,j),utmp)
      end do

      return
      end

*     **********************************************
*     *                                            *
*     *                vdw_DF_poly                 *
*     *                                            *
*     **********************************************
      subroutine vdw_DF_poly(i,x,p,dp)
      implicit none
      integer i
      real*8 x
      real*8 p,dp

#include "bafdecls.fh"
#include "errquit.fh"
#include "vdw-DF.fh"

      !**** local variables ****
      integer nx

      !**** external functions ****
      real*8   nwpw_splint,nwpw_dsplint
      external nwpw_splint,nwpw_dsplint


      if ((x.ge.qmin).and.(x.le.qmax)) then
         nx = dsqrt(x/dbl_mb(qmesh(1)+Nqs-1)) * Nqs
         if (nx.lt.1)       nx = 1
         if (nx.gt.(Nqs-1)) nx = Nqs-1

         p = nwpw_splint(dbl_mb(qmesh(1)),
     >                   dbl_mb(ya(1)+Nqs*(i-1)),
     >                   dbl_mb(y2a(1)+Nqs*(i-1)),
     >                   Nqs,nx,x)
         dp = nwpw_dsplint(dbl_mb(qmesh(1)),
     >                     dbl_mb(ya(1)+Nqs*(i-1)),
     >                     dbl_mb(y2a(1)+Nqs*(i-1)),
     >                     Nqs,nx,x)
      else
         p = 0.0d0
         dp = 0.0d0
      end if
      return
      end 


*     **********************************************
*     *                                            *
*     *                vdw_DF_exist                *
*     *                                            *
*     **********************************************
      logical function vdw_DF_exist()
      implicit none

#include "vdw-DF.fh"

      vdw_DF_exist = has_vdw
      return
      end

*     **********************************************
*     *                                            *
*     *                vdw_DF                      *
*     *                                            *
*     **********************************************
*      input:  n2ft3d                  grid
*              dn                     density
*              agr                     absolute gradient of density,|grad n|
*      output: exc                     exchange correlation energy density
*              fn                      d(n*exc)/dnup,
*              fdn                     d(n*exc)/d(|grad n|)

      subroutine vdw_DF(n2ft3d0,ispin,dn,agr,exc,fn,fdn)
      implicit none
      integer n2ft3d0,ispin
      real*8 dn(n2ft3d0,ispin)
      real*8 agr(n2ft3d0,ispin)
      real*8 exc(n2ft3d0)
      real*8 fn(n2ft3d0,ispin)
      real*8 fdn(n2ft3d0,ispin)

#include "bafdecls.fh"
#include "errquit.fh"
#include "vdw-DF.fh"


      !*** generate LDA results ***
      call vxc(n2ft3d,ispin,dn,dbl_mb(xcp(1)),dbl_mb(xce(1)),
     >         dbl_mb(rho(1)))
      if (ispin.eq.1) 
     >   call v_dirac(n2ft3d,ispin,dn,dbl_mb(xxp(1)),dbl_mb(xxe(1)),
     >                dbl_mb(rho(1)))

      !**** generate rho ***
      call vdw_DF_Generate_rho(ispin,n2ft3d,dn,dbl_mb(rho(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(rho(1)))

      !*** Generate theta(G), ptheta(r), dthetadrho(r,ms), dthetaddrho(r) ****
      call vdw_DF_Generate_thetag(Nqs,nfft3d,ispin,n2ft3d,
     >                            Zab,qmin,qmax,dbl_mb(rho(1)),agr,
     >                            dbl_mb(xcp(1)),dbl_mb(xce(1)),
     >                            dbl_mb(xxp(1)),dbl_mb(xxe(1)),
     >                            dbl_mb(theta(1)))

      !*** compute ufunc(G,i) = Sum(j) theta(G,j)*phi(G,i,j) ***
      call vdw_DF_Generate_ufunc(nk1,Nqs,
     >                           dbl_mb(gphi(1)),dbl_mb(phi(1)),
     >                           npack0,nfft3d,
     >                           dbl_mb(Gpack(1)),int_mb(nxpack(1)),
     >                           dbl_mb(theta(1)),dbl_mb(ufunc(1)))

         
      !*** compute contributions to xce(r), fn(r,ms), fdn ****
      call vdw_DF_Generate_potentials(Nqs,nfft3d,ispin,n2ft3d,
     >             dbl_mb(ufunc(1)),
     >             dbl_mb(xce(1)),dbl_mb(xcp(1)),dbl_mb(xxe(1)),
     >             dbl_mb(xce(1)+n2ft3d),dbl_mb(xxp(1)),dbl_mb(rho(1)),
     >             exc,fn,fdn)


      return
      end 

*     ************************************************
*     *                                              *
*     *         vdw_DF_Generate_rho                  *
*     *                                              *
*     ************************************************
      subroutine vdw_DF_Generate_rho(ispin,n2ft3d,dn,rho)
      implicit none
      integer ispin,n2ft3d
      real*8 dn(n2ft3d,ispin)
      real*8 rho(n2ft3d)

*     **** local variables ****
      integer i,tid,nthr
      real*8 dncut
      parameter (dncut=1.0d-30)

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      do i=tid+1,n2ft3d,nthr
         rho(i) = dn(i,1)+dn(1,ispin)+dncut
      end do
!$OMP BARRIER

      return
      end

*     ************************************************
*     *                                              *
*     *         vdw_DF_Generate_thetag               *
*     *                                              *
*     ************************************************
*
      subroutine vdw_DF_Generate_thetag(Nqs,nfft3d,ispin,n2ft3d,
     >                                 Zab,qmin,qmax,
     >                                 rho,agr,vxc,exc,vxx,exx,
     >                                 theta)
      implicit none
      integer Nqs,nfft3d,ispin,n2ft3d
      real*8  Zab,qmin,qmax
      real*8  rho(n2ft3d),agr(n2ft3d)
      real*8  vxc(n2ft3d,ispin), exc(n2ft3d)
      real*8  vxx(n2ft3d,ispin), exx(n2ft3d)
      real*8  theta(n2ft3d,Nqs)

*     **** local variables ****
      integer i,j,k,jstart,nj,taskid_j,np_j,tid,nthr
      integer nx,ny,nz,r,ms
      real*8 A,Cf,pi,pj,dpj,tsum,dtsum,q0,q0sat,xi,dxi,Cxi,scal1
      real*8 dq0drho(2),dq0ddrho,dq0satdq0
      real*8 vxgga,exgga
      real*8 onethird,frthrd,elthrd,dncut,check
      parameter (onethird=1.0d0/3.0d0)
      parameter (frthrd=4.0d0/3.0d0)
      parameter (elthrd=11.0d0/3.0d0)
      parameter (dncut=1.0d-12)

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      nx = Nqs/np_j
      r = mod(Nqs,np_j)
      if (taskid_j.lt.r) then
         nj = nx + 1
         jstart = nx*taskid_j + taskid_j + 1
      else
         nj = nx
         jstart = nx*taskid_j + r + 1
      end if

      call D3dB_r_nZero(1,Nqs,theta)

      if (nj.gt.0) then
         pi = 4.0d0*datan(1.0d0)
         Cf = (3.0d0*pi*pi)**onethird
         A = -4.0d0*pi/3.0d0
         Cxi = Zab/(36.0d0*Cf)

         !*** compute theta(r) ***
         do i=tid+1,n2ft3d,nthr
            if (rho(i).ge.dncut) then
               xi  = Cxi*(agr(i)/rho(i)**frthrd)**2
               dxi = 2.0d0*Cxi*agr(i)*(1.0d0/rho(i)**frthrd)**2
               exgga = exx(i) - xi*exx(i)
               vxgga = vxx(i,1) - xi*(vxx(i,1)-elthrd*exx(i))
               if ((vxgga-exgga).lt.0.0d0) then
                  xi = 0.0d0
                  dxi = 0.0d0
               end if

               q0 = A*(exc(i) - exx(i)*xi)

               if (q0.lt.qmin) then
                  q0 = qmin
                  do ms=1,ispin
                     dq0drho(ms)  = 0.0d0
                  end do
                  dq0ddrho = 0.0d0
               else
                  do ms=1,ispin
                     dq0drho(ms) = A*((vxc(i,ms)-exc(i))
     >                           - xi*(vxx(i,ms)-elthrd*exx(i)))/rho(i)
                  end do
                  dq0ddrho = -A*exx(i)*dxi
               end if

               tsum = 0.0d0
               dtsum = 0.0d0
               do k=1,12
                  tsum = tsum + (q0/qmax)**k / dble(k)
                  dtsum = dtsum + (q0/qmax)**(k-1) 
               end do
               q0sat     = qmax*(1.0d0-dexp(-tsum))
               dq0satdq0 = dexp(-tsum)*dtsum

               exc(i)  =  q0sat
               do ms=1,ispin
                  vxc(i,ms) = rho(i)*dq0satdq0*dq0drho(ms)
               end do
               exx(i)  =  rho(i)*dq0satdq0*dq0ddrho
            else
               xi = 0.0d0
               q0sat    = qmax
               exc(i) = qmax
               do ms=1,ispin
                  vxc(i,ms) = 0.0d0
               end do
               exx(i) = 0.0d0
            end if

            do j=jstart,jstart-1+nj
               call vdw_DF_poly(j,q0sat,pj,dpj)
               theta(i,j)  = rho(i)*pj
            end do
         end do
!$OMP BARRIER

         !*** compute theta(g) ***
         do j=jstart,jstart-1+nj
            call D3dB_r_Zero_Ends(1,theta(1,j))
         end do
         call Grsm_hg_fftf(nfft3d,nj,theta(1,jstart))
         call Grsm_gg_dScale1(nfft3d,nj,scal1,theta(1,jstart))
      end if

      call D1dB_Vector_SumAll(Nqs*n2ft3d,theta)

      return
      end


*     ************************************************
*     *                                              *
*     *         vdw_DF_Generate_ufunc                *
*     *                                              *
*     ************************************************
*
*   This routine computes the 
*
*   ufunc(k,i) = Sum(j=1,Nqs) theta(k,j) * phi(k,i,j)
*
*   where phi is interpolated using a cubic spline.
*
      subroutine vdw_DF_Generate_ufunc(nk1,Nqs,gphi,phi,
     >                                 npack0,nfft3d,
     >                                 Gpack,nxpack,theta,ufunc)
      implicit none
      integer nk1,Nqs
      real*8  gphi(nk1),phi(nk1,2,Nqs*(Nqs+1)/2)
      integer    npack0,nfft3d
      real*8     Gpack(npack0)
      integer    nxpack(npack0)
      complex*16 theta(nfft3d,*)
      complex*16 ufunc(nfft3d,*)

*     **** local variables ****
      integer tid,nthr,taskid_j,np_j,pcount
      integer i,j,k,indx,klo,khi
      real*8 a,b,f,g,h

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      
c      write(*,*) "tid,nthr=",tid,nthr
c      write(*,*) "taskid_j,np_j=",taskid_j,np_j
      !*** compute ufunc(g) ***
      call D3dB_c_nZero(1,Nqs,ufunc)
      pcount = 0 
      indx = 0
      do j=1,Nqs    !*** assuming phi(:,1,1) = 0 ***
         do i=1,j
            indx = indx + 1
c            write(*,*) "i,j=",i,j, indx, Nqs,npack0,nfft3d
            if (pcount.eq.taskid_j) then

               do k=tid+1,npack0,nthr
                  g  = Gpack(k)
                  klo  = nxpack(k)
                  khi  = klo + 1
                  !write(*,*) "k,g,klo,khi=",k,g,klo,khi,nk1

                  h = gphi(khi)-gphi(klo)
                  a = (gphi(khi)-g)/h
                  b = (g-gphi(klo))/h
                  f = a*phi(klo,1,indx)
     >              + b*phi(khi,1,indx)
     >              + ((a**3-a)*phi(klo,2,indx)
     >                +(b**3-b)*phi(khi,2,indx))*h**2/6.0d0

                  ufunc(k,i) = ufunc(k,i) + theta(k,j) * f
                  if (i.ne.j)  ufunc(k,j) = ufunc(k,j) + theta(k,i)*f
               end do
            end if
            pcount = mod(pcount+1,np_j)
         end do
      end do
!$OMP BARRIER
      call D1dB_Vector_SumAll(2*Nqs*nfft3d,ufunc)

      return
      end 



*     ************************************************
*     *                                              *
*     *         vdw_DF_Generate_potentials           *
*     *                                              *
*     ************************************************
*
      subroutine vdw_DF_Generate_potentials(Nqs,nfft3d,ispin,n2ft3d,
     >                            ufunc,
     >                            q0,drho,ddrho,tmpexc,tmpfn,tmpfdn,
     >                            exc,fn,fdn)
      implicit none
      integer Nqs,nfft3d,ispin,n2ft3d
      real*8  ufunc(n2ft3d,Nqs)
      real*8  q0(n2ft3d)
      real*8  drho(n2ft3d,ispin)
      real*8  ddrho(n2ft3d,ispin)
      real*8  tmpexc(n2ft3d)
      real*8  tmpfn(n2ft3d,ispin),tmpfdn(n2ft3d,ispin)
      real*8  exc(n2ft3d),fn(n2ft3d,ispin),fdn(n2ft3d,ispin)

*     **** local variables ****
      integer i,j,jstart,nj,taskid_j,np_j,tid,nthr
      integer nx,r,ms
      real*8  pj,dpj

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

      nx = Nqs/np_j
      r = mod(Nqs,np_j)
      if (taskid_j.lt.r) then
         nj = nx + 1
         jstart = nx*taskid_j + taskid_j + 1
      else
         nj = nx
         jstart = nx*taskid_j + r + 1
      end if
      call D3dB_r_Zero(1,tmpexc)
      call D3dB_r_nZero(1,ispin,tmpfn)
      call D3dB_r_nZero(1,ispin,tmpfdn)

      if (nj.gt.0) then

         !*** compute ufunc(r) ***
         call Grsm_gh_fftb(nfft3d,nj,ufunc(1,jstart))

         !*** generate tmpexc,tmpfn,tmpfdn ***
         do i=tid+1,n2ft3d,nthr
            do j=jstart,jstart-1+nj
               call vdw_DF_poly(j,q0(i),pj,dpj)
               tmpexc(i) = tmpexc(i) + 0.5d0*pj*ufunc(i,j)
               do ms=1,ispin
                  tmpfn(i,ms) = tmpfn(i,ms) 
     >                        + (pj + dpj*drho(i,ms))*ufunc(i,j)
                  tmpfdn(i,ms)=tmpfdn(i,ms) + dpj*ddrho(i,ms)*ufunc(i,j)
               end do
            end do
         end do

      end if
!$OMP BARRIER
      call D1dB_Vector_SumAll(n2ft3d,tmpexc)
      call D1dB_Vector_SumAll(ispin*n2ft3d,tmpfn)
      call D1dB_Vector_SumAll(ispin*n2ft3d,tmpfdn)

      call D3dB_rr_Sum2(1,tmpexc,exc)
      do ms=1,ispin
         call D3dB_rr_Sum2(1,tmpfn(1,ms),fn(1,ms))
         call D3dB_rr_Sum2(1,tmpfdn(1,ms),fdn(1,ms))
      end do

      return
      end
