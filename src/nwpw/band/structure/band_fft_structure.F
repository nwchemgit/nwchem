      subroutine band_fft_structure(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"


*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

      integer i,j,k,ms,n
      integer ne(2),ne_excited(2),neall(2)
      integer eigs(2),beigs(2),tmp(2),ntmp
      integer mgrid(3),ngrid(3)
      real*8  ks(3),kv(3),ee
      character*255 full_filename
 

*     **** external functions ****
      logical  cpsi_spin_orbit,control_spin_orbit,control_print
      integer  cpsi_ispin,cpsi_ne,control_num_kvectors_structure
      integer  control_ksvector_index
      integer  control_excited_ne
      real*8   lattice_unitg
      external cpsi_spin_orbit,control_spin_orbit,control_print
      external cpsi_ispin,cpsi_ne,control_num_kvectors_structure
      external control_ksvector_index
      external control_excited_ne
      external lattice_unitg


      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (cpsi_ispin().eq.2) ne_excited(2) = control_excited_ne(2)
      neall(1) = ne(1) + ne_excited(1)
      neall(2) = ne(2) + ne_excited(2)

      call control_monkhorst_pack_grid(mgrid)
      call control_fft_dos_grid(ngrid)

      if ((mgrid(1).ne.0) .and. (ngrid(1).ne.0) .and.
     >    (mgrid(2).ne.0) .and. (ngrid(2).ne.0) .and.
     >    (mgrid(3).ne.0) .and. (ngrid(3).ne.0)) then

      call Parallel_taskid(taskid)

*     **** allocate eigs ****
      if (.not.BA_push_get(mt_dcpl,mgrid(1)*mgrid(2)*mgrid(3),
     >                     'eigs',eigs(2),eigs(1)))
     >  call errquit('band_interpolate structure:out of stack memory',
     >               0,MA_ERR)

      if (.not.BA_push_get(mt_dcpl,ngrid(1)*ngrid(2)*ngrid(3),
     >                     'beigs',beigs(2),beigs(1)))
     >  call errquit('band_interpolate structure:out of stack memory',
     >               0,MA_ERR)

      ntmp = ngrid(1)*ngrid(2)*ngrid(3)+5*(ngrid(1)+ngrid(2)+ngrid(3))
      if (.not.BA_push_get(mt_dcpl,ntmp,'tmp',tmp(2),tmp(1)))
     >  call errquit('band_interpolate structure:out of stack memory',
     >               0,MA_ERR)

*     **** produce eigenvalue band file(s) ****
      call util_file_name('scatter_band_fft_interpolate',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
      end if


      if (taskid.eq.MASTER) then
         write(luout,*) "FFT interpolate band structure calculation" 
         write(luout,*) "mgrid=",mgrid
         write(luout,*) "ngrid=",ngrid
      end if



      do ms=1,cpsi_ispin()
      do n=1,neall(ms)

         call band_fft_interpolate_fetch_mp(mgrid,ne,ms,n,
     >                                      dcpl_mb(eigs(1)))

         call band_fft_interpolate_expand(mgrid,dcpl_mb(eigs(1)),
     >                                    ngrid,dcpl_mb(beigs(1)),
     >                                    ntmp,dcpl_mb(tmp(1)))
         do k=0,ngrid(3)-1
         do j=0,ngrid(2)-1
         do i=0,ngrid(1)-1
           ks(1) = dble(i)/dble(ngrid(1))
           ks(2) = dble(j)/dble(ngrid(2))
           ks(3) = dble(k)/dble(ngrid(3))
           kv(1) = lattice_unitg(1,1)*ks(1)
     >           + lattice_unitg(1,2)*ks(2)
     >           + lattice_unitg(1,3)*ks(3)
           kv(2) = lattice_unitg(2,1)*ks(1)
     >           + lattice_unitg(2,2)*ks(2)
     >           + lattice_unitg(2,3)*ks(3)
           kv(3) = lattice_unitg(3,1)*ks(1)
     >           + lattice_unitg(3,2)*ks(2)
     >           + lattice_unitg(3,3)*ks(3)
           ee = dcpl_mb(beigs(1) + i 
     >                           + j*ngrid(1) 
     >                           + k*ngrid(1)*ngrid(2))

           if (taskid.eq.MASTER) write(58,*) ms,n,i,j,k,ks,kv,ee
         end do
         end do
         end do

      end do
      end do
       
      if (taskid.eq.MASTER) then
        close(58)
      end if

      if (.not.BA_pop_stack(tmp(2))) 
     > call errquit('band_fft_interpolate:popping stack',0,MA_ERR)

      if (.not.BA_pop_stack(beigs(2))) 
     >  call errquit('band_fft_interpolate:popping stack',0,MA_ERR)

      if (.not.BA_pop_stack(eigs(2))) 
     >  call errquit('band_fft_interpolate:popping stack',0,MA_ERR)


      end if

      return
      end 

*     ********************************************************
*     *                                                      *
*     *             band_fft_interpolate_fetch_mp            *
*     *                                                      *
*     ********************************************************
      subroutine band_fft_interpolate_fetch_mp(mgrid,ne,ms,n,eigs)
      implicit none
      integer mgrid(3),ne(2),ms,n
      complex*16 eigs(mgrid(1),mgrid(2),mgrid(3))

*     **** local variables ****
      integer i,j,k,nb
      real*8 xx,yy,zz,xxx,yyy,zzz,ks(3),f

*     **** external functions ****
      integer   control_ksvector_index
      external  control_ksvector_index
      real*8    cpsi_eigenvalue_brdcst,cpsi_virtual_brdcst
      externale cpsi_eigenvalue_brdcst,cpsi_virtual_brdcst

      xxx = 1.0d0/(2.0d0*mgrid(1))
      yyy = 1.0d0/(2.0d0*mgrid(2))
      zzz = 1.0d0/(2.0d0*mgrid(3))

      do k=0,mgrid(3)-1
      do j=0,mgrid(2)-1
      do i=0,mgrid(1)-1
         xx = (1.0d0 + 2*i - mgrid(1))
         yy = (1.0d0 + 2*j - mgrid(2))
         zz = (1.0d0 + 2*k - mgrid(3))
         ks(1) = xx*xxx
         ks(2) = yy*yyy
         ks(3) = zz*zzz
         nb = control_ksvector_index(ks)
         if (n.le.ne(ms)) then
            f = cpsi_eigenvalue_brdcst(nb,ms,n)
         else
            f = cpsi_virtual_brdcst(nb,ms,n-ne(ms))
         end if
         eigs(i+1,j+1,k+1) = f
      end do
      end do
      end do

      return
      end

*     ***************************************************************
*     *                                                             *
*     *                band_fft_interpolate_expand                  *
*     *                                                             *
*     ***************************************************************
      subroutine band_fft_interpolate_expand(mgrid,fsmall,ngrid,fbig,
     >                                       ntmp,tmp)
      implicit none
      integer mgrid(3)
      complex*16 fsmall(mgrid(1),mgrid(2),mgrid(3))
      integer ngrid(3)
      complex*16 fbig(ngrid(1),ngrid(2),ngrid(3))
      integer ntmp
      complex*16 tmp(*)

*     **** local variables ****
      integer i,j,k,i2,j2,k2,idiff,jdiff,kdiff
      real*8  xshift(3),scal
      complex*16 phfac(3)

      scal = 1.0d0/(mgrid(1)*mgrid(2)*mgrid(3))
      
      do i=1,3
         !xshift(i) = 0.0d0
         !xshift(i) = -0.5d0*(1.0d0-1.0d0/dble(mgrid(i)))
         !xshift(i) = 0.5d0*(1.0d0/dble(mgrid(i)))
         !xshift(i) = -0.5d0*(1.0d0-1.0d0/dble(mgrid(i))) - best?
         !xshift(i) = 0.5d0*(1.0d0-1.0d0/dble(mgrid(i)))
         !xshift(i) = 0.5d0*(1.0d0)
         !xshift(i) = -0.5d0*(1.0d0+1.0d0/dble(mgrid(i)))
         !xshift(i) = -0.5d0*(1.0d0)
c         xshift(i) = 0.5d0*(1.0d0-1.0d0/dble(mgrid(i))
c     >                     +1.0d0/dble(ngrid(i)))
         xshift(i) = -(0.5d0 - 0.5d0/dble(mgrid(i)))
      end do

      call nwpw_cc_fft3f(fsmall,
     >                    mgrid(1),mgrid(1)*mgrid(2),
     >                    mgrid(1),mgrid(2),mgrid(3),
     >                    tmp,ntmp)

      call band_fft_interpolate_translate(xshift,mgrid,fsmall,
     >                                    tmp,
     >                                    tmp(mgrid(1)+1),
     >                                    tmp(mgrid(1)+mgrid(2)+2))

      idiff = ngrid(1) - mgrid(1)
      jdiff = ngrid(2) - mgrid(2)
      kdiff = ngrid(3) - mgrid(3)

      call dcopy(2*ngrid(1)*ngrid(2)*ngrid(3),0.0d0,0,fbig,1)

      do k=0,mgrid(3)-1
      do j=0,mgrid(2)-1
      do i=0,mgrid(1)-1
         if (k.lt.mgrid(3)/2) then
            k2 = k
         else
            k2 = kdiff + k
         end if
         if (j.lt.mgrid(2)/2) then
            j2 = j
         else
            j2 = jdiff + j
         end if
         if (i.lt.mgrid(1)/2) then
            i2 = i
         else
            i2 = idiff + i
         end if
         fbig(i2+1,j2+1,k2+1) = fsmall(i+1,j+1,k+1)*scal
      end do
      end do
      end do

      call nwpw_cc_fft3b(fbig,
     >                   ngrid(1),ngrid(1)*ngrid(2),
     >                   ngrid(1),ngrid(2),ngrid(3),
     >                   tmp,ntmp)

      return
      end 


*     *****************************************************
*     *                                                   *
*     *           band_fft_interpolate_translate          *
*     *                                                   *
*     *****************************************************
      subroutine band_fft_interpolate_translate(trans,ngrid,f,
     >                                          ex1,ex2,ex3)
      implicit none
      real*8 trans(3)
      integer ngrid(3)
      complex*16 f(ngrid(1),ngrid(2),ngrid(3))
      complex*16 ex1(*),ex2(*),ex3(*)
   
*     **** local variables ****
      integer i,j,k,taskid
      real*8     pi,sw1,sw2,sw3
      complex*16 cw1,cw2,cw3

      pi = 4.0d0*datan(1.0d0)
      call Parallel_taskid(taskid)
      
      !sw1 = 2.0d0*pi*trans(1) + pi
      !sw2 = 2.0d0*pi*trans(2) + pi
      !sw3 = 2.0d0*pi*trans(3) + pi
      sw1 = 2.0d0*pi*trans(1)
      sw2 = 2.0d0*pi*trans(2)
      sw3 = 2.0d0*pi*trans(3)
      cw1 = dcmplx(dcos(sw1),-dsin(sw1))
      cw2 = dcmplx(dcos(sw2),-dsin(sw2))
      cw3 = dcmplx(dcos(sw3),-dsin(sw3))
      ex1(1) = dcmplx(1.0d0,0.0d0)
      ex2(1) = dcmplx(1.0d0,0.0d0)
      ex3(1) = dcmplx(1.0d0,0.0d0)
      do i=1,ngrid(1)/2
         ex1(i+1)          = ex1(i)*cw1
         ex1(ngrid(1)-i+1) = dconjg(ex1(i+1))
      end do
      !if (mod(ngrid(1),2).eq.0) ex1(ngrid(1)/2+1)=dcmplx(0.0d0,0.0d0)

      do j=1,ngrid(2)/2
         ex2(j+1)          = ex2(j)*cw2
         ex2(ngrid(2)-j+1) = dconjg(ex2(j+1))
      end do
      !if (mod(ngrid(2),2).eq.0) ex2(ngrid(2)/2+1)=dcmplx(0.0d0,0.0d0)

      do k=1,ngrid(3)/2
         ex3(k+1)          = ex3(k)*cw3
         ex3(ngrid(3)-k+1) = dconjg(ex3(k+1))
      end do
      !if (mod(ngrid(3),2).eq.0) ex3(ngrid(3)/2+1)=dcmplx(0.0d0,0.0d0)
      if (taskid.eq.0) then
      write(*,*) "ex1(ngrid(1)/2)   = ", ex1(ngrid(1)/2) 
      write(*,*) "ex1(ngrid(1)/2+1) = ", ex1(ngrid(1)/2+1) 

      write(*,*) "ex2(ngrid(2)/2)   = ", ex1(ngrid(2)/2) 
      write(*,*) "ex2(ngrid(2)/2+1) = ", ex1(ngrid(2)/2+1) 

      write(*,*) "ex3(ngrid(3)/2)   = ", ex1(ngrid(3)/2) 
      write(*,*) "ex3(ngrid(3)/2+1) = ", ex1(ngrid(3)/2+1) 
      end if

      do k=1,ngrid(3)
      do j=1,ngrid(2)
      do i=1,ngrid(1)
         f(i,j,k) = f(i,j,k)*ex1(i)*ex2(j)*ex3(k)
      end do
      end do
      end do

      return
      end 
