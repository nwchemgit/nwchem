      subroutine band_interpolate_structure(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"


*     **** parallel variables ****
      integer  taskid,taskid_k,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

      integer nbrillall,k,isize,NN,i,j,n1,n2,n3,ii,i1,i2,i3
      integer ne(2),ne_excited(2),neall(2)
      integer itet(3,10),eigs(2),mgrid(3),tindx(10),ix,iy,iz
      real*8  ktet(3,10),kstet(3,10),kstetm(3,10),rx,ry,rz
      real*8  ks(3),kv(3),kold(3),pathlength,dist,e1,e2,EV
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
      NN=neall(1)-neall(2)
      EV=27.2116d0

      nbrillall = control_num_kvectors_structure()
      call control_monkhorst_pack_grid(mgrid)

      if ((nbrillall.gt.1).and.
     >    (mgrid(1).ne.0) .and.
     >    (mgrid(2).ne.0) .and.
     >    (mgrid(3).ne.0)) then

      call Parallel_taskid(taskid)

*     **** allocate eigs ****
      isize = neall(1)+neall(2)
      if (.not.BA_push_get(mt_dbl,isize,'eigs',eigs(2),eigs(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',0,MA_ERR)



*     **** produce eigenvalue band file(s) ****
      if (cpsi_ispin().eq.1) then
        call util_file_name('restricted_band_interpolated',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('restricted_band_interpolated_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if

      else
        if (cpsi_spin_orbit()) then
        call util_file_name('spinor_band_interpolated',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('spinor_band_interpolated_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if


        else 
        call util_file_name('alpha_band_interpolated',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band_interpolated',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=59,file=full_filename,form='formatted')
        end if
        call util_file_name('alpha_band_interpolated_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band_interpolated_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=69,file=full_filename,form='formatted')
        end if


        end if
      end if

      if (taskid.eq.MASTER) 
     >   write(luout,*) "interpolated band structure calculation" 

      do k=1,nbrillall
         call control_ksvector_structure(k,ks)
         call control_kvector_structure(k,kv)
         rx = 0.5d0*(ks(1)*2*mgrid(1) + mgrid(1) - 1)
         ry = 0.5d0*(ks(2)*2*mgrid(2) + mgrid(2) - 1)
         rz = 0.5d0*(ks(3)*2*mgrid(3) + mgrid(3) - 1)
         i1 = nint(rx)
         i2 = nint(ry)
         i3 = nint(rz)
         if (rx.ge.i1)  then
            ix= 1
         else
            ix=-1
         end if
         if (ry.ge.i2)  then
            iy= 1
         else
            iy=-1
         end if
         if (rz.ge.i3)  then
            iz= 1
         else
            iz=-1
         end if
         itet(1,1) = i1
         itet(2,1) = i2
         itet(3,1) = i3

         itet(1,2) = i1 + ix
         itet(2,2) = i2
         itet(3,2) = i3

         itet(1,3) = i1
         itet(2,3) = i2 + iy
         itet(3,3) = i3

         itet(1,4) = i1
         itet(2,4) = i2
         itet(3,4) = i3 + iz

         itet(1,5) = i1 - ix
         itet(2,5) = i2
         itet(3,5) = i3

         itet(1,6) = i1
         itet(2,6) = i2 - iy
         itet(3,6) = i3

         itet(1,7) = i1
         itet(2,7) = i2
         itet(3,7) = i3 - iz

         itet(1,8) = i1 + ix
         itet(2,8) = i2 + iy
         itet(3,8) = i3

         itet(1,9) = i1 + ix
         itet(2,9) = i2 
         itet(3,9) = i3 + iz

         itet(1,10) = i1 
         itet(2,10) = i2 + iy
         itet(3,10) = i3 + iz

c         itet(1,11) = i1 + ix
c         itet(2,11) = i2 - iy
c         itet(3,11) = i3
c         itet(1,12) = i1 - ix
c         itet(2,12) = i2 + iy
c         itet(3,12) = i3 
c         itet(1,13) = i1 - ix
c         itet(2,13) = i2 - iy
c         itet(3,13) = i3
c         itet(1,14) = i1 + ix
c         itet(2,14) = i2 
c         itet(3,14) = i3 - iz
c         itet(1,15) = i1 - ix
c         itet(2,15) = i2 
c         itet(3,15) = i3 + iz
c         itet(1,16) = i1 - ix
c         itet(2,16) = i2 
c         itet(3,16) = i3 - iz
c         itet(1,17) = i1 
c         itet(2,17) = i2 + iy
c         itet(3,17) = i3 - iz
c         itet(1,18) = i1 
c         itet(2,18) = i2 - iy
c         itet(3,18) = i3 + iz
c         itet(1,19) = i1 
c         itet(2,19) = i2 - iy
c         itet(3,19) = i3 - iz
c         itet(1,20) = i1 + ix
c         itet(2,20) = i2 + iy
c         itet(3,20) = i3 + iz


         do j=1,10
         do i=1,3
            kstet(i,j) = (2*itet(i,j)+1-mgrid(i))/(2.0d0*mgrid(i))
         end do
         end do

         do j=1,10
         do i=1,3
            ktet(i,j) = lattice_unitg(i,1)*kstet(1,j)
     >                + lattice_unitg(i,2)*kstet(2,j)
     >                + lattice_unitg(i,3)*kstet(3,j)
         end do
         end do

         do j=1,10
            do i=1,3
               itet(i,j)   = mod(itet(i,j),mgrid(i))
               kstetm(i,j) = (2*itet(i,j)+1-mgrid(i))/(2.0d0*mgrid(i))
            end do
            tindx(j) = control_ksvector_index(kstetm(1,j))
         end do

         if (taskid.eq.MASTER) then
            if (k.eq.1) then
               pathlength = 0.0d0
            else
               dist =        (kv(1)-kold(1))**2
               dist = dist + (kv(2)-kold(2))**2
               dist = dist + (kv(3)-kold(3))**2
               dist = dsqrt(dist)
               pathlength = pathlength + dist
            end if
            kold(1) = kv(1)
            kold(2) = kv(2)
            kold(3) = kv(3)
            write(luout,1508) k,pathlength,ks,kv
c            do i=1,4
c               write(luout,1509) i,tindx(i),
c     >                             kstet(1,i),kstet(2,i),kstet(3,i),
c     >                              ktet(1,i), ktet(2,i), ktet(3,i)
c            end do
           
            write(luout,1500)
         end if

         call band_interpolate_eigs(kv,tindx,ktet,dbl_mb(eigs(1)))

         !*** spin-orbit ****
         if (taskid.eq.MASTER) then
            if (cpsi_spin_orbit()) then
               do i=0,neall(1)-1
                  e1 = dbl_mb(eigs(1)+neall(1)-i-1)
                  write(luout,1510)  e1,e1*EV
               end do
            !*** not spin-orbit ****
            else
               do i=0,NN-1
                  e1 = dbl_mb(eigs(1)+neall(1)-i-1)
                  write(luout,1510) e1,e1*EV
               end do
               do i=0,neall(2)-1
                  e1 = dbl_mb(eigs(1)+neall(1)-i-NN-1)
                  e2 = dbl_mb(eigs(1)+neall(1)+neall(2)-i-1)
                  write(luout,1510)  e1,e1*EV,e2,e2*EV
               end do
            end if
            write(58,'(1000E14.6)') pathlength,
     >          (dbl_mb(eigs(1)+i-1),i=1,neall(1))
            write(68,'(1000E14.6)') pathlength,ks,kv,
     >          (dbl_mb(eigs(1)+i-1),i=1,neall(1))
            if ((.not.cpsi_spin_orbit()).and.(cpsi_ispin().eq.2)) then
               write(59,'(1000E14.6)') pathlength,
     >             (dbl_mb(eigs(1)+i-1),i=1,neall(2))
               write(69,'(1000E14.6)') pathlength,ks,kv,
     >             (dbl_mb(eigs(1)+i-1),i=1,neall(2))
            end if
         end if

      end do

      if (taskid.eq.MASTER) then
        close(58)
        close(68)
        if (cpsi_ispin().eq.2) then
           close(59)
           close(69)
        end if
      end if

      if (.not.BA_pop_stack(eigs(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)


      end if

 1500 FORMAT(/' interpolated orbital energies:')
 1508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1509 FORMAT(/' Brillouin zone tet-point: ',i3,i5
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))


      return
      end 

      subroutine band_interpolate_eigs(kv,tindx,ktet,eigs)
      implicit none
      real*8  kv(3)
      integer tindx(*)
      real*8  ktet(3,*)
      real*8  eigs(*)

#include "bafdecls.fh"
#include "errquit.fh"

      integer  taskid,MASTER
      parameter(MASTER=0)

*     **** local variables ****
      logical done,noarea,novolume
      integer nbrill,nb,i,imax,n1,n2,n3,ierr,i1,i2,i3,i4
      integer dist(2),kvv(2),indx(2),itet(4),ipiv(10)
      integer ne(2),ne_excited(2),neall(2)
      real*8  vandermonde(10,10),A(10,10),f(10),kmin(3),kf(3)
      real*8 dd,dp,dm,dx,dy,dz
      real*8 xp,yp,zp,xm,ym,zm

*     **** external functions ****
      logical  cpsi_spin_orbit
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      integer  control_excited_ne
      real*8   brillioun_k_brdcst,lattice_unitg,cpsi_eigenvalue_brdcst
      real*8   cpsi_virtual_brdcst
      external cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external control_excited_ne
      external brillioun_k_brdcst,lattice_unitg,cpsi_eigenvalue_brdcst
      external cpsi_virtual_brdcst

      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (cpsi_ispin().eq.2) ne_excited(2) = control_excited_ne(2)
      neall(1) = ne(1) + ne_excited(1)
      neall(2) = ne(2) + ne_excited(2)
      
      call band_interpolate_setvandermonde(ktet,vandermonde)

      do i=1,neall(1)
         do i1 = 1,10
            if (i.le.ne(1)) then
               f(i1) = cpsi_eigenvalue_brdcst(tindx(i1),1,i)
            else
               f(i1) = cpsi_virtual_brdcst(tindx(i1),1,i-ne(1))
            end if
         end do
         call dcopy(100,vandermonde,1,A,1)
         call DGESV(10,1,A,10,ipiv,f,10,ierr)
         eigs(i) = f(1)*kv(1)+f(2)*kv(2)+f(3)*kv(3)
     >           + f(4)*kv(1)*kv(1)
     >           + f(5)*kv(1)*kv(2)
     >           + f(6)*kv(1)*kv(3)
     >           + f(7)*kv(2)*kv(2)
     >           + f(8)*kv(2)*kv(3)
     >           + f(9)*kv(3)*kv(3)
     >           + f(10)
c     >           + f(10)*kv(1)*kv(1)*kv(1)
c     >           + f(11)*kv(1)*kv(1)*kv(2)
c     >           + f(12)*kv(1)*kv(1)*kv(3)
c     >           + f(13)*kv(1)*kv(2)*kv(2)
c     >           + f(14)*kv(1)*kv(2)*kv(3)
c     >           + f(15)*kv(1)*kv(3)*kv(3)
c     >           + f(16)*kv(2)*kv(2)*kv(2)
c     >           + f(17)*kv(2)*kv(2)*kv(3)
c     >           + f(18)*kv(2)*kv(3)*kv(3)
c     >           + f(19)*kv(3)*kv(3)*kv(3)
c     >           + f(20)
      end do

      if ((.not.cpsi_spin_orbit()).and.(cpsi_ispin().eq.2)) then
         do i=1,neall(2)
            do i1 = 1,10
               if (i.le.ne(2)) then
                  f(i1) = cpsi_eigenvalue_brdcst(tindx(i1),2,i)
               else
                  f(i1) = cpsi_virtual_brdcst(tindx(i1),2,i-ne(2))
               end if
            end do
            call dcopy(100,vandermonde,1,A,1)
            call DGESV(10,1,A,10,ipiv,f,10,ierr)
            eigs(i+neall(1))=f(1)*kv(1)+f(2)*kv(2)+f(3)*kv(3)
     >           + f(4)*kv(1)*kv(1)
     >           + f(5)*kv(1)*kv(2)
     >           + f(6)*kv(1)*kv(3)
     >           + f(7)*kv(2)*kv(2)
     >           + f(8)*kv(2)*kv(3)
     >           + f(9)*kv(3)*kv(3)
     >           + f(10)
c     >           + f(10)*kv(1)*kv(1)*kv(1)
c     >           + f(11)*kv(1)*kv(1)*kv(2)
c     >           + f(12)*kv(1)*kv(1)*kv(3)
c     >           + f(13)*kv(1)*kv(2)*kv(2)
c     >           + f(14)*kv(1)*kv(2)*kv(3)
c     >           + f(15)*kv(1)*kv(3)*kv(3)
c     >           + f(16)*kv(2)*kv(2)*kv(2)
c     >           + f(17)*kv(2)*kv(2)*kv(3)
c     >           + f(18)*kv(2)*kv(3)*kv(3)
c     >           + f(19)*kv(3)*kv(3)*kv(3)
c     >           + f(20)
         end do
      end if         

      return
      end



      subroutine band_interpolate_setvandermonde(ktet,vandermonde)
      implicit none
      real*8 ktet(3,*),vandermonde(10,10)

      integer i

      do i=1,10
         vandermonde(i,1) = ktet(1,i)
         vandermonde(i,2) = ktet(2,i)
         vandermonde(i,3) = ktet(3,i)
         vandermonde(i,4) = ktet(1,i)*ktet(1,i)
         vandermonde(i,5) = ktet(1,i)*ktet(2,i)
         vandermonde(i,6) = ktet(1,i)*ktet(3,i)
         vandermonde(i,7) = ktet(2,i)*ktet(2,i)
         vandermonde(i,8) = ktet(2,i)*ktet(3,i)
         vandermonde(i,9) = ktet(3,i)*ktet(3,i)
         vandermonde(i,10) = 1.0d0

c         vandermonde(i,10) = ktet(1,i)*ktet(1,i)*ktet(1,i)
c         vandermonde(i,11) = ktet(1,i)*ktet(1,i)*ktet(2,i)
c         vandermonde(i,12) = ktet(1,i)*ktet(1,i)*ktet(3,i)
c         vandermonde(i,13) = ktet(1,i)*ktet(2,i)*ktet(2,i)
c         vandermonde(i,14) = ktet(1,i)*ktet(2,i)*ktet(3,i)
c         vandermonde(i,15) = ktet(1,i)*ktet(3,i)*ktet(3,i)
c         vandermonde(i,16) = ktet(2,i)*ktet(2,i)*ktet(2,i)
c         vandermonde(i,17) = ktet(2,i)*ktet(2,i)*ktet(3,i)
c         vandermonde(i,18) = ktet(2,i)*ktet(3,i)*ktet(3,i)
c         vandermonde(i,19) = ktet(3,i)*ktet(3,i)*ktet(3,i)
c         vandermonde(i,20) = 1.0d0
      end do

      return
      end

c $Id$
