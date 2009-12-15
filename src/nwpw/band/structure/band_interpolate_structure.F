      subroutine band_interpolate_structure(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"


*     **** parallel variables ****
      integer  taskid,taskid_k,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

      integer nbrillall,k,isize,NN,i
      integer eigs(2)
      real*8  ks(3),kv(3),kold(3),pathlength,dist,e1,e2,EV
      character*255 full_filename
 

*     **** external functions ****
      logical     cpsi_spin_orbit,control_spin_orbit,control_print
      integer     cpsi_ispin,cpsi_ne,control_num_kvectors_structure
      external    cpsi_spin_orbit,control_spin_orbit,control_print
      external    cpsi_ispin,cpsi_ne,control_num_kvectors_structure


      NN=cpsi_ne(1)-cpsi_ne(2)
      EV=27.2116d0

      nbrillall = control_num_kvectors_structure()
      if (nbrillall.gt.0) then
      call Parallel_taskid(taskid)

*     **** allocate eigs ****
      isize = (cpsi_ne(1)+cpsi_ne(2))
      if (.not.MA_push_get(mt_dbl,isize,'eigs',eigs(2),eigs(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',0,MA_ERR)



*     **** produce eigenvalue band file(s) ****
      if (cpsi_ispin().eq.1) then
        call util_file_name('restricted_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
      else
        if (cpsi_spin_orbit()) then
        call util_file_name('spinor_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        else 
        call util_file_name('alpha_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=59,file=full_filename,form='formatted')
        end if
        end if
      end if

      if (taskid.eq.MASTER) 
     >   write(luout,*) "band structure calculation" 

      do k=1,nbrillall
         call control_ksvector_structure(k,ks)
         call control_kvector_structure(k,kv)

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
            write(luout,1500)
         end if

         call band_interpolate_eigs(kv,dbl_mb(eigs(1)))

         !*** spin-orbit ****
         if (taskid.eq.MASTER) then
            if (cpsi_spin_orbit()) then
               do i=0,cpsi_ne(1)-1
                  e1 = dbl_mb(eigs(1)+cpsi_ne(1)-i-1)
                  write(luout,1510)  e1,e1*EV
               end do
            !*** not spin-orbit ****
            else
               do i=0,NN-1
                  e1 = dbl_mb(eigs(1)+cpsi_ne(1)-i-1)
                  write(luout,1510) e1,e1*EV
               end do
               do i=0,cpsi_ne(2)-1
                  e1 = dbl_mb(eigs(1)+cpsi_ne(1)-i-NN-1)
                  e2 = dbl_mb(eigs(1)+cpsi_ne(1)+cpsi_ne(2)-i-1)
                  write(luout,1510)  e1,e1*EV,e2,e2*EV
               end do
            end if
            write(58,'(1000E14.6)') pathlength,
     >          (dbl_mb(eigs(1)+i-1),i=1,cpsi_ne(1))
            if ((.not.cpsi_spin_orbit()).and.(cpsi_ispin().eq.2)) then
               write(59,'(1000E14.6)') pathlength,
     >             (dbl_mb(eigs(1)+i-1),i=1,cpsi_ne(2))
            end if
         end if

      end do

      if (taskid.eq.MASTER) then
        close(58)
        if (cpsi_ispin().eq.2) close(59)
      end if

      if (.not.MA_pop_stack(eigs(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)


      end if

 1500 FORMAT(/' interpolated orbital energies:')
 1508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))


      return
      end 

      subroutine band_interpolate_eigs(kv,eigs)
      implicit none
      real*8 kv(3)
      real*8 eigs(*)

#include "mafdecls.fh"
#include "errquit.fh"

      integer  taskid,MASTER
      parameter(MASTER=0)

*     **** local variables ****
      logical done,noarea,novolume
      integer nbrill,nb,i,imax,n1,n2,n3,ierr,i1,i2,i3,i4
      integer dist(2),kvv(2),indx(2),itet(4),ipiv(4)
      real*8  vandermonde(4,4),A(4,4),f(4),kmin(3),kf(3)
      real*8 dd,dp,dm,dx,dy,dz
      real*8 xp,yp,zp,xm,ym,zm

*     **** external functions ****
      logical  cpsi_spin_orbit
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      real*8   brillioun_k_brdcst,lattice_unitg,cpsi_eigenvalue_brdcst
      external cpsi_spin_orbit
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external brillioun_k_brdcst,lattice_unitg,cpsi_eigenvalue_brdcst

      call Parallel_taskid(taskid)
      
      nbrill = brillioun_nbrillioun()

      if (.not.MA_push_get(mt_dbl,nbrill,'dist',dist(2),dist(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',0,MA_ERR)
      if (.not.MA_push_get(mt_dbl,3*nbrill,'kvv',kvv(2),kvv(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',1,MA_ERR)
      if (.not.MA_push_get(mt_int,nbrill,'indx',indx(2),indx(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',2,MA_ERR)


      do nb=1,nbrill
         kf(1) = brillioun_k_brdcst(1,nb)
         kf(2) = brillioun_k_brdcst(2,nb)
         kf(3) = brillioun_k_brdcst(3,nb)
      
         xp = kf(1)-kv(1)
         yp = kf(2)-kv(2)
         zp = kf(3)-kv(3)
         xm = -kf(1)-kv(1)
         ym = -kf(2)-kv(2)
         zm = -kf(3)-kv(3)

         dd = 1.0d10
         do n1=-1,1
         do n2=-1,1
         do n3=-1,1
            dx = n1*lattice_unitg(1,1)
     >         + n2*lattice_unitg(1,2)
     >         + n3*lattice_unitg(1,3)
            dy = n1*lattice_unitg(2,1)
     >         + n2*lattice_unitg(2,2)
     >         + n3*lattice_unitg(2,3)
            dz = n1*lattice_unitg(3,1)
     >         + n2*lattice_unitg(3,2)
     >         + n3*lattice_unitg(3,3)
            dp = dsqrt((xp+dx)**2 + (yp+dy)**2 + (zp+dz)**2)
            dm = dsqrt((xm+dx)**2 + (ym+dy)**2 + (zm+dz)**2)
            if (dp.lt.dd) then
               dd=dp
               kmin(1) = xp+dx
               kmin(2) = yp+dy
               kmin(3) = zp+dz
            end if
            if (dm.lt.dd) then
               dd=dm
               kmin(1) = xm+dx
               kmin(2) = ym+dy
               kmin(3) = zm+dz
            end if
         end do
         end do
         end do
         dbl_mb(dist(1)+nb-1)      = dd
         dbl_mb(kvv(1)+3*(nb-1))   = kmin(1)
         dbl_mb(kvv(1)+3*(nb-1)+1) = kmin(2)
         dbl_mb(kvv(1)+3*(nb-1)+2) = kmin(3)
         int_mb(indx(1)+nb-1) = nb
         if (taskid.eq.MASTER) write(*,*) "nb,kmin,dd=",nb,kmin,dd
      end do
      call band_interpolate_sort(nbrill,int_mb(indx(1)),dbl_mb(dist(1)))

      imax = 4
      itet(1) = 1
      itet(2) = 2
      itet(3) = 3
      itet(4) = 4
      done = .false.
      do while ((.not.done).and.(imax.le.nb))
         call band_interpolate_checktet(itet,int_mb(indx(1)),
     >                                  dbl_mb(kvv(1)),
     >                                  noarea,novolume)
         if (noarea) then
            itet(3) = imax+1
            imax    = imax+1
         else
            if (novolume) then
               itet(4) = imax+1
               imax    = imax+1
            else
               done = .true.
            end if
         end if
      end do
      if (done) then
         call band_interpolate_setvandermonde(itet,int_mb(indx(1)),
     >                                        dbl_mb(kvv(1)),
     >                                        vandermonde)

         i1 = int_mb(indx(1)+itet(1)-1) 
         i2 = int_mb(indx(1)+itet(2)-1) 
         i3 = int_mb(indx(1)+itet(3)-1) 
         i4 = int_mb(indx(1)+itet(4)-1) 
         do i=1,cpsi_ne(1)
            f(1) = cpsi_eigenvalue_brdcst(i1,1,i)
            f(2) = cpsi_eigenvalue_brdcst(i2,1,i)
            f(3) = cpsi_eigenvalue_brdcst(i3,1,i)
            f(4) = cpsi_eigenvalue_brdcst(i4,1,i)
            if (taskid.eq.MASTER) then
               write(*,*) "kv=",kv
               write(*,*) "f =",f
               call dcopy(16,vandermonde,1,A,1)
               write(*,*) "A=",A(1,1),A(1,2),A(1,3),A(1,4)
               write(*,*) "  ",A(2,1),A(2,2),A(2,3),A(2,4)
               write(*,*) "  ",A(3,1),A(3,2),A(3,3),A(3,4)
               write(*,*) "  ",A(4,1),A(4,2),A(4,3),A(4,4)

               call DGESV(4,1,A,4,ipiv,f,4,ierr)
               write(*,*) "ierr=",ierr
               eigs(i) = f(1)*kv(1)+f(2)*kv(2)+f(3)*kv(3)+f(4)
               write(*,*) "fout=",f
               write(*,*) "eig=",eigs(i)
            end if
         end do

         if ((.not.cpsi_spin_orbit()).and.(cpsi_ispin().eq.2)) then
            do i=1,cpsi_ne(2)
               f(1) = cpsi_eigenvalue_brdcst(i1,2,i)
               f(2) = cpsi_eigenvalue_brdcst(i2,2,i)
               f(3) = cpsi_eigenvalue_brdcst(i3,2,i)
               f(4) = cpsi_eigenvalue_brdcst(i4,2,i)
               call dcopy(16,vandermonde,1,A,1)
               call DGESV(4,1,A,4,ipiv,f,4,ierr)
               eigs(i+cpsi_ne(1))=f(1)*kv(1)+f(2)*kv(2)+f(3)*kv(3)+f(4)
            end do
         end if         

      else
         i1 = int_mb(indx(1)+itet(1)-1) 
         do i=1,cpsi_ne(1)
            eigs(i) = cpsi_eigenvalue_brdcst(i1,1,i)
         end do
         if ((.not.cpsi_spin_orbit()).and.(cpsi_ispin().eq.2)) then
            do i=1,cpsi_ne(2)
               eigs(i+cpsi_ne(1))= cpsi_eigenvalue_brdcst(i1,2,i)
            end do
         end if         
      end if

*     *** deallocate stack ***
      if (.not.MA_pop_stack(indx(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)
      if (.not.MA_pop_stack(kvv(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)
      if (.not.MA_pop_stack(dist(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)
     
      
      return
      end


      subroutine band_interpolate_checktet(itet,indx,point,
     >                                     noarea,novolume)
      implicit none
      integer itet(4),indx(*)
      real*8 point(3,*)
      logical noarea,novolume

      real*8 a(3,3),a1a2(3),vv

      a(1,1) = point(1,indx(itet(2)))-point(1,indx(itet(1)))
      a(2,1) = point(2,indx(itet(2)))-point(2,indx(itet(1)))
      a(3,1) = point(3,indx(itet(2)))-point(3,indx(itet(1)))
      a(1,2) = point(1,indx(itet(3)))-point(1,indx(itet(1)))
      a(2,2) = point(2,indx(itet(3)))-point(2,indx(itet(1)))
      a(3,2) = point(3,indx(itet(3)))-point(3,indx(itet(1)))
      a(1,3) = point(1,indx(itet(4)))-point(1,indx(itet(1)))
      a(2,3) = point(2,indx(itet(4)))-point(2,indx(itet(1)))
      a(3,3) = point(3,indx(itet(4)))-point(3,indx(itet(1)))
      a1a2(1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      a1a2(2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      a1a2(3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      vv = a(1,3)*a1a2(1)+a(2,3)*a1a2(2)+a(3,3)*a1a2(3)
      novolume = (dabs(vv).lt.1.0d-6)
      noarea   = (dabs(a1a2(1)).lt.1.0d-6).and.
     >           (dabs(a1a2(2)).lt.1.0d-6).and.
     >           (dabs(a1a2(3)).lt.1.0d-6)
      return
      end


      subroutine band_interpolate_setvandermonde(itet,indx,point,
     >                                           vandermonde)
      implicit none
      integer itet(4),indx(*)
      real*8 point(3,*),vandermonde(4,4)

      vandermonde(1,1) = point(1,indx(itet(1)))
      vandermonde(1,2) = point(2,indx(itet(1)))
      vandermonde(1,3) = point(3,indx(itet(1)))
      vandermonde(1,4) = 1.0d0
      vandermonde(2,1) = point(1,indx(itet(2)))
      vandermonde(2,2) = point(2,indx(itet(2)))
      vandermonde(2,3) = point(3,indx(itet(2)))
      vandermonde(2,4) = 1.0d0
      vandermonde(3,1) = point(1,indx(itet(3)))
      vandermonde(3,2) = point(2,indx(itet(3)))
      vandermonde(3,3) = point(3,indx(itet(3)))
      vandermonde(3,4) = 1.0d0
      vandermonde(4,1) = point(1,indx(itet(4)))
      vandermonde(4,2) = point(2,indx(itet(4)))
      vandermonde(4,3) = point(3,indx(itet(4)))
      vandermonde(4,4) = 1.0d0

      return
      end

      subroutine band_interpolate_sort(N,indx,key)
      implicit none
      integer N,indx(*)
      real*8 key(*)
      integer i,j,ii
      do i=1,N-1
      do j=i+1,N
         if (key(indx(j)).lt.key(indx(i))) then
           ii = indx(i)
           indx(i) = indx(j)
           indx(j) = ii
         end if
      end do
      end do
      return
      end
