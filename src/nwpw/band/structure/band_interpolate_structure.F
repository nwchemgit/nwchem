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

      integer nbrillall,k,isize
      integer eigs(2)
      real*8  ks(3),kv(3),kold(3),pathlength,dist
      character*255 full_filename
 

*     **** external functions ****
      logical     cpsi_spin_orbit,control_spin_orbit,control_print
      integer     cpsi_ispin,cpsi_ne,control_num_kvectors_structure
      external    cpsi_spin_orbit,control_spin_orbit,control_print
      external    cpsi_ispin,cpsi_ne,control_num_kvectors_structure


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
         if (k.eq.1) then
            pathlength = 0.0d0
         else
            dist=dsqrt((kv(1)-kold(1))**2
     >                +(kv(2)-kold(2))**2
     >                +(kv(3)-kold(3))**2)
            pathlength = pathlength + dist
         end if
         kold(1) = kv(1)
         kold(2) = kv(2)
         kold(3) = kv(3)


         if (taskid.eq.MASTER) then
            write(luout,1508) k,pathlength,ks,kv
            write(luout,1500)
         end if

         call band_interpolate_eigs(kv,dbl_mb(eigs(1)))

      end do


      if (taskid.eq.MASTER) then
        close(58)
        if (cpsi_ispin().eq.2) close(59)
      end if

      if (.not.MA_push_stack(eigs(2))) 
     >  call errquit(
     > 'band_interpolate structure:popping stack memory',0,MA_ERR)



      end if

 1500 FORMAT(/' interpolated orbital energies:')
 1508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')

      return
      end 

      subroutine band_interpolate_eigs(ks,kv,eigs)
      implicit none
      real*8 ks(3),kv(3)
      real*8 eigs(*)

      integer near,nb
      real*8 knear(3,10),dist(10),kf(3)

      
      if (.not.MA_push_get(mt_dbl,isize,'dist',dist(2),dist(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',0,MA_ERR)
      if (.not.MA_push_get(mt_int,isize,'indx',indx(2),indx(1)))
     >  call errquit(
     > 'band_interpolate structure:out of stack memory',0,MA_ERR)

      nbrill = brillioun_nbrillioun()

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
            if (dp.lt.dist) then
               dd=dp
               kmin(1) = xp+dx
               kmin(2) = yp+dy
               kmin(3) = zp+dz
            end if
            if (dm.lt.dist) then
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
      end do
      call band_interpolate_sort(nbrill,int_mb(indx(1)),dbl_mb(dist(1)))

      
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

      subroutine band_interpolate_prune(N,indx,key) 
      implicit none
      integer N,indx(*)
      real*8 key(*)

      return
      end
