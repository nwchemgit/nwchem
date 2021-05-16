c
c     $Id$
c

*     *********************************
*     *                               *
*     *   velocity_autocorrelation    *
*     *                               *
*     *********************************
*
*   This routine calculated the velocity auto correlation
* using the MOTION file.  
*
*  - This routine needs to be rewritten so that it
*    doesn't use so much memory
*
*  Updated: 5-8-2002

      subroutine velocity_autocorrelation(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


      real*8 c,auTOsec
      parameter (c=2.996d+10,auTOsec=2.41889d-17)
 
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer i,k,natoms,natoms1,nframes,iii
      integer vx(2),vy(2),vz(2)
      integer vt(2),vk(2)
      integer tfft(2)
      real*8  time1,time2,t,freq,dt,volume,volume1
      real*8  tmax
      real*8  x,y,z
      real*8  v1,v2,v3
      real*8  unita(3,3)
      character*80  tfilename
      character*255 full_filename1
      character*255 full_filename2
      character*2 tag2
      character*4 tag4


*     *** synchronization ***
      call ga_sync()

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:velocity_correlation_filename',
     >                   1,tfilename))
     >   call util_file_prefix('velocity_correlation',tfilename)

      call util_file_name_noprefix(tfilename,.false.,
     >                             .false.,
     >                             full_filename1)

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:power_spectrum_filename',
     >                   1,tfilename))
     >   call util_file_prefix('power_spectrum',tfilename)
      call util_file_name_noprefix(tfilename,.false.,
     >                             .false.,
     >                             full_filename2)

      call Parallel_taskid(taskid)
     
*     **********************************************
*     ****** Velocity autocorrelation routine ******
*     *********************************************
      if (taskid.eq.MASTER) then
     
      write(*,*)
      write(*,*) '===== Calculating velocity auto correlation ===='
      write(*,*)    
      call MOTION_nframes(filename,nframes,natoms,volume,unita)

      open(unit=19,file=filename,form='formatted')

*     **** try to exit in a nice way if too not enough memory available ****
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vx',vx(2),vx(1))) 
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vy',vy(2),vy(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value = BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',1,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vz',vz(2),vz(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',2,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vt',vt(2),vt(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',3,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vk',vk(2),vk(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',4,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
       if (.not.BA_push_get(mt_dcpl,(8*nframes+15),
     >    'tfft',tfft(2),tfft(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vk(2))
         value = value.and.BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',5,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      

      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vx(1)),1)
      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vy(1)),1)
      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vz(1)),1)
      call dcopy(4*nframes,       0.0d0,0,dcpl_mb(vt(1)),1)
      call dcopy(4*nframes,        0.0d0,0,dcpl_mb(vk(1)),1)
     
      time1   = 0.0d0
      time2   = 0.0d0
      do k=1,nframes
        
         time1    = time2
         read(19,*) time2,natoms1,volume1,unita
         do i=1,natoms
            read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3
            dcpl_mb(vx(1)+(k-1) + (i-1)*2*nframes) = v1
            dcpl_mb(vy(1)+(k-1) + (i-1)*2*nframes) = v2
            dcpl_mb(vz(1)+(k-1) + (i-1)*2*nframes) = v3        
         end do

      end do
      close(unit=19)

      dt   = time2-time1
      write(6,*)
      write(6,*) 'frames used: ',nframes
      write(6,*) 'natoms:      ',natoms
      write(6,*) 'volume:      ',volume
      write(6,*) 'time interval (au)       : ', dt
      write(6,*) 'Max time      (au)       : ', dt*(nframes-1)
      write(6,*) 'frequency interval (cm-1): ', 
     >           1.0d0/(dt*nframes*auTOsec*c)
      write(6,*) 'Max frequency (cm-1)     : ', 
     >           1.0d0/(2*dt*auTOsec*c)
      write(6,*) 'MOTION filename: ',filename
      write(6,*) '<v*v>  filename: ',full_filename1
      write(6,*) 's(k)   filename: ',full_filename2
      write(6,*)

      dt   = time2-time1
      tmax = dt*(nframes-1)
      call dcffti(2*nframes,dcpl_mb(tfft(1)))
      do i=1,natoms
         call dcfftf(2*nframes,dcpl_mb(vx(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
         call dcfftf(2*nframes,dcpl_mb(vy(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
         call dcfftf(2*nframes,dcpl_mb(vz(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
      end do

      do k=1,(2*nframes)
         dcpl_mb(vt(1)+(k-1)) = dcmplx(0.0d0,0.0d0)
         do i=1,natoms
           dcpl_mb(vt(1)+(k-1)) = dcpl_mb(vt(1)+(k-1))
     >           + dconjg(dcpl_mb(vx(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vx(1)+(k-1)+(i-1)*2*nframes)
     >           + dconjg(dcpl_mb(vy(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vy(1)+(k-1)+(i-1)*2*nframes)
     >           + dconjg(dcpl_mb(vz(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vz(1)+(k-1)+(i-1)*2*nframes)

         end do
      end do
      x = 1.0d0/dble(2*nframes)
      call dscal(4*nframes,x,dcpl_mb(vt(1)),1)
      call dcfftb(2*nframes,dcpl_mb(vt(1)),dcpl_mb(tfft(1)))

      !**** average changes because there are less number of samples for higher t ****
      do k=1,nframes
         x = 1.0d0/dble(nframes-(k-1))
         dcpl_mb(vt(1)+k-1) = dcpl_mb(vt(1)+k-1)*x
      end do

      !**** normalize the velocity autocorrelation function by <v(0)*v(0)> ****
      x =  dble(dcpl_mb(vt(1)))
      x = 1.0d0/x
      do k=1,nframes
         dcpl_mb(vt(1)+k-1) = dcpl_mb(vt(1)+k-1)*x
      end do

      

*     ***** write out velocity autocorrelation ****
      open(unit=19,file=full_filename1,form='formatted')
      do k=1,nframes
         t = dt*(k-1)
c         x = 1.0d0/dble(nframes-(k-1))
c         dcpl_mb(vt(1)+k-1) = dcpl_mb(vt(1)+k-1)*x
         write(19,*) t,dble(dcpl_mb(vt(1)+k-1)),
     >                dimag(dcpl_mb(vt(1)+k-1))
      end do
      close(unit=19)


*     ************************************************
*     ****                                        ****
*     **** calculate the power spectrum using fft ****
*     ****                                        ****
*     ************************************************
      call dcffti(nframes,dcpl_mb(tfft(1)))
      call dcfftf(nframes,dcpl_mb(vt(1)),dcpl_mb(tfft(1)))
      !x = 1.0d0/dble(nframes)
      dcpl_mb(vk(1)) = dconjg(dcpl_mb(vt(1)))
     >                *dcpl_mb(vt(1))
      dcpl_mb(vk(1)+nframes/2) = dconjg(dcpl_mb(vt(1)+nframes/2))
     >                          *dcpl_mb(vt(1)+nframes/2)
      do k=1,(nframes/2)-1
        dcpl_mb(vk(1)+k) =   dconjg(dcpl_mb(vt(1)+k))
     >                      *dcpl_mb(vt(1)+k)
     >                   +   dconjg(dcpl_mb(vt(1)+nframes-k))
     >                      *dcpl_mb(vt(1)+nframes-k)
      end do

      !*** scale the power spectrum ***
c      x = 0.0d0
c      do k=0,nframes/2
c        x = x + dble(dcpl_mb(vk(1)+k))
c      end do
c      write(*,*) "powernorm=",x,x/(dt*nframes),
c     >           " dw=",1.0d0/(dt*nframes*auTOsec*c)

c      if (x.gt.1.0d-12) then
c         x = 1.0d0/x 
c      else 
c         x = 1.0d0
c      end if
      x = auTOsec*c
      call dscal(4*nframes,x,dcpl_mb(vk(1)),1)

      open(unit=19,file=full_filename2,form='formatted')
      do k=0,nframes/2
         if(k.ne.0) then
            t = nframes*dt/dble(k)
            freq = 1.0d0/(t*auTOsec*c)
         else
            freq=0.0d0
         endif
         write(19,*) freq, dble(dcpl_mb(vk(1)+k)),
     >                    dimag(dcpl_mb(vk(1)+k))
      end do
      close(unit=19)


      value =           BA_pop_stack(tfft(2))
      value = value.and.BA_pop_stack(vk(2))
      value = value.and.BA_pop_stack(vt(2))
      value = value.and.BA_pop_stack(vz(2))
      value = value.and.BA_pop_stack(vy(2))
      value = value.and.BA_pop_stack(vx(2))
      if (.not. value) 
     >    call errquit('cpmd_properties: cant release stack memory',0,
     >       MA_ERR)
     
      end if

*     *** synchronization ***
  999 continue
      call ga_sync()

      return
      end



*     **************************************
*     *                                    *
*     *   dipole_motion_autocorrelation    *
*     *                                    *
*     **************************************
*
*   This routine calculated the dipole auto correlation
* using the MOTION file.  
*
*  - This routine needs to be rewritten so that it
*    doesn't use so much memory
*
*  Updated: 5-8-2002

      subroutine dipole_motion_autocorrelation(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


      real*8 c,auTOsec,CAconv
      parameter (c=2.996d+10,auTOsec=2.41889d-17)
      parameter (CAconv=3.058900977e9)
 
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer i,k,natoms,natoms1,nframes,iii
      integer vx(2),vy(2),vz(2),dx(2),dy(2),dz(2)
      integer vt(2),vk(2)
      integer tfft(2)
      real*8  time1,time2,t,freq,dt,volume,volume1
      real*8  tmax
      real*8  x,y,z
      real*8  v1,v2,v3,sv1,sv2,sv3
      real*8  unita(3,3),omega
      character*80  tfilename
      character*255 full_filename1
      character*255 full_filename2

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      omega = lattice_omega()

*     *** synchronization ***
      call ga_sync()

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:dipole_motion_correlation_filename',
     >                   1,tfilename))
     >   call util_file_prefix('dipole_correlation',tfilename)

      call util_file_name_noprefix(tfilename,.false.,.false.,
     >                             full_filename1)

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:dipole_motion_powerspectrum_filename',
     >                   1,tfilename))
     >   call util_file_prefix('dipole_powerspectrum',tfilename)
      call util_file_name_noprefix(tfilename,.false.,.false.,
     >                             full_filename2)

      call Parallel_taskid(taskid)
     
*     **********************************************
*     ****** Velocity autocorrelation routine ******
*     *********************************************
      if (taskid.eq.MASTER) then
     
      write(*,*)
      write(*,*) '===== Calculating dipole auto correlation ===='
      write(*,*)    
      call dipole_motion_nframes(filename,nframes)

      open(unit=36,file=filename,form='formatted')

*     **** try to exit in a nice way if too not enough memory available ****
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vx',vx(2),vx(1))) 
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole motion autocorrelation not performed'
         close(unit=36)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vy',vy(2),vy(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole motion autocorrelation not performed'
         value = BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',1,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vz',vz(2),vz(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole motion autocorrelation not performed'
         value =           BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',2,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if

      call dcopy(4*nframes,0.0d0,0,dcpl_mb(vx(1)),1)
      call dcopy(4*nframes,0.0d0,0,dcpl_mb(vy(1)),1)
      call dcopy(4*nframes,0.0d0,0,dcpl_mb(vz(1)),1)
     
      !**** allocate dx,dy,dz ****
      if (.not.BA_push_get(mt_dbl,nframes,'dx',dx(2),dx(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole autocorrelation not performed'
         value =           BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',6,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
      if (.not.BA_push_get(mt_dbl,nframes,'dy',dy(2),dy(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole autocorrelation not performed'
         value =           BA_pop_stack(dx(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value)
     >    call errquit('cpmd_properties: cannot release stack memory',7,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
      if (.not.BA_push_get(mt_dbl,nframes,'dz',dz(2),dz(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole autocorrelation not performed'
         value =           BA_pop_stack(dy(2))
         value = value.and.BA_pop_stack(dx(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value)
     >    call errquit('cpmd_properties: cannot release stack memory',8,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if


      time1   = 0.0d0
      time2   = 0.0d0
      sv1 = 0.0d0
      sv2 = 0.0d0
      sv3 = 0.0d0
      do k=1,nframes
        
         time1    = time2
         read(36,*) time2,v1,v2,v3
         dbl_mb(dx(1)+(k-1)) = v1
         dbl_mb(dy(1)+(k-1)) = v2
         dbl_mb(dz(1)+(k-1)) = v3        
         sv1 = sv1 + v1
         sv2 = sv2 + v2
         sv3 = sv3 + v3
      end do
      close(unit=36)
      sv1 = sv1/dble(nframes)
      sv2 = sv2/dble(nframes)
      sv3 = sv3/dble(nframes)
      do k=1,nframes
         dbl_mb(dx(1)+k-1) = dbl_mb(dx(1)+k-1) - sv1
         dbl_mb(dy(1)+k-1) = dbl_mb(dy(1)+k-1) - sv2
         dbl_mb(dz(1)+k-1) = dbl_mb(dz(1)+k-1) - sv3
      end do

      !**** compute the dipole velocities ****
      dt   = time2-time1
      v1 = (dbl_mb(dx(1)+1) - dbl_mb(dx(1)))/(dt)
      v2 = (dbl_mb(dy(1)+1) - dbl_mb(dy(1)))/(dt)
      v3 = (dbl_mb(dz(1)+1) - dbl_mb(dz(1)))/(dt)
      dcpl_mb(vx(1)) = dcmplx(v1,0.0d0)
      dcpl_mb(vy(1)) = dcmplx(v2,0.0d0)
      dcpl_mb(vz(1)) = dcmplx(v3,0.0d0)
      do k=2,nframes-1
         v1 = (dbl_mb(dx(1)+k) - dbl_mb(dx(1)+k-2))/(2.0d0*dt)
         v2 = (dbl_mb(dy(1)+k) - dbl_mb(dy(1)+k-2))/(2.0d0*dt)
         v3 = (dbl_mb(dz(1)+k) - dbl_mb(dz(1)+k-2))/(2.0d0*dt)
         dcpl_mb(vx(1)+k-1) = dcmplx(v1,0.0d0)
         dcpl_mb(vy(1)+k-1) = dcmplx(v2,0.0d0)
         dcpl_mb(vz(1)+k-1) = dcmplx(v3,0.0d0)
      end do
      v1 = (dbl_mb(dx(1)+nframes-1) - dbl_mb(dx(1)+nframes-2))/(dt)
      v2 = (dbl_mb(dy(1)+nframes-1) - dbl_mb(dy(1)+nframes-2))/(dt)
      v3 = (dbl_mb(dz(1)+nframes-1) - dbl_mb(dz(1)+nframes-2))/(dt)
      dcpl_mb(vx(1)+nframes-1) = dcmplx(v1,0.0d0)
      dcpl_mb(vy(1)+nframes-1) = dcmplx(v2,0.0d0)
      dcpl_mb(vz(1)+nframes-1) = dcmplx(v3,0.0d0)

      !**** deallocate dx,dy,dz ****
      value =           BA_pop_stack(dz(2))
      value = value.and.BA_pop_stack(dy(2))
      value = value.and.BA_pop_stack(dx(2))
      if (.not. value)
     > call errquit('cpmd_properties: cannot release stack memory',9,
     >              MA_ERR)
   
      !**** allocate vt,vk,tfft ****
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vt',vt(2),vt(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole motion autocorrelation not performed'
         value =           BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value)
     >    call errquit('cpmd_properties: cannot release stack memory',3,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vk',vk(2),vk(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole autocorrelation not performed'
         value =           BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value)
     >    call errquit('cpmd_properties: cannot release stack memory',4,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
       if (.not.BA_push_get(mt_dcpl,(8*nframes+15),
     >    'tfft',tfft(2),tfft(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + dipole autocorrelation not performed'
         value =           BA_pop_stack(vk(2))
         value = value.and.BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value)
     >    call errquit('cpmd_properties: cannot release stack memory',5,
     >                 MA_ERR)
         close(unit=36)
         goto 999
      end if
      call dcopy(4*nframes,0.0d0,0,dcpl_mb(vt(1)),1)
      call dcopy(4*nframes, 0.0d0,0,dcpl_mb(vk(1)),1)



      write(6,*)
      write(6,*) 'frames used: ',nframes
      write(6,*) 'time interval (au)       : ', dt
      write(6,*) 'Max time      (au)       : ', dt*(nframes-1)
      write(6,*) 'frequency interval (cm-1): ', 
     >           1.0d0/(dt*nframes*auTOsec*c)
      write(6,*) 'Max frequency (cm-1)     : ', 
     >           1.0d0/(2*dt*auTOsec*c)
      write(6,*) 'dipole_motion filename: ',filename
      write(6,*) '<mu*mu>       filename: ',full_filename1
      write(6,*) 's(k)          filename: ',full_filename2
      write(6,*)

      dt   = time2-time1
      tmax = dt*(nframes-1)
      call dcffti(2*nframes,dcpl_mb(tfft(1)))

      call dcfftf(2*nframes,dcpl_mb(vx(1)),dcpl_mb(tfft(1)))
      call dcfftf(2*nframes,dcpl_mb(vy(1)),dcpl_mb(tfft(1)))
      call dcfftf(2*nframes,dcpl_mb(vz(1)),dcpl_mb(tfft(1)))

      do k=1,(2*nframes)
         dcpl_mb(vt(1)+(k-1)) = 
     >           + dconjg(dcpl_mb(vx(1)+(k-1)))
     >                   *dcpl_mb(vx(1)+(k-1))
     >           + dconjg(dcpl_mb(vy(1)+(k-1)))
     >                   *dcpl_mb(vy(1)+(k-1))
     >           + dconjg(dcpl_mb(vz(1)+(k-1)))
     >                   *dcpl_mb(vz(1)+(k-1))
      end do
      x = 1.0d0/dble(2*nframes)
      call dscal(4*nframes,x,dcpl_mb(vt(1)),1)
      call dcfftb(2*nframes,dcpl_mb(vt(1)),dcpl_mb(tfft(1)))


*     ***** write out dipole motion autocorrelation ****
      open(unit=36,file=full_filename1,form='formatted')
      do k=1,nframes
         t = dt*(k-1)
         x = 1.0d0/dble(nframes-(k-1))
         dcpl_mb(vt(1)+k-1) = dcpl_mb(vt(1)+k-1)*x
         write(36,*) t,dble(dcpl_mb(vt(1)+k-1)),
     >                dimag(dcpl_mb(vt(1)+k-1))
      end do
      close(unit=36)


*     ************************************************
*     ****                                        ****
*     **** calculate the power spectrum using fft ****
*     ****                                        ****
*     ************************************************
      call dcffti(nframes,dcpl_mb(tfft(1)))
      call dcfftf(nframes,dcpl_mb(vt(1)),dcpl_mb(tfft(1)))
      !x = 1.0d0/dble(nframes)
      dcpl_mb(vk(1)) = dconjg(dcpl_mb(vt(1)))
     >                *dcpl_mb(vt(1))
      dcpl_mb(vk(1)+nframes/2) = dconjg(dcpl_mb(vt(1)+nframes/2))
     >                          *dcpl_mb(vt(1)+nframes/2)
      do k=1,(nframes/2)-1
        dcpl_mb(vk(1)+k) =   dconjg(dcpl_mb(vt(1)+k))
     >                      *dcpl_mb(vt(1)+k)
     >                   +   dconjg(dcpl_mb(vt(1)+nframes-k))
     >                      *dcpl_mb(vt(1)+nframes-k)
      end do

c      !*** scale the power spectrum ***
c      x = 0.0d0
c      do k=0,nframes/2
c        x = x + dble(dcpl_mb(vk(1)+k))
c      end do
c      if (x.gt.1.0d-12) then
c         x = 1.0d0/x 
c      else 
c         x = 1.0d0
c      end if
c      call dscal(4*nframes,x,dcpl_mb(vk(1)),1)
c      x = auTOsec*c
c      call dscal(4*nframes,x,dcpl_mb(vk(1)),1)

      open(unit=36,file=full_filename2,form='formatted')
      do k=0,nframes/2
         if(k.ne.0) then
            t = nframes*dt/dble(k)
            freq = 1.0d0/(t*auTOsec*c)
         else
            freq=0.0d0
         endif
         write(36,*) freq, dble(dcpl_mb(vk(1)+k)),
     >                    dimag(dcpl_mb(vk(1)+k)),
     >                     dble(dcpl_mb(vk(1)+k))*CAconv/omega,
     >                    dimag(dcpl_mb(vk(1)+k))*CAconv/omega
      end do
      close(unit=36)

      value =           BA_pop_stack(tfft(2))
      value = value.and.BA_pop_stack(vk(2))
      value = value.and.BA_pop_stack(vt(2))
      value = value.and.BA_pop_stack(vz(2))
      value = value.and.BA_pop_stack(vy(2))
      value = value.and.BA_pop_stack(vx(2))
      if (.not. value) 
     >    call errquit('cpmd_properties: cant release stack memory',0,
     >       MA_ERR)
     
      end if

*     *** synchronization ***
  999 continue
      call ga_sync()

      return
      end



*     *******************************************
*     *                                         *
*     *       dipole_motion_nframes             *
*     *                                         *
*     *******************************************
      subroutine dipole_motion_nframes(filename,nframes)
      implicit none
      character*(*) filename
      integer nframes

      real*8    time,dipole(3)

      nframes = 0

      open(unit=36,file=filename,form='formatted',status='old')
      do while(.true.)
        read(36,*,ERR=30,END=30) time,dipole
        nframes = nframes + 1
      end do
 30   continue
      close(unit=36)

      return
      end




*     *********************************
*     *                               *
*     *       velocity_spectrum       *
*     *                               *
*     *********************************
*
*   This routine calculated the velocity spectrum         
* using the MOTION file.  
*
*  - This routine needs to be rewritten so that it
*    doesn't use so much memory
*

      subroutine velocity_spectrum(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"

      real*8 c,auTOsec
      parameter (c=2.996d+10,auTOsec=2.41889d-17)
 
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer i,k,natoms,natoms1,nframes,iii
      integer vx(2),vy(2),vz(2)
      integer vt(2),vk(2)
      integer tfft(2)
      real*8  time1,time2,t,freq,dt,volume,volume1
      real*8  tmax,mass
      real*8  x,y,z
      real*8  v1,v2,v3
      real*8  unita(3,3)
      character*80  tfilename
      !character*255 full_filename1
      character*255 full_filename2
      character*2 tag2
      character*4 tag4

*     **** external functions ****
      integer  ion_katm
      external ion_katm
      real*8   ion_amass
      external ion_amass


*     *** synchronization ***
      call ga_sync()

c      call util_file_name_noprefix('VEL_CORRELATION',.false.,
c     >                             .false.,
c     >                             full_filename1)

      if (.not.btdb_cget(rtdb,
     >                   'nwpw:velocity_spectrum_filename',
     >                   1,tfilename))
     >   call util_file_prefix('velocity_spectrum',tfilename)
      call util_file_name_noprefix(tfilename,.false.,
     >                             .false.,
     >                             full_filename2)

      call Parallel_taskid(taskid)
     
*     ***************************************
*     ****** Velocity spectrum routine ******
*     ***************************************
      if (taskid.eq.MASTER) then
     
      write(*,*)
      write(*,*) '===== Calculating velocity spectrum ====='
      write(*,*)    
      call MOTION_nframes(filename,nframes,natoms,volume,unita)

      open(unit=19,file=filename,form='formatted')

*     **** try to exit in a nice way if too not enough memory available ****
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vx',vx(2),vx(1))) 
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vy',vy(2),vy(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value = BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',1,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*natoms*nframes,'vz',vz(2),vz(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',2,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vt',vt(2),vt(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',3,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      if (.not.BA_push_get(mt_dcpl,2*nframes,'vk',vk(2),vk(1)))
     >   then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',4,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
       if (.not.BA_push_get(mt_dcpl,(8*nframes+15),
     >    'tfft',tfft(2),tfft(1))) then
         write(*,*) 'cpmd_properties: out of stack memory'
         write(*,*) ' + velocity autocorrelation not performed'
         value =           BA_pop_stack(vk(2))
         value = value.and.BA_pop_stack(vt(2))
         value = value.and.BA_pop_stack(vz(2))
         value = value.and.BA_pop_stack(vy(2))
         value = value.and.BA_pop_stack(vx(2))
         if (.not. value) 
     >    call errquit('cpmd_properties: cannot release stack memory',5,
     >                 MA_ERR)
         close(unit=19)
         goto 999
      end if
      

      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vx(1)),1)
      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vy(1)),1)
      call dcopy(4*natoms*nframes,0.0d0,0,dcpl_mb(vz(1)),1)
      call dcopy(4*nframes,       0.0d0,0,dcpl_mb(vt(1)),1)
      call dcopy(4*nframes,        0.0d0,0,dcpl_mb(vk(1)),1)
     
      time1   = 0.0d0
      time2   = 0.0d0
      do k=1,nframes
        
         time1    = time2
         read(19,*) time2,natoms1,volume1,unita
         do i=1,natoms
            read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3
            dcpl_mb(vx(1)+(k-1) + (i-1)*2*nframes) = v1
            dcpl_mb(vy(1)+(k-1) + (i-1)*2*nframes) = v2
            dcpl_mb(vz(1)+(k-1) + (i-1)*2*nframes) = v3        
         end do

      end do
      close(unit=19)

      dt   = time2-time1
      write(luout,*)
      write(luout,*) 'frames used: ',nframes
      write(luout,*) 'natoms:      ',natoms
      write(luout,*) 'volume:      ',volume
      write(luout,*) 'time interval (au)       : ', dt
      write(luout,*) 'Max time      (au)       : ', dt*(nframes-1)
      write(luout,*) 'frequency interval (cm-1): ', 
     >           1.0d0/(dt*nframes*auTOsec*c)
      write(luout,*) 'Max frequency (cm-1)     : ', 
     >           1.0d0/(2*dt*auTOsec*c)
      write(luout,*) 'MOTION filename: ',filename
      write(luout,*) 'S(k)   filename: ',full_filename2
      write(luout,*)

      dt   = time2-time1
      tmax = dt*(nframes-1)
      call dcffti(2*nframes,dcpl_mb(tfft(1)))
      do i=1,natoms
         call dcfftf(2*nframes,dcpl_mb(vx(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
         call dcfftf(2*nframes,dcpl_mb(vy(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
         call dcfftf(2*nframes,dcpl_mb(vz(1)+ (i-1)*2*nframes),
     >                         dcpl_mb(tfft(1)))
      end do

      do k=1,(2*nframes)
         dcpl_mb(vt(1)+(k-1)) = dcmplx(0.0d0,0.0d0)
         do i=1,natoms
           mass = ion_amass(ion_katm(i))
           dcpl_mb(vt(1)+(k-1)) = dcpl_mb(vt(1)+(k-1)) +
     >           mass*(
     >           + dconjg(dcpl_mb(vx(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vx(1)+(k-1)+(i-1)*2*nframes)
     >           + dconjg(dcpl_mb(vy(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vy(1)+(k-1)+(i-1)*2*nframes)
     >           + dconjg(dcpl_mb(vz(1)+(k-1)+(i-1)*2*nframes))
     >                   *dcpl_mb(vz(1)+(k-1)+(i-1)*2*nframes))

         end do
      end do
      x = 1.0d0/dble(2*nframes)
      call dscal(4*nframes,x,dcpl_mb(vt(1)),1)


*     ************************************************
*     ****                                        ****
*     ****    generate the velocity spectrum      ****
*     ****                                        ****
*     ************************************************
      open(unit=19,file=full_filename2,form='formatted')
      !do k=0,nframes/2
      do k=0,nframes
         if(k.ne.0) then
            t = 2*nframes*dt/dble(k)
            freq = 1.0d0/(t*auTOsec*c)
         else
            freq=0.0d0
         endif
         write(19,*) freq, dble(dcpl_mb(vt(1)+k))
      end do
      close(unit=19)


      value =           BA_pop_stack(tfft(2))
      value = value.and.BA_pop_stack(vk(2))
      value = value.and.BA_pop_stack(vt(2))
      value = value.and.BA_pop_stack(vz(2))
      value = value.and.BA_pop_stack(vy(2))
      value = value.and.BA_pop_stack(vx(2))
      if (.not. value) 
     >    call errquit('cpmd_properties: cant release stack memory',0,
     >       MA_ERR)
     
      end if

*     *** synchronization ***
  999 continue
      call ga_sync()

      return
      end

