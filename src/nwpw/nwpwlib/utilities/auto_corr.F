c
c     $Id: auto_corr.F,v 1.7 2007-04-02 15:38:38 bylaska Exp $
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

      subroutine velocity_autocorrelation(filename)
      implicit none
#include "errquit.fh"
      character*(*) filename

#include "mafdecls.fh"


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
      character*255 full_filename1
      character*255 full_filename2
      character*2 tag2
      character*4 tag4


*     *** synchronization ***
      call ga_sync()

      call util_file_name_noprefix('VEL_CORRELATION',.false.,
     >                             .false.,
     >                             full_filename1)
      call util_file_name_noprefix('POWER_SPECTRUM',.false.,
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

      value = MA_push_get(mt_dcpl,2*natoms*nframes,'vx',vx(2),vx(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*natoms*nframes,'vy',vy(2),vy(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*natoms*nframes,'vz',vz(2),vz(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*nframes,'vt',vt(2),vt(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*nframes,'vk',vk(2),vk(1))
       value = value.and.
     >        MA_push_get(mt_dcpl,(8*nframes+15),'tfft',tfft(2),tfft(1))
      
*     **** try to exit in a nice way ****
      if (.not. value) then
        write(*,*) 'cpmd_properties: out of stack memory'
        write(*,*) ' + velocity autocorrelation not performed'

        value = MA_pop_stack(tfft(2))
        value = MA_pop_stack(vk(2))
        value = MA_pop_stack(vt(2))
        value = MA_pop_stack(vz(2))
        value = MA_pop_stack(vy(2))
        value = MA_pop_stack(vx(2))

        return
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


*     ***** write out velocity autocorrelation ****
      open(unit=19,file=full_filename1,form='formatted')
      do k=1,nframes
         t = dt*(k-1)
         x = 1.0d0/dble(nframes-(k-1))
         dcpl_mb(vt(1)+k-1) = dcpl_mb(vt(1)+k-1)*x
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
      x = 0.0d0
      do k=0,nframes/2
        x = x + dble(dcpl_mb(vk(1)+k))
      end do
      if (x.gt.1.0d-12) then
         x = 1.0d0/x 
      else 
         x = 1.0d0
      end if
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

      value =           MA_pop_stack(tfft(2))
      value = value.and.MA_pop_stack(vk(2))
      value = value.and.MA_pop_stack(vt(2))
      value = value.and.MA_pop_stack(vz(2))
      value = value.and.MA_pop_stack(vy(2))
      value = value.and.MA_pop_stack(vx(2))
      if (.not. value) 
     >    call errquit('cpmd_properties: cant release stack memory',0,
     &       MA_ERR)
     
      end if

*     *** synchronization ***
      call ga_sync()

      return
      end

