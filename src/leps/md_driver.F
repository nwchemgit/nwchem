      logical function md_driver(rtdb,grad,thr,trajnum)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
 
      integer geom 
      integer natom
      integer trajnum
      integer md_data_get_i_c, md_data_get_i_m, md_data_get_i_q
      integer md_data_get_i_t, md_data_get_i_v, md_data_get_i_a
      integer md_data_get_i_b, md_data_get_i_d, md_data_get_i_f
      integer md_data_get_natom, md_data_get_nsteps
      double precision md_data_get_timestep, md_data_get_temp 
      logical status,md_data_set_temp,md_set

      logical grad
      external grad
      integer coor,mass,bq,tag,vel,acc,bterm,dterm,frc
      integer step, nsteps,naccept
      double precision dt,ek,temp 
      double precision gks_gasdev,gasdev,seed
      external gks_gasdev

      character*64 theory,blanks
      character*32 thr,task,blanks32
      data blanks(1:64) /' '/
      data blanks32(1:32)/' '/
      logical restart
      integer i 
      integer nxing,side_prev
      logical forward,backward
      logical mc_data_get_forward,mc_data_get_backward,mc_data_get_ifirc
      seed = -1234
      write(*,*) 'inside md_driver'
      gasdev=gks_gasdev(seed)


      restart=.false.
      theory=blanks
c      thr=blanks32
      task=':gradient'
      call concat(thr,task,theory) 
      dt = 0.5  ! time step
      ek=0.0 ! kinetic energy
      nsteps=200 

c      double precision eleps, md_driver1
      md_driver=.true.
      forward= mc_data_get_forward()
      backward=mc_data_get_backward()

      nsteps = md_data_get_nsteps()
      dt = md_data_get_timestep()   
      temp =md_data_get_temp()
      write(*,*) 'Md_driver:input params', nsteps, dt, temp
c
      if (.not.restart.and.mc_data_get_ifirc().and.forward) then
       write(*,*) 'Opening xyz and rest files on traj ',trajnum  
       open(unit=10,file='traj.xyz',status='new')
       open(unit=11,file='restart.pun',status='new')
      else
       open(unit=10,file='traj.xyz',status='old')
       open(unit=11,file='restart.pun',status='old')
c       rewind(unit=10)
c       rewind(unit=11)

      endif

cc     create geometry object 
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('leps_getgeom: failed creating geometry',
     & 0,GEOM_ERR)

C     Get the handles for everything
      natom=md_data_get_natom()
      coor=md_data_get_i_c()
      mass=md_data_get_i_m()
      bq=md_data_get_i_q()
      tag=md_data_get_i_t()
      vel=md_data_get_i_v()
      acc=md_data_get_i_a()
      bterm=md_data_get_i_b()
      dterm=md_data_get_i_d()
      frc=md_data_get_i_f()

c     end of initialization
c      write(*,*) 'Calling MD test:'
c      call md_test1(rtdb)
      write(*,*) 'Starting MD'
      if (restart) then
      write(*,*) 'Restarting from punch file'
             call read_punch (11,natom,nsteps,step,dt,dbl_mb(mass),
     &            dbl_mb(coor),dbl_mb(vel),dbl_mb(acc),dbl_mb(bterm),
     &            dbl_mb(dterm),byte_mb(tag))
      endif
      nxing=0
      side_prev=0

C-------- print the first point
      call print_coord(10,natom,dbl_mb(coor),byte_mb(tag))
C-------- this is where the MD loop starts

      do step=1,nsteps
       write(*,*) 'Predicting step :',step

         call predict(natom,dt,dbl_mb(coor),dbl_mb(vel),dbl_mb(acc),
     &             dbl_mb(bterm),dbl_mb(dterm))

C---- Add the predicted coordinates to the rtdb

       call md_data_geom_set(rtdb,geom)

C --- calculate forces for predicted coords             
       if (.not. grad(rtdb)) 
     & call errquit('md_gradient: no grad calculated ', 0, RTDB_ERR)

         if (.not. rtdb_get(rtdb,theory,mt_dbl,3*natom,
     &     dbl_mb(frc)))
     &   call errquit('md_gradient: no grad found ', 0, RTDB_ERR)

         write(*,*) 'Correcting step :',step

         call correct(natom,dt,dbl_mb(mass),dbl_mb(coor),dbl_mb(vel),
     &        dbl_mb(acc), dbl_mb(bterm),dbl_mb(dterm),dbl_mb(frc),ek)
         call print_coord(10,natom,dbl_mb(coor),byte_mb(tag)) 

         if (mod(step,10).eq.0.0) then
             rewind(11)  
             call punch (11,natom,nsteps,step,dt,dbl_mb(mass),
     &            dbl_mb(coor),dbl_mb(vel),dbl_mb(acc),dbl_mb(bterm),
     &            dbl_mb(dterm),byte_mb(tag))
         endif    
      enddo
C ---- End of MD_loop      
      call crossing(natom,dbl_mb(coor),trajnum,side_prev)

      if (.not. geom_destroy(geom))
     & call errquit('leps_getgeom:  geom_destroy failed',
     & 911,GEOM_ERR)
      
c     close(unit=10)
      if (backward) then 
        write(10,*) 'end of traj foreard/backward' 
      endif
      close(unit=11)
      write(*,*) 'END OF TRAJECTORY, side is :',side_prev
         return
      end

      subroutine predict (n,dt,r,v,a,b,c )
      implicit none
        integer     n
        double precision r(3,n), v(3,n), a(3,n), b(3,n), c(3,n)
        double precision dt
        integer     i,j
C ---- dt# are the Taylor expansion series terms dt1=dt, dt2=dt^2/2, dt3=dt^3/6 etc.
        real        dt1, dt2, dt3, dt4

c    *******************************************************************

c       write(*,*) 'In predictor:'
c       do i=1,n
c          do j=1,3
c           write (*,*) r(i,j),v(i,j),a(i,j),b(i,j),c(i,j)
c          enddo
c       enddo   
       dt1 = dt
       dt2 = dt1 * dt / 2.0
       dt3 = dt2 * dt / 3.0
       dt4 = dt3 * dt / 4.0
       do i = 1,n
        do j=1,3
         r(j,i) = r(j,i) + dt1*v(j,i) + dt2*a(j,i)
     1          +  dt3*b(j,i) + dt4*c(j,i)
         v(j,i) = v(j,i) + dt1*a(j,i) + dt2*b(j,i) + dt3*c(j,i)
         a(j,i) = a(j,i) + dt1*b(j,i) + dt2*c(j,i)
         b(j,i) = b(j,i) + dt1*c(j,i)
        enddo
       enddo
       return
      end
C ----- GEAR CORRECTOR
      subroutine correct (n,dt,mass,r,v,a,b,c,force,ek )
      implicit none
        integer     n
        double precision r(3,n),v(3,n),a(3,n),b(3,n),c(3,n),force(3,n)
        double precision dt
        integer     i,j
C ---- dt# are the Taylor expansion series terms dt1=dt, dt2=dt^2/2, dt3=dt^3/6 etc.
        double precision        dt1, dt2, dt3, dt4
        double precision       mass(n), ek
        double precision ai(3),corr(3)
        double precision  cr, cv, cb, cc
        double precision  c0, c1, c3, c4
        parameter ( c0 = 19.0 / 120.0, c1 = 3.0 / 4.0,
     1              c3 = 1.0 / 2.0,    c4 = 1.0 / 12.0 )

c    *******************************************************************

        dt1 = dt
        dt2 = dt1 * dt / 2.0
        dt3 = dt2 * dt / 3.0
        dt4 = dt3 * dt / 4.0

        cr = c0 * dt2
        cv = c1 * dt2 / dt1
        cb = c3 * dt2 / dt3
        cc = c4 * dt2 / dt4

        ek = 0.0

        do i = 1, n
         do j=1,3

           ai(j) = -force(j,i) / mass(i)
           corr(j) = ai(j) - a(j,i)

           r(j,i) = r(j,i) + cr * corr(j)
           v(j,i) = v(j,i) + cv * corr(j)
           a(j,i) = ai(j)
           b(j,i) = b(j,i) + cb * corr(j)
           c(j,i) = c(j,i) + cc * corr(j)

           ek = ek + v(j,i) ** 2
         enddo
       enddo

       ek = 0.5 * mass(i) * ek

       return
       end

       subroutine print_coord(rfile,n,c,t)
       implicit none
       integer n,i,rfile
       double precision c(3,n),cnv
       character*16 t(n)
       cnv=1.889725989
       if (rfile.eq.10) then
        write(rfile,*) '3'
        write(rfile,*) 
       endif  
       do i=1,n
          write(rfile,*) t(i),c(1,i)/cnv,c(2,i)/cnv,c(3,i)/cnv
       enddo  
       return
       end

      subroutine concat(str1,str2,str)
      implicit none
      character*32 str1
      character*32 str2
      character*64 str,blanks
      integer i,j,i1 
      data blanks(1:64) /' '/
C
      str = blanks
 
C Find the first nonblank in the file name.
C
      write(*,*) 'str1: ',str1
      write(*,*) 'str2: ',str2
      write(*,*) 'str: ',str
c      i1=1
c      j=1
c      Do while ((j .le. len(str1)) .and. (str1(j:j) .eq. ' '))
c        j=j+1
c      enddo
C
C Now concatenate the file name.
C
      i=1
     
      do while ((i .le. len(str1)) .and. (str1(i:i) .ne. ' '))
         str (i:i) = str1 (i:i)
c         i1 = i1 + 1
         i=i+1
      enddo
C Find the first nonblank in the file type.
C
c      j=1
c      Do while ((j .le. len(str2)) .and. (str2(j:j) .eq. ' '))
c        j=j+1
c      enddo
      i1=i 
      j=1
      do while ((j .le. len(str2)) .and. (str2(j:j) .ne. ' '))
            str (i1:i1) = str2 (j:j)
            i1 = i1 + 1
            j=j+1
      enddo

      write (*,*) 'Long string : ', str
      write (*,*)
      end
      subroutine punch (rfile,n,nsteps,step,dt,mass,c,v,a,b,d,t)
      implicit none
        integer     rfile,n,nsteps,step 
        double precision mass(n),c(3,n), v(3,n), a(3,n), b(3,n), d(3,n)
        double precision dt
        character*16 t(n)  
         write(rfile,*) n,nsteps,step,dt,mass,c,v,a,b,d,t
c        write(rfile,*) 'natom, md steps, curr. step, time step:'
c        write(rfile,*)  n,nsteps,step,dt
c        write(rfile,*) 'mass' 
c        do i =1,n
c          write(rfile,*) t(i),mass(i)
c        enddo
c        write(rfile,*) 'coordinates'
c        call print_coord(rfile,n,c,t)
c        write(rfile,*) 'velocities'
c        call print_coord(rfile,n,v,t)
c        write(rfile,*) 'accelerations'
c        call print_coord(rfile,n,a,t)
c        write(rfile,*) 'b-term'
c        call print_coord(rfile,n,b,t)
c        write(rfile,*) 'd-terms'
c        call print_coord(rfile,n,d,t)
c        write(rfile,*) 'End of restart file'
      end

      subroutine read_punch(rfile,n,nsteps,step,dt,mass,c,v,a,b,d,t)
      implicit none
        integer     rfile,n,nsteps,step
        double precision mass(n),c(3,n), v(3,n), a(3,n), b(3,n), d(3,n)
        double precision dt
        character*16 t(n)
         open(unit=12,file='restart.pun',status='old')
         read(12,*) n,nsteps,step,dt,mass,c,v,a,b,d,t
      end
        
