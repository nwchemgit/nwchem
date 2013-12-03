      logical function md_driver(rtdb,grad,thr)
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
#include "const_data.fh" 
      integer geom 
      integer natom
      integer trajnum
      integer md_data_get_i_c, md_data_get_i_m, md_data_get_i_q
      integer md_data_get_i_t, md_data_get_i_v, md_data_get_i_a
      integer md_data_get_i_b, md_data_get_i_d, md_data_get_i_f
      integer md_data_get_natom, md_data_get_nsteps,mc_data_get_trajnum
      double precision md_data_get_timestep, md_data_get_temp 

      logical grad, task_energy
      external grad, task_energy
      integer coor,mass,bq,tag,vel,acc,bterm,dterm,frc
      integer step, nsteps 
      double precision dt,ek,temp,E1,Etot 

      character*64 theorye,theoryg,blanks
      character*32 thr,taskg,taske
      data blanks(1:64) /' '/
      logical restart
      logical forward,backward
      logical mc_data_get_forward,mc_data_get_backward,mc_data_get_ifirc
      integer i, iw, iwtr
      logical dbug, status
      dbug=.false.
      iwtr=12
      iw=6
      md_driver=.false.

       write(iw,*) 'MD-DRIVER PRED-CORR: ',mc_data_get_trajnum()
      trajnum=mc_data_get_trajnum()
      restart=.false.
      theorye=blanks
      theoryg=blanks
      taskg=':gradient'
      taske=':energy'
      call concat(thr,taskg,theoryg) 
      call concat(thr,taske,theorye)
c      write(*,*) 'after concat'
      dt = 0.5  ! time step
      ek=0.0 ! kinetic energy
      nsteps=100 

      forward= mc_data_get_forward()
      backward=mc_data_get_backward()


      nsteps = md_data_get_nsteps()
      if (mc_data_get_ifirc()) nsteps=1
      dt = md_data_get_timestep()   
      temp =md_data_get_temp()
      write(iw,*) '\nInput parameters: \nNumber of steps, stepsize, 
     $ temperature:   \n',
     $   nsteps, dt, temp
       write (iw,*) 'Theory: ', thr
C  ----- may want to add MD restart, but not priority ----
      goto 100

      if (.not.restart.and.mc_data_get_ifirc().and.forward) then
       write(*,*) 'Opening xyz and rest files on traj ',trajnum  
       open(unit=10,file='traj.xyz',status='new')
       open(unit=11,file='restart.pun',status='new')
      else
#ifdef USE_OPEN_ACCESS
       open(unit=10,file='traj.xyz',status='old',ACCESS='append')
#else
       open(unit=10,file='traj.xyz',status='old',POSITION='append')
#endif
       open(unit=11,file='restart.pun',status='old')
      endif
  100 continue

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
c      write(*,*) 'Starting MD:'
      if (restart) then
      write(iw,*) 'Restarting from punch file'
             call read_punch (11,natom,nsteps,step,dt,dbl_mb(mass),
     &            dbl_mb(coor),dbl_mb(vel),dbl_mb(acc),dbl_mb(bterm),
     &            dbl_mb(dterm),byte_mb(tag))
      endif

C-------- print the first point
       if(dbug) then
         write(iw,*) ' Before predicting,coord:'
         call print_coord(iw,natom,dbl_mb(coor),byte_mb(tag))
         write(iw,*) ' vel:'
         call print_coord(iw,natom,dbl_mb(vel),byte_mb(tag))
         write(iw,*) ' acc:'
         call print_coord(iw,natom,dbl_mb(acc),byte_mb(tag))
         write(iw,*) ' b:'
         call print_coord(iw,natom,dbl_mb(bterm),byte_mb(tag))
         write(iw,*) ' d:'
         call print_coord(iw,natom,dbl_mb(dterm),byte_mb(tag))
       endif
C-------- i set initial params for higher derivatives 
      call md_init_trj(natom,dt,dbl_mb(mass),dbl_mb(coor),dbl_mb(vel),
     &     dbl_mb(acc), dbl_mb(bterm),dbl_mb(dterm),grad,rtdb,theoryg)
       if(dbug) then
         write(iw,*) ' Before predicting,after init,coord:'
         call print_coord(iw,natom,dbl_mb(coor),byte_mb(tag))
         write(iw,*) ' vel:'
         call print_coord(iw,natom,dbl_mb(vel),byte_mb(tag))
         write(iw,*) ' acc:'
         call print_coord(iw,natom,dbl_mb(acc),byte_mb(tag))
         write(iw,*) ' b:'
         call print_coord(iw,natom,dbl_mb(bterm),byte_mb(tag))
         write(iw,*) ' d:'
         call print_coord(iw,natom,dbl_mb(dterm),byte_mb(tag))
      endif
c      call  md_init_trj(n,dt,mass,r,v,a,b,c,grad,rtdb,theoryg)
C-------- this is where the MD loop starts
      do step=1,nsteps

         call predict1(natom,dt,dbl_mb(coor),dbl_mb(vel),dbl_mb(acc),
     &             dbl_mb(bterm),dbl_mb(dterm))

       if(dbug) then
         write(iw,*) 'Predicted coord:'
         call print_coord(iw,natom,dbl_mb(coor),byte_mb(tag))
         write(iw,*) 'Predicted vel:'
         call print_coord(iw,natom,dbl_mb(vel),byte_mb(tag))
         write(iw,*) 'Predicted acc:'
         call print_coord(iw,natom,dbl_mb(acc),byte_mb(tag))
         write(iw,*) ' b:'
         call print_coord(iw,natom,dbl_mb(bterm),byte_mb(tag))
         write(iw,*) ' d:'
         call print_coord(iw,natom,dbl_mb(dterm),byte_mb(tag))
      endif
         call md_data_geom_set(rtdb)
C --- calculate forces for predicted coords             
         if (.not. grad(rtdb)) 
     & call errquit('md_gradient: no grad calculated ', 0, RTDB_ERR)
         if (.not. rtdb_get(rtdb,theoryg,mt_dbl,3*natom,
     &     dbl_mb(frc)))
     &   call errquit('md_gradient: no grad found ', 0, RTDB_ERR)

         do i=1,3*natom
            if (abs(dbl_mb(frc+i-1)).lt.tolg) then
               dbl_mb(frc+i-1)=0.0
            endif
         enddo

         call correct1(natom,dt,dbl_mb(mass),dbl_mb(coor),dbl_mb(vel),
     &        dbl_mb(acc), dbl_mb(bterm),dbl_mb(dterm),dbl_mb(frc),ek)
       if (dbug) then 
         write(iw,*) 'Corrected coord:'
         call print_coord(iw,natom,dbl_mb(coor),byte_mb(tag))
         write(iw,*) 'Corrected vel:'
         call print_coord(iw,natom,dbl_mb(vel),byte_mb(tag))
         write(iw,*) 'Corrected acc:'
         call print_coord(iw,natom,dbl_mb(acc),byte_mb(tag))
         write(iw,*) ' b:'
         call print_coord(iw,natom,dbl_mb(bterm),byte_mb(tag))
         write(iw,*) ' d:'
         call print_coord(iw,natom,dbl_mb(dterm),byte_mb(tag))
       endif
c --- calculate energy for the corrected geometry         
         call md_data_geom_set(rtdb)

         status = task_energy(rtdb)
        
         if (.not. rtdb_get(rtdb,theorye,mt_dbl,1,E1))
     &   call errquit('md_driver:  no energy ', 0,
     &                 RTDB_ERR)
         Etot=ek+E1
         if (ga_nodeid().eq.0) then
           write(iwtr,1001) natom
           write(iwtr,1000) '@MD trajectory ',trajnum, ' MD step ',
     &   step,' V: ',E1,' Ek: ',ek,' Tot E: ',Etot
 
           call print_coord(iwtr,natom,dbl_mb(coor),byte_mb(tag))
         endif

         write(iw,1000) '@MD trajectory ',trajnum, ' MD step ', 
     &   step,' V: ',E1,' Ek: ',ek,' Tot E: ',Etot

          if (ga_nodeid().eq.0) then
             call crossing(natom,dbl_mb(coor),step)
          endif      
      enddo
C ---- End of MD_loop      
C      if (ga_nodeid().eq.0) then
C        call crossing(natom,dbl_mb(coor),step)
C      endif 
      call ga_sync()


      write(iw,*) 'END OF TRAJECTORY :',mc_data_get_trajnum()
      md_driver=.true. 
      return
 1000 format(A,i4,a,i4,a,F14.8,a,F14.8,a,F14.8)
 1001 format(i4)
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
C ********** NEW PREDICTOR
      subroutine predict1 (n,dt,r,v,a,b,c )
      implicit none
#include "global.fh"
        integer     n
        double precision r(3*n), v(3*n), a(3*n), b(3*n), c(3*n)
        double precision dt,sum0,sum1,sum2,sum3,sum4
        integer     i,j,iw
C ----  time -scaled variables, the expansion is Pascal triangle

c    *******************************************************************
        iw=6
       
c       write(iw,*) 'In predictor1,dt:',dt
       do i = 1,3*n
         r(i) = r(i) + dt*v(i) + (a(i)*(dt**2))/2.0
     1          +  (b(i)*(dt**3))/6.0 + (c(i)*(dt**4))/24.0
       enddo
       do i = 1,3*n
         v(i)= v(i)*dt + (dt**2)*a(i) +(b(i)*(dt**3))/2.0 
     1          + (c(i)*(dt**4))/6.0
       enddo
       do i = 1,3*n
         a(i) = a(i)*dt**2/2.0 + (dt**3)/2.0*b(i) + (dt**4)/4.0*c(i)
       enddo
       do i = 1,3*n
         b(i) = b(i)*dt**3/6.0 + (dt**4)/6.0*c(i)
       enddo
       do i = 1,3*n
         c(i) = (dt**4)/24.0*c(i)
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
c           write(6,*) 'corr i',corr(j)
           r(j,i) = r(j,i) + cr * corr(j)
           v(j,i) = v(j,i) + cv * corr(j)
           a(j,i) = ai(j)
           b(j,i) = b(j,i) + cb * corr(j)
           c(j,i) = c(j,i) + cc * corr(j)

           ek = ek + 0.5 * mass(i)*v(j,i) ** 2
         enddo
       enddo


       return
       end
C  ----- GEAR CORRECTOR NEW
      subroutine correct1 (n,dt,mass,r,v,a,b,c,force,ek )
      implicit none
#include "global.fh"
#include "const_data.fh"
        integer     n
        double precision r(3*n),v(3*n),a(3*n),b(3*n),c(3*n),force(3*n)
        double precision dt
        integer     i,j,iw
C ---- dt# are the Taylor expansion series terms dt1=dt, dt2=dt^2/2, dt3=dt^3/6 etc.
        double precision       mass(n), ek
        double precision ai(3*n),corr(3*n)
        double precision  c0, c1, c2, c3, c4
        parameter ( c0 = 19.0 / 120.0, c1 = 3.0 / 4.0, c2 = 1.0,
     1              c3 = 1.0 / 2.0,    c4 = 1.0 / 12.0 )

c    *******************************************************************
        iw=6
      
c        write(iw,*) 'In corrector 1, dt:',dt
        ek = 0.0

         do i=1,3*n
            if (abs(force(i)).lt.tolg) then
              force(i)=0.0
            endif
         enddo
c        do i=1,3*n
c           write(iw,*) 'In corrector 1, i:',i,r(i)
c        enddo
        do i = 1, n
         do j=1,3
           ai(j+(i-1)*3) = -force(j+(i-1)*3) / mass(i)
           corr(j+(i-1)*3) = (ai(j+(i-1)*3)*dt**2/2.0 - a(j+(i-1)*3))
c           write(6,*) 'corr i',corr(j+(i-1)*3)
           r(j+(i-1)*3) =  r(j+(i-1)*3) + c0 * corr(j+(i-1)*3)
           v(j+(i-1)*3) = (v(j+(i-1)*3) + c1 * corr(j+(i-1)*3)) * 1.0/dt
           a(j+(i-1)*3) = (a(j+(i-1)*3)+ c2*corr(j+(i-1)*3)) * 2.0/dt**2
           b(j+(i-1)*3) = (b(j+(i-1)*3)+ c3*corr(j+(i-1)*3)) * 6.0/dt**3
           c(j+(i-1)*3) = (c(j+(i-1)*3)+c4*corr(j+(i-1)*3)) * 24.0/dt**4
           ek = ek + 0.5 * mass(i)*v(j+(i-1)*3) ** 2
         enddo
       enddo


       return
       end
C ---- END OF GEAR CORRECTOR NEW       
       subroutine md_init_trj(n,dt,mass,r,v,a,b,c,grad,rtdb,theoryg)
       implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "const_data.fh"
C --- initialize higher derivatives for points dt/3 and 2dt/3
C --- r,v and a are non-zero, b and c are zero
C --- at the end of it, b and c will not be 0
       integer n,i,j,rtdb
       double precision r(3*n)
       double precision v(3*n)
       double precision a(3*n)
       double precision b(3*n)
       double precision c(3*n)
       double precision dt,dt1,dt2 
       double precision mass(n)
       integer tag,md_data_get_i_c,md_data_get_i_t
       integer coord
       integer h_r0,i_r0,h_r1,i_r1,h_r2,i_r2,h_f1,i_f1,h_f2,i_f2
       integer h_a1,i_a1,h_a2,i_a2
       logical grad 
       external grad
       character*64 theoryg
       integer iw
       logical dbug
       dbug=.false. 

       iw=6
       if (ga_nodeid().eq.0) then
         iw=6
       endif   

       if (dbug) then 
         write(iw,*) 'Starting with:'
         do i=1,3*n
           write(iw,9010) 'init values are:i,r,b,c:',i,r(i),b(i),c(i)
         enddo
         do i=1,n 
            write(iw,*) 'i,mass',i,mass(i)
         enddo
       endif
       if (.not. ma_alloc_get(mt_dbl, 3*n, 'r0' ,
     & h_r0,i_r0))
     & call errquit('md_data_getgeom: r0 not allocated',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'r1' ,
     & h_r1,i_r1))
     & call errquit('md_data_getgeom: r1 not allocated',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'r2' ,
     & h_r2,i_r2))
     & call errquit('md_data_getgeom: r2 not allocated',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'f1' ,
     & h_f1,i_f1))
     & call errquit('md_data_getgeom: f1',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'f2' ,
     & h_f2,i_f2))
     & call errquit('md_data_getgeom: f2',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'a1' ,
     & h_a1,i_a1))
     & call errquit('md_data_getgeom: a1',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'a2' ,
     & h_a2,i_a2))
     & call errquit('md_data_getgeom: a2',
     & 3*n,MA_ERR)

       tag=md_data_get_i_t()
       coord=md_data_get_i_c()
       dt1=dt/3.0
       dt2=2.0*dt/3.0
       write(iw,*) ' Initializing MD trajcory: points t/3 and 2t/3 '
C ---  take step dt/3 into r_one, and 2dt/3 into r_two
       if (dbug) then
        call print_coord(6,n,r,byte_mb(tag)) 
        call print_coord(6,n,v,byte_mb(tag))  
        call print_coord(6,n,a,byte_mb(tag))
        call print_coord(6,n,b,byte_mb(tag)) 
        call print_coord(6,n,c,byte_mb(tag))    
       endif
       do i=1,3*n
          dbl_mb(i_r0+i-1) = r(i) 
          dbl_mb(i_r1+i-1) = r(i)+v(i)*dt1+a(i)*dt1**2/2.0
     &                    + b(i)*dt1**3/6.0 + c(i)*dt1**4/24.0
          dbl_mb(i_r2+i-1) = r(i)+v(i)*dt2+a(i)*dt2**2/2.0
     &                    + b(i)*dt2**3/6.0 + c(i)*dt2**4/24.0
       enddo
       if(dbug) then 
        write(iw,*) 'Entering MD _init_traj: set coords'
        call print_coord(6,n,dbl_mb(i_r0),byte_mb(tag))
        write(iw,*) 'r1'
        call print_coord(6,n,dbl_mb(i_r1),byte_mb(tag))
        write(iw,*) 'r2'
        call print_coord(6,n,dbl_mb(i_r2),byte_mb(tag))
       endif
C  --- find the gradient for point dt/3
       do i=1,3*n
          r(i) = dbl_mb(i_r1+i-1)
       enddo
C ---  r(i) is really I-C data so, sould be set properlly
       call md_data_geom_set(rtdb)
C --- calculate forces for predicted coords             
       call ga_sync()
         if (.not. grad(rtdb)) 
     & call errquit('md_gradient: no grad calculated ', 0, RTDB_ERR)


         if (.not. rtdb_get(rtdb,theoryg,mt_dbl,3*n,dbl_mb(i_f1)))
     &   call errquit('md_gradient: no grad found ', 0, RTDB_ERR)
C  --- find the gradient for point 2*dt/3
         do i=1,3*n
            if (abs(dbl_mb(i_f1+i-1)).lt.tolg) then
              dbl_mb(i_f1+i-1)=0.0
            endif
         enddo
       do i=1,3*n
          r(i) = dbl_mb(i_r2+i-1)
       enddo
C ---  r(i) is really I-C data so, sould be set properlly
       call md_data_geom_set(rtdb)
       call ga_sync()
C --- calculate forces for predicted coords             
         if (.not. grad(rtdb))
     & call errquit('md_gradient: no grad calculated ', 0, RTDB_ERR)

         if (.not. rtdb_get(rtdb,theoryg,mt_dbl,3*n,dbl_mb(i_f2)))
     &   call errquit('md_gradient: no grad found ', 0, RTDB_ERR)
c       call print_coord(6,n,dbl_mb(i_f2),byte_mb(tag)) 
         do i=1,3*n
            if (abs(dbl_mb(i_f2+i-1)).lt.tolg) then
              dbl_mb(i_f2+i-1)=0.0
            endif
         enddo
c       write(iw,*) 'Entering MD _init_traj: idid two grad'
C --- need to transform force into acceleration
       do i=1,n
         do j=1,3
           dbl_mb(i_a1+(i-1)*3+j-1) = -dbl_mb(i_f1+(i-1)*3+j-1)/ mass(i)
           dbl_mb(i_a2+(i-1)*3+j-1) = -dbl_mb(i_f2+(i-1)*3+j-1)/ mass(i)
         enddo
       enddo
       if (dbug) then
        write(iw,*) 'f1/mass'
        call print_coord(6,n,dbl_mb(i_a1),byte_mb(tag))
        write(iw,*) 'f2/mass'
        call print_coord(6,n,dbl_mb(i_a2),byte_mb(tag))
       endif
C ---- set the equations and the initial coordinates
       do i=1,3*n
         b(i) = (4*dbl_mb(i_a1+i-1)-dbl_mb(i_a2+i-1)-3*a(i))/(dt1*2.0)
         c(i) = (dbl_mb(i_a2+i-1)-2*(dbl_mb(i_a1+i-1))+a(i))/(dt1*dt1)
         r(i) = dbl_mb(i_r0+i-1) 
c         write(iw,9010) 'init values are:i,r,b,c:',i,r(i),b(i),c(i)
       enddo
      write(iw,*) '**************************'
      write(iw,*) 'Done initializing MD driver'
      if(.not.ma_free_heap(h_a1))
     &     call errquit('memory deallocation:a1',0,RTDB_ERR)        
      if(.not.ma_free_heap(h_a2))
     &     call errquit('memory deallocation:a2',0,RTDB_ERR)  
      if(.not.ma_free_heap(h_r1))
     &     call errquit('memory deallocation:r1',0,RTDB_ERR)  
      if(.not.ma_free_heap(h_r2))
     &     call errquit('memory deallocation:r2',0,RTDB_ERR)  
      if(.not.ma_free_heap(h_f1))
     &     call errquit('memory deallocation:f1',0,RTDB_ERR)  
      if(.not.ma_free_heap(h_f2))
     &     call errquit('memory deallocation:f2',0,RTDB_ERR)  
      if(.not.ma_free_heap(h_r0))
     &     call errquit('memory deallocation:r0',0,RTDB_ERR)  
 9010 format (A30,I4,F11.6,F11.6,F11.6)
       return

       end

       subroutine print_coord(rfile,n,c,t)
       implicit none
       integer n,i,rfile
       double precision c(3,n),cnv
       character*16 t(n)
       cnv=1.889725989
c       if (rfile.eq.10) then
c        write(rfile,*) n 
c        write(rfile,*) 
c       endif  
       do i=1,n
          write(rfile,9000) t(i),c(1,i)/cnv,c(2,i)/cnv,c(3,i)/cnv
       enddo 
 9000 format (A4,F18.14,F18.14,F18.14) 
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
c      write(*,*) 'str1: ',str1
c      write(*,*) 'str2: ',str2
c      write(*,*) 'str: ',str
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
      i1=i 
      j=1
      do while ((j .le. len(str2)) .and. (str2(j:j) .ne. ' '))
            str (i1:i1) = str2 (j:j)
            i1 = i1 + 1
            j=j+1
      enddo

c      write (*,*) 'Long string : ', str
      write (*,*)
      end
      subroutine punch (rfile,n,nsteps,step,dt,mass,c,v,a,b,d,t)
      implicit none
        integer     rfile,n,nsteps,step 
        double precision mass(n),c(3,n), v(3,n), a(3,n), b(3,n), d(3,n)
        double precision dt
        character*16 t(n)  
         write(rfile,*) n,nsteps,step,dt,mass,c,v,a,b,d,t
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
        
c $Id$
