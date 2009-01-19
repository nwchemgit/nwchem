      function mc_main(rtdb)
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

      logical mc_main
      logical md_driver, leps_gradient,mc_driver,mc_init
      integer mcsteps
      external leps_gradient      
      logical status
      integer geom
      integer natom 
      integer c_start
      integer vel_init,acc_init
      logical forward,backward   ! tells you which part of the traj is being called
      logical mc_data_set_ifirc, mc_data_get_ifirc
      integer mc_data_get_forside, mc_data_get_forside0
      integer mc_data_get_backside, mc_data_get_backside0
      logical mc_data_set_backside0, mc_data_set_forside0
      logical mc_data_set_backside, mc_data_set_forside 
      integer i,j,k,l,naccept 
      integer c_array, vel, acc, prp 
      integer mc_data_get_i_c_array, md_data_get_i_v, md_data_get_i_a 
      integer md_data_get_i_c, mc_data_get_i_in_vel,mc_data_get_i_in_acc
      integer mc_data_get_natom, mc_data_get_i_prp
      integer tag,charge,md_data_get_i_t,md_data_get_i_q
      character*64 thr
      logical ircflag,flag,mc_data_set_forward,mc_data_set_backward
      integer i_sc, mc_data_get_i_sc, nonreact,reactive,reverse
      integer i_s, mc_data_get_i_s,mc_data_get_mcsteps
      double precision avg,rc,temp,mc_data_get_temp
      logical md_data_set_temp
      double precision mdtemp,mdtemp10, md_data_get_temp

      forward=.true.
      backward=.false.
      status=mc_data_set_ifirc(.true.)  
      status=mc_data_set_forward(forward)
      status=mc_data_set_backward(backward)
      status=mc_driver(rtdb,naccept)

      write(*,*) 'In MC_main, after mc_driver: naccept: ', naccept  
      ircflag=.true.
c      status=mc_init(rtdb,leps_gradient,ircflag)
c      write(*,*) 'In MC_main, after mc_init' 
      thr='leps'
      mdtemp10 =10.0 
      mdtemp = 298.15
      c_array=mc_data_get_i_c_array()
      vel_init=mc_data_get_i_in_vel()
      acc_init=mc_data_get_i_in_acc()
      natom=mc_data_get_natom()
c      do i=1,naccept
c        write(*,*) 'MC main accelerate: '
c        do k=1,natom*3
c           write(*,*) dbl_mb(acc_init+(i-1)*3*natom+k-1)
c        enddo
c      enddo

      call md_data_set_natom(natom)
      call md_data_allocate()
      call md_set(rtdb) ! sets md_temp, md_steps and mc_stepsize
      mdtemp=md_data_get_temp()
      tag=md_data_get_i_t()
      charge=md_data_get_i_q()
      c_start=md_data_get_i_c()
      vel=md_data_get_i_v()
      acc=md_data_get_i_a()
      i_sc=mc_data_get_i_sc()
      i_s=mc_data_get_i_s()
      prp=mc_data_get_i_prp()
      mcsteps=mc_data_get_mcsteps()

      write(*,*) 'In MC main coordinates'
      do i=1,naccept+1
         write(*,*) 'point ',i
         do k=1,3*natom
            write(*,*) dbl_mb(c_array+(i-1)*3*natom+k-1)
         enddo
      enddo  
      do i=1,naccept+1
        write(*,*) 'traj, prp: ',i,dbl_mb(prp+i-1) 
      enddo
      write(*,*) 'In MC main mdtemp, mdtemp10:', mdtemp, mdtemp10
      write(*,*) 'STARTING TRAJECTORIES '

C  ---- do the IRC trajectory first,  at low temp
      status=md_data_set_temp(mdtemp10)
C --- set vel and acc      
      status=mc_init(rtdb,leps_gradient,ircflag)
C  --- run the IRC traj forward 
      do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+k-1)
          dbl_mb(vel+k-1)=dbl_mb(vel_init+k-1)
          dbl_mb(acc+k-1)=dbl_mb(acc_init+k-1)
      enddo
      status=md_driver(rtdb,leps_gradient,thr,0)      
C --- run the IRC traj backward
      call md_set(rtdb)
      do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+k-1)
          dbl_mb(vel+k-1)=-dbl_mb(vel_init+k-1)
          dbl_mb(acc+k-1)=-dbl_mb(acc_init+k-1)
      enddo
      status=md_driver(rtdb,leps_gradient,thr,0) 
C --- get geometru from the big array in mc_driver
C--- the first one is always the TS
      write(*,*) 'Done with IRC'
      ircflag=.false.
      status=mc_init(rtdb,leps_gradient,ircflag)
      write(*,*) 'In MC_main, after mc_init' 
      flag =.true.
      do i=1,naccept+1
       if (i.gt.1) then
         status=mc_data_set_ifirc(.false.)
          status= mc_data_set_forside(0)
          status=mc_data_set_backside(0)
       else   
          status=mc_data_set_forside0(0)
          status=mc_data_set_backside0(0)
       endif 

       write(*,*) 'TRAJECTORY ',i
C ---- set the appropriate coordinate             
        do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+(i-1)*3*natom+k-1)
          dbl_mb(vel+k-1)=dbl_mb(vel_init+(i-1)*3*natom+k-1)
          dbl_mb(acc+k-1)=dbl_mb(acc_init+(i-1)*3*natom+k-1)
        enddo
        write(*,*) 'Coordinates:'
        do k=1,3*natom
         write(*,*) dbl_mb(c_start+k-1)
        enddo
        write(*,*) 'Velocieties:'
        do k=1,3*natom
         write(*,*) dbl_mb(vel+k-1)
        enddo
        write(*,*) 'Acceleration:'
        do k=1,3*natom
         write(*,*) dbl_mb(acc+k-1) 
        enddo


        do  j=1,3
          if (forward.and. .not.backward) then
           write(*,*) 'Forward: traj',i  
            status=md_driver(rtdb,leps_gradient,thr,i)
            forward=.false.
            backward=.true. 
            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)

          elseif(backward.and. .not.forward) then
            write(*,*) 'Backward: traj',i
C          reset the geometry to c_start
C          get the initial conditions generated in the forward run  
            call md_set(rtdb)
            do k=1,3*natom
              dbl_mb(c_start+k-1)=dbl_mb(c_array+(i-1)*3*natom+k-1)
              dbl_mb(vel+k-1)=-dbl_mb(vel_init+(i-1)*3*natom+k-1)
              dbl_mb(acc+k-1)=-dbl_mb(acc_init+(i-1)*3*natom+k-1)
            enddo
            status=md_driver(rtdb,leps_gradient,thr,i)
            forward=.false.
            backward=.false.
        write(*,*) 'traj ifirc forside0 backside0 forside backside ', i,
     &  mc_data_get_ifirc(),
     &  mc_data_get_forside0(), mc_data_get_backside0(),
     &  mc_data_get_forside(),mc_data_get_ backside()

            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)
          else
            write(*,*) 'Next traj '
            forward=.true.
            backward=.false.
            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)
            
c      frajectory finished, moved to the other trajectory  
C      here you would also check if the trajectory was reactive or not 
          endif 
        enddo 
c        write(*,*) 'traj ifirc forside0 backside0 forside backside ', i,
c     &  mc_data_get_ifirc(),
c     &  mc_data_get_forside0(), mc_data_get_backside0(), 
c     &  mc_data_get_forside(),mc_data_get_ backside()  
      enddo
      write(*,*) ' Recrossing Results:'
      nonreact=0
      reactive=0
      reverse=0
      do i=1,naccept+1
        write(*,*) 'traj  sidef sideb ',i, int_mb(i_sc+2*i-2), 
     &               int_mb(i_sc+2*i-1) 
        if (int_mb(i_sc+2*i-2).eq.int_mb(i_sc+2*i-1)) then
           nonreact=nonreact+1
           dbl_mb(prp+i-1)=0.0 
           write(*,*) 'non-reactive traj: ', i 
        else if (int_mb(i_sc+2*i-2).eq.mc_data_get_forside0()) then
           reactive=reactive+1
           dbl_mb(prp+i-1)=0.5*dbl_mb(prp+i-1)
        else
           reverse=reverse+1
           dbl_mb(prp+i-1)=-0.5*dbl_mb(prp+i-1)
        endif
      enddo

      j=0 
      avg=0 
      do i=0,mcsteps
        if (int_mb(i_s+i).eq.1) then
           j=j+1
        endif 
        write(*,*) 'adding prp :',j,' val: ',dbl_mb(prp+j-1)    
        avg=avg+dbl_mb(prp+j-1)
      enddo
      avg=avg/(mcsteps+1)
      temp = mc_data_get_temp()
      rc =avg/sqrt((temp*boltz)/(twopi))

      write(*,*) 'non-reactive, reactive,reverse traj: ', nonreact,
     &            reactive,reverse
      write(*,*)  ' Average prp value: ', avg
      write(*,*)  ' Recrossing coef: ',rc
      return    
      end








