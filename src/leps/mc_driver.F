      function mc_driver(rtdb,naccept)
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

      logical mc_driver
      integer geom
      integer mcstep,natom !number of Monte Carlo step
c      integer h_in,i_in
      integer i_c,i_evec,i_eval,i_in,i_s,i_q,i_t 
      integer mc_data_get_i_c, mc_data_get_i_evec, mc_data_get_i_eval
      integer mc_data_get_i_c_array, mc_data_get_i_s, mc_data_get_i_q
      integer mc_data_get_i_t
      integer mc_data_get_mcsteps
      double precision mc_data_get_stepsize, mc_data_get_ratio 
      double precision mc_data_get_temp
      double precision temp,stepsize,ratio
      logical  mc_data_set_natom, mc_data_set_naccept 
      double precision E0,E1,R,ediff,dir,dynrnd,urand
      integer i,j,step,ivib  
      integer naccept,nreject
      character*64 theory,blanks
      character*32 thr,task
      data blanks(1:64) /' '/
      double precision seed,nseed
      integer seed1
      COMMON/DYSEED/SEED,NSEED
 
      logical md_driver
      external md_driver
      logical leps_gradient
      external leps_gradient
      logical status 
      theory=blanks
      task=':energy'
c      thr="leps"
      call mc_data_inp_set(rtdb)
      natom=0
      seed1=152
      seed= 52.34
      mcstep=mc_data_get_mcsteps() ! number of monte carlo steps  
      temp=mc_data_get_temp()      ! temperature, should be the same as in md
      stepsize=mc_data_get_stepsize() ! stepsize scale for the geometry displacement
      ratio = mc_data_get_ratio()     ! the ratio of accepted and rejected steps in the run
                                      ! if the ratio is below or above the given value, the stepsize should be 
                                      ! increased or reduced NOT YET IMPLEMENTED 
      write(*,*) 'mymc input parameters:_mcsteps, temp, stepsize, ratio'
      write(*,*) mcstep,temp,stepsize,ratio
      naccept=0
      nreject=0
C will need an array of coordinaates      
C allocate memory for the init_points array
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, thr))
     $     call errquit('mc_driver: theory not specified',0, RTDB_ERR)
      call concat(thr,task,theory)
      write(*,*) 'Theory is:  ',theory

      if (.not. geom_create(geom, 'geometry'))
     & call errquit('mc_driver_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
      write(*,*) 'Created Geom object'
C     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('mc_driver_getgeom: no geometry',
     & 0,RTDB_ERR)
      write(*,*) 'Loaded Geom object'

      if (.not. geom_ncent(geom,natom))
     & call errquit('mc_data: ncent error',
     & 0,GEOM_ERR)

      status = mc_data_set_natom(natom)
      write(*,*) 'Entering mc_driver'
      write(*,*) 'Allocated big array'
 
      call mc_data_allocate(natom)
      call mc_data_set(rtdb,geom)
      write(*,*) 'Mc data set'

      i_c = mc_data_get_i_c()
      i_evec = mc_data_get_i_evec()
      i_eval =  mc_data_get_i_eval()
      i_in = mc_data_get_i_c_array() 
      i_s = mc_data_get_i_s()
      i_t = mc_data_get_i_t()
      i_q = mc_data_get_i_q() 
c      write(*,*) 'got handles'
c      do i=1, 3*natom
c          write(*,*) 'TST 0: ', dbl_mb(i_in+i-1)
c          write(*,*) 'TST 1: ', dbl_mb(i_c+i-1)
c      enddo

c      call mc_data_geom_copy(rtdb,geom)


      call mc_data_init(rtdb)
      write(*,*) 'mc_data_init'
      call printeg(3*natom,dbl_mb(i_eval),dbl_mb(i_evec))

C  ! get ihe energy for the initial point
C      R=util_random(seed1) ! init random number generator
      R=dynrnd(1) 
C      do i=1,10
C          write(*,*) 'dynrnd',dynrnd(0)
C      enddo 
C      return 

      int_mb(i_s)=1  ! always the first step is accepted cause it is the irc 
      do step=1,mcstep
c     write(*,*) 'Starting step',step

      call task_energy(rtdb) 
      
      if (.not. rtdb_get(rtdb,theory,mt_dbl,1,E0))
     & call errquit('mc_driver:  no energy ', 0,
     & RTDB_ERR)

      write(*,*)'In MC-driver:current step :',step,' from ',mcstep 
      do i=1, 3*natom
         write(*,*) dbl_mb(i_c+i-1)
      enddo
c      write(*,*)'In MC-driver: end of temp array C'
      do ivib=2, 3*natom
       if (dbl_mb(i_eval+ivib-1).gt.0.1) then  
c       write(*,*)'In MC-driver:mode',ivib,'is ',(dbl_mb(i_eval+ivib-1))
c       write(*,*)'In MC-driver:mode', ivib, ' vector is:'
c       do i=1, 3*natom
c        write(*,*)dbl_mb(i_evec+(ivib-1)*3*natom+i-1)
c       enddo
C            R=util_random(0)
C            dir=util_random(0)
            R=(dynrnd(0))
            dir=dynrnd(0)
            write(*,*) 'dynrnd1 ', R
            write(*,*) 'dynrnd2 ', dir  
            if(dir.le.0.5) then
               dir=-1.0
            else
               dir=1.0
            endif
            write(*,*)'In MC-driver:dir , rnd ', dir,R*dir*stepsize
            do i=1, 3*natom
              dbl_mb(i_c+i-1)=dbl_mb(i_c+i-1)+
     $              dbl_mb(i_evec+(ivib-1)*3*natom+i-1)*R*dir*stepsize
            enddo
       endif
      enddo

      write(*,*) 'In MC-driver the new state is:'
      do i=1, 3*natom
         write(*,*) dbl_mb(i_c+i-1)
      enddo

C store the new coords in rtdb
        if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &   dbl_mb(i_c), dbl_mb(i_q)))
     &   call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

        if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &   call errquit('mc_driver: Failed to store geom to rtdb',0, 
     &               RTDB_ERR)

C call the energy calculation 
C        call leps_energy(rtdb)
        call task_energy(rtdb)
        if (.not. rtdb_get(rtdb,theory,mt_dbl,1,E1))
     &   call errquit('mc_driver:  no energy ', 0,
     &                 RTDB_ERR)

       
        if (E1.lt.E0) then
          int_mb(i_s+step)=1
          naccept=naccept+1
          write(*,*) 'naccept1= ',naccept
          do i=1, 3*natom
             dbl_mb(i_in+naccept*3*natom+i-1)=dbl_mb(i_c+i-1)
          enddo
c          status=md_driver(rtdb,leps_gradient,thr,naccept)
c          call crossing(natom)
        else
C          R=util_random(0)
          R=(dynrnd(0))
          write(*,*) 'dynrnd3 ', R
          ediff=exp(-(E1-E0)/(temp*boltzmann*JToHartree))
          if(ediff.ge.R) then
             int_mb(i_s+step)=1
             naccept=naccept+1
             write(*,*) 'naccept2= ',naccept
             do i=1, 3*natom
                dbl_mb(i_in+naccept*3*natom+i-1)=dbl_mb(i_c+i-1)
             enddo
c             status=md_driver(rtdb,leps_gradient,thr,naccept)  
c             call crossing(natom)   
          else
             int_mb(i_s+step)=0
             nreject=nreject+1
             write(*,*) 'nreject= ',nreject
C     ---and return the old coordinates in the rtdb  
             do i=1,3*natom
                dbl_mb(i_c+i-1)=dbl_mb(i_in+naccept*3*natom+i-1)
             enddo 
             if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &            dbl_mb(i_c), dbl_mb(i_q)))
     &       call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

             if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &       call errquit('mc_driver: Failed to store geom to rtdb',0,
     &               RTDB_ERR)
             write(*,*) 'Returning to previous state :'
             do i=1,3*natom
                write(*,*) dbl_mb(i_c+i-1)
             enddo


         endif 
      write(*,*) 'In MC-driver E0,E1 :',E0,E1
      write(*,*) 'k*T :', (temp*boltzmann*JToHartree)
      write(*,*) 'Ediff :', ediff
      write(*,*) 'Rand :',R  
        endif
      enddo
      do i=0,mcstep
          write(*,*) 'property is ',int_mb(i_s+i) 
      enddo 
      write(*,*) 'number of accepted steps: ',naccept
      write(*,*) 'number of rejected steps: ', nreject
      status=mc_data_set_naccept(naccept)
      mc_driver=.true.
      return 
      END
      Subroutine printeg(nat3,eigval,eigvec)

       integer nat3, i,nhess
       double precision eigval(nat3)
       double precision eigvec(nat3*nat3)

       nhess = nat3*nat3  
       write(*,*) 'Printing Eigenvalues in mc_driver:',nat3
       do i=1, nat3
        write(*,*) eigval(i)
       enddo
       write(*,*) 'Printing Eigenvectors in mc_driver:',nhess
       do i=1, nat3
         write(*,*) 'Vector: ',i
         do j=1,nat3
          write(*,*) eigvec((i-1)*nat3+j)
         enddo
       enddo



      END

