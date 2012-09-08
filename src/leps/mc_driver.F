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
#include "util_sgroup.fh"  !for subgroup support
c#include "subgr.fh"

      logical mc_driver, task_energy
      external task_energy
      integer geom
      integer mcstep,natom !number of Monte Carlo step
      integer i_c,i_evec,i_eval,i_in,i_s,i_q,i_t 
      integer mc_data_get_i_c, mc_data_get_i_evec, mc_data_get_i_eval
      integer mc_data_get_i_c_array, mc_data_get_i_s, mc_data_get_i_q
      integer mc_data_get_i_t
      integer mc_data_get_mcsteps
      double precision mc_data_get_stepsize, mc_data_get_ratio 
      double precision mc_data_get_temp
      double precision temp,stepsize,ratio
      logical  mc_data_set_natom, mc_data_set_naccept 
      double precision E0,E1,R,ediff,dir,dynrnd
      integer i,step,ivib  
      integer naccept,nreject
      character*64 theory,blanks
      character*32 thr,task
      data blanks(1:64) /' '/
      double precision seed,nseed,mc_data_get_seed
      integer iworld,idef
      COMMON/DYSEED/SEED,NSEED
 
      logical md_driver
      external md_driver
      logical status, task_qmmm 
      integer iw
      integer cntnrg
      theory=blanks
      task=':energy'
      call mc_data_inp_set(rtdb)
      natom=0
      seed= 523.4
      seed = mc_data_get_seed()
      mcstep=mc_data_get_mcsteps() ! number of monte carlo steps  
      temp=mc_data_get_temp()      ! temperature, should be the same as in md
      stepsize=mc_data_get_stepsize() ! stepsize scale for the geometry displacement
      ratio = mc_data_get_ratio()     ! the ratio of accepted and rejected steps in the run
                                      ! if the ratio is below or above the given value, the stepsize should be 
                                      ! increased or reduced NOT YET IMPLEMENTED 
      cntnrg=0
      iw=6
      mc_driver=.false.      
      task_qmmm=.false.                
      write(iw,*) 'MC driver input parameters:\n
     $  number of steps, temperature, stepsize, ratio',
     &  mcstep,temp,stepsize,ratio
      naccept=0
      nreject=0
C will need an array of coordinaates      
C allocate memory for the init_points array
      if( .not. rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &    task_qmmm = .false.

      if (.not. rtdb_cget(rtdb, 'task:theory', 1, thr))
     $     call errquit('mc_driver: theory not specified',0, RTDB_ERR)
      call concat(thr,task,theory)
C      write(*,*) 'Theory is:  ',theory

      if(task_qmmm) then
        theory='qmmm:energy'
      endif
   
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('mc_driver_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
C      write(*,*) 'Created Geom object'
C     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('mc_driver_getgeom: no geometry',
     & 0,RTDB_ERR)
C      write(*,*) 'Loaded Geom object'

      if (.not. geom_ncent(geom,natom))
     & call errquit('mc_data: ncent error',
     & 0,GEOM_ERR)


      status = mc_data_set_natom(natom)
C      write(*,*) 'Entering mc_driver'
C      write(*,*) 'Allocated big array'
 
      call mc_data_allocate(natom,rtdb)
      call mc_data_set(rtdb,geom)
C      write(*,*) 'Mc data set'

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
c      write(iw,*) 'ieigen values and vectors'
c      call printeg(3*natom,dbl_mb(i_eval),dbl_mb(i_evec))
C    find out if it is a TS
      if (dbl_mb(i_eval).gt.-0.1.or.abs(dbl_mb(i_eval+1)).gt.0.1)  then
        write (iw,*)  'The initial structure is not a transition state:'
        write (iw,*) dbl_mb(i_eval)
        write (iw,*) dbl_mb(i_eval+1)
        return
      endif      

C*********** This is where the restart would come in **********
C *********** skip the rest of the subroutine *****************   
C**************************************************************   
C  ! get ihe energy for the initial point
      R=dynrnd(1) 

      int_mb(i_s)=1  ! always the first step is accepted cause it is the irc 

C ---- get the energy of the TST (0th step)
      status = task_energy(rtdb)
      cntnrg=cntnrg+1
c      write (iw,*) 'cntnrg = ',cntnrg 
      if (.not. rtdb_get(rtdb,theory,mt_dbl,1,E0))
     & call errquit('mc_driver:  no energy ', 0,
     & RTDB_ERR)

      do step=1,mcstep
c      write(*,*) 'Starting step',step
      write(iw,*)'In MC-driver:current step :',step,' from ',mcstep
c      write(*,*)'In MC-driver:current geometry '
C      call printcoord(natom,dbl_mb(i_c),byte_mb(i_t))

C --- if the previous value is 0, it means the previous step was rejacted
C --- so E0 is E0, but if it is 1, the previous step was accepted, so E0=E1 
c      if (int_mb(i_s+step-1).ne.0) then
c     &   status = task_energy(rtdb) 
      
c      if (.not. rtdb_get(rtdb,theory,mt_dbl,1,E0))
c     & call errquit('mc_driver:  no energy ', 0,
c     & RTDB_ERR)
c      endif
c      call printcoord(natom,dbl_mb(i_c),byte_mb(i_t))

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
c            write(*,*) 'dynrnd1 ', R
c            write(*,*) 'dynrnd2 ', dir  
            if(dir.le.0.5) then
               dir=-1.0
            else
               dir=1.0
            endif
c            write(*,*)'In MC-driver:dir , rnd ', dir,R*dir*stepsize
            do i=1, 3*natom
              dbl_mb(i_c+i-1)=dbl_mb(i_c+i-1)+
     $              dbl_mb(i_evec+(ivib-1)*3*natom+i-1)*R*dir*stepsize
            enddo
       endif
      enddo
C store the new coords in rtdb
        if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &   dbl_mb(i_c), dbl_mb(i_q)))
     &   call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

        if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &    call errquit('mc_driver: Failed to store geom to rtdb',0, 
     &               RTDB_ERR)
        call ga_sync()
C call the energy calculation 
c        write(*,*)'In MC-driver: new geometry :'
        status = task_energy(rtdb)
      cntnrg=cntnrg+1
c      write (iw,*) 'cntnrg = ',cntnrg
        if (.not. rtdb_get(rtdb,theory,mt_dbl,1,E1))
     &   call errquit('mc_driver:  no energy ', 0,
     &                 RTDB_ERR)

        write(iw,*) ' In MC-driver E0,E1 :',E0,E1       
        if (E1.lt.E0) then
          int_mb(i_s+step)=1
          naccept=naccept+1
c          write(*,*) 'naccept1= ',naccept
          do i=1, 3*natom
             dbl_mb(i_in+naccept*3*natom+i-1)=dbl_mb(i_c+i-1)
          enddo
          E0=E1
          write(iw,*) 'In MC-driver the new state is:'
          call printcoord(iw,natom,dbl_mb(i_c),byte_mb(i_t))
        else
C          R=util_random(0)
          R=(dynrnd(0))
c          write(*,*) 'dynrnd3 ', R
          ediff=exp(-(E1-E0)/(temp*boltzmann*JToHartree))
          if(ediff.ge.R) then
             int_mb(i_s+step)=1
             naccept=naccept+1
c             write(*,*) 'naccept2= ',naccept
             do i=1, 3*natom
                dbl_mb(i_in+naccept*3*natom+i-1)=dbl_mb(i_c+i-1)
             enddo
             E0=E1
             write(iw,*) 'In MC-driver the new state is:'
             call printcoord(iw,natom,dbl_mb(i_c),byte_mb(i_t))
          else
             int_mb(i_s+step)=0
             nreject=nreject+1
             E0=E0
c             write(*,*) 'nreject= ',nreject
C     ---and return the old coordinates in the rtdb  
             do i=1,3*natom
                dbl_mb(i_c+i-1)=dbl_mb(i_in+naccept*3*natom+i-1)
             enddo 
             if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &            dbl_mb(i_c), dbl_mb(i_q)))
     &        call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

             if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &         call errquit('mc_driver: Failed to store geom to rtdb',0,
     &               RTDB_ERR)
             call ga_sync()
              write(iw,*) 'Step rejected. Returning to previous state :'
C ----       no need to recalculate energy at previous step
C ----       just make a step  

c             do i=1,3*natom
c                write(*,*) dbl_mb(i_c+i-1)
c             enddo


         endif 
C      write(*,*) 'In MC-driver the new state is:'
C      call printcoord(natom,dbl_mb(i_c),byte_mb(i_t))
c      write(*,*) 'k*T :', (temp*boltzmann*JToHartree)
c      write(*,*) 'Ediff :', ediff
c      write(*,*) 'Rand :',R  
        endif
      write(iw,*) 'In MC-driver E0,E1 -two:',E0,E1
      enddo
       call ga_sync()
       do i=0,mcstep
          write(iw,*) 'property is ',int_mb(i_s+i) 
       enddo 
       write(iw,*) 'number of accepted steps: ',naccept
       write(iw,*) 'number of rejected steps: ', nreject
      status=mc_data_set_naccept(naccept)
c       if(ga_nodeid().eq.0.and.util_sgroup_mygroup().eq.1) then
c         call print_div(dbl_mb(i_in),natom,naccept,byte_mb(i_t))
c       endif
       if (.not. geom_destroy(geom, 'geometry'))
     & call errquit('mc_driver: failed destroing geometry',
     & 0,GEOM_ERR)
      mc_driver=.true.
      return 
      END

      subroutine print_div(coords,nat,naccept,tags)
        double precision coord(3*nat)
        character*16 tags(nat)
        integer i,j,nat,nfile

        nfile=13
        open (nfile,file='dividing_surface.xyz', status='REPLACE') 
        do i=0,naccept-1
           do j=1,3*nat
              coord(j)=coords(i*3*nat+j) 
           enddo 
           call printcoord(nfile,nat,coord,tags)
        enddo
        close(13)

      end
      subroutine printcoord(nfile,nat,coord,tags)
#include "global.fh"
         integer i,nat,iw,nfile
         double precision coord(3,nat)
         character*16 tags(nat) 
         nfile=6 
         iw=nfile
         do i=1,nat
            write(iw,9000) tags(i),coord(1,i),coord(2,i),coord(3,i)
         enddo 
 9000 format (A4,F14.9,F14.9,F14.9) 
      END 

      Subroutine printeg(nat3,eigval,eigvec)
#include "global.fh"
 
       integer nat3, i,nhess
       double precision eigval(nat3)
       double precision eigvec(nat3*nat3)
       integer iw
       iw=6
       
       nhess = nat3*nat3  
       write(*,*) 'Printing Eigenvalues in mc_driver:',nat3
       do i=1, nat3
        write(iw,*) eigval(i)
       enddo
       write(iw,*) 'Printing Eigenvectors in mc_driver:',nhess
       do i=1, nat3
         write(iw,*) 'Vector: ',i
         do j=1,nat3
          write(iw,*) eigvec((i-1)*nat3+j)
         enddo
       enddo



      END

c $Id$
