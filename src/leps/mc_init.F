      subroutine mc_init(rtdb,grad,ircflag)
      implicit none
      integer rtdb
#include "const_data.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "geom.fh"

      logical     ircflag 
      integer     natom,geom, i,j,k,l, ivib, seed
      integer     c_start,c_array
      integer     tag,charge,mass,force
      integer     mc_data_get_i_c_array,mc_data_get_i_c
      integer     mc_data_get_i_t, mc_data_get_i_q, mc_data_get_natom
      integer     mc_data_get_i_m, mc_data_get_i_f
      integer     naccept, mc_data_get_naccept
      integer     mc_data_get_i_in_vel, mc_data_get_i_in_acc
      integer     mc_data_get_i_prp 
      integer     prp,vel,acc  
      integer     i_eval,mc_data_get_i_eval,i_evec,mc_data_get_i_evec
      integer     md_data_get_i_v, md_data_get_i_a, mdvel,mdacc
      double precision kB,temp,temp10,vau,pimg      
      double precision mc_data_get_temp,gks_gasdev
      logical     grad
      external    grad       
      character*32 thr     ! this is just the theory ex: leps
      character*64 theory, blanks64  !theory should be leps:gradient 
      character*32 task, blanks32    ! this is ":gradient"
      data blanks64(1:64) /' '/
      data blanks32(1:32)/' '/
      double precision unitt,unitl,uma  
      logical dbug 
      theory=blanks64
      thr=blanks32
      task=blanks32
      task=':gradient'
      dbug=.true. 
      c_array=mc_data_get_i_c_array()  ! get the entire array of starting geometries
      c_start=mc_data_get_i_c()
      tag=mc_data_get_i_t()
      charge=mc_data_get_i_q()
      mass=mc_data_get_i_m()
      force=mc_data_get_i_f()   ! this is not set up yet 
      natom=mc_data_get_natom()
      naccept=mc_data_get_naccept()
      i_eval=mc_data_get_i_eval()
      i_evec=mc_data_get_i_evec()

      call mc_data_allocate_vel_acc(naccept)
      vel=mc_data_get_i_in_vel()
      acc=mc_data_get_i_in_acc()
      prp = mc_data_get_i_prp() 
C -- get the handles for MD to be used for IRC traj
      unitt   = 1.0D-15
      unitl   = 5.29177249D-09
C     FKC     = 627.509451D+00
      uma     = 1.6605655D-24
C                             ( = ATOMIC MASS UNIT IN GRAMS )
C      boltz   = 1.380662D-16
C                             ( = BOLTZMANN CONSTANT IN ERG/DEG. )


      seed = -1234
      kB = boltz
      temp = mc_data_get_temp()
      temp10=10.0   

      if (.not. geom_create(geom, 'geometry'))
     & call errquit('mc_driver_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
      write(*,*) 'MC_init created geom object'

      if (.not. rtdb_cget(rtdb, 'task:theory', 1, thr))
     $     call errquit('mc_driver: theory not specified',0, RTDB_ERR)
      call concat(thr,task,theory)
c      write(*,*) 'Theory is:  ',thr
C --- get geometry from the big array in mc_driver
C--- the first one is always the TS
C ---- set the appropriate coordinate
      if(ircflag) then
c         mdvel=md_data_get_i_v()
c         mdacc=md_data_get_i_a()
         temp = 10.0

         do k=1,3*natom
             dbl_mb(vel+k-1)=0.0
             dbl_mb(acc+k-1)=0.0
         enddo
         do k=1,3*natom
           dbl_mb(c_start+k-1)=dbl_mb(c_array+k-1)
         enddo
         if (.not. geom_cart_set(geom, natom, byte_mb(tag),
     &       dbl_mb(c_start), dbl_mb(charge)))
     &       call errquit('mc_init_putgeom: coord error',0,GEOM_ERR)

         if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &       call errquit('mc_init: Failed to store geom to rtdb',0,
     &       RTDB_ERR)

C---- calculate the gradient to get the force for acceleration
         if (.not.grad(rtdb))
     &       call errquit('mc_init: no grad calculated ', 0, RTDB_ERR)

         if (.not. rtdb_get(rtdb,theory,mt_dbl,3*natom, dbl_mb(force)))
     &       call errquit('mc_grad: no grad stored ', 0, RTDB_ERR)

         write(*,*) "MC-init IRC: temperature is:",temp
         do ivib=1,3*natom
            if (abs(dbl_mb(i_eval+ivib-1)).lt.0.1) then
                vau=0.0
            else
                vau=gks_gasdev(seed)*sqrt(kb*temp)
                vau=vau*(unitt/unitl/sqrt(uma))
            endif
            write(*,*) 'For ivib IRC:',ivib,' vau is: ', vau
            do k=1,3*natom
               dbl_mb(vel+k-1) =
     &            dbl_mb(vel+k-1)
     &          + dbl_mb(i_evec+(ivib-1)*3*natom+k-1)*vau
            enddo
         enddo
         write(*,*) 'IRC Velocities:'
         do k=1,3*natom
            write(*,*) dbl_mb(vel+k-1)
         enddo

c----- set the  acceleration
         do k=1,natom
            do l=1,3
            dbl_mb(acc+(k-1)*natom+l-1) =
     &          -dbl_mb(force+(k-1)*natom+l-1)/dbl_mb(mass+k-1)
            enddo
         enddo

         write(*,*) 'IRC Acceleration:'
         do k=1,3*natom
            write(*,*) dbl_mb(acc+k-1)
         enddo
      else 
      write(*,*) 'in mc_init:naccept',naccept
      do i=1,naccept+1   ! loop over all of the initial structures
         pimg=0.0 
         do k=1,3*natom
           dbl_mb(c_start+k-1)=dbl_mb(c_array+(i-1)*3*natom+k-1)
         enddo
         if (.not. geom_cart_set(geom, natom, byte_mb(tag),
     &       dbl_mb(c_start), dbl_mb(charge)))
     &       call errquit('mc_init_putgeom: coord error',0,GEOM_ERR)

         if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &       call errquit('mc_init: Failed to store geom to rtdb',0,
     &       RTDB_ERR)

C---- calculate the gradient to get the force for acceleration
         if (.not.grad(rtdb))
     &       call errquit('mc_init: no grad calculated ', 0, RTDB_ERR)

         if (.not. rtdb_get(rtdb,theory,mt_dbl,3*natom, dbl_mb(force)))
     &       call errquit('mc_grad: no grad stored ', 0, RTDB_ERR)

C----- set the velocity and acc to zero
         do k=1,3*natom
             dbl_mb(vel+(i-1)*3*natom+k-1)=0.0
             dbl_mb(acc+(i-1)*3*natom+k-1)=0.0 
         enddo
         write(*,*) "MC-init: temperature is:",temp  
         do ivib=1,3*natom
            if (abs(dbl_mb(i_eval+ivib-1)).lt.0.1) then
                vau=0.0
            else
                vau=gks_gasdev(seed)*sqrt(kb*temp)
                if (ivib.eq.1) then
                    pimg=vau
                    dbl_mb(prp+i-1)=pimg
                endif
                vau=vau*(unitt/unitl/sqrt(uma))
            endif
            write(*,*) 'For ivib:',ivib,' vau is: ', vau
            do k=1,3*natom
               dbl_mb(vel+(i-1)*3*natom+k-1) = 
     &            dbl_mb(vel+(i-1)*3*natom+k-1)
     &          + dbl_mb(i_evec+(ivib-1)*3*natom+k-1)*vau
            enddo
         enddo
         write(*,*) 'Velocities:'
         do k=1,3*natom
            write(*,*) dbl_mb(vel+(i-1)*3*natom+k-1)
         enddo 
         write(*,*) 'PIMG= ',pimg
c----- set the  acceleration
         do k=1,natom
            do l=1,3
            dbl_mb(acc+(i-1)*3*natom+(k-1)*natom+l-1) = 
     &          -dbl_mb(force+(k-1)*natom+l-1)/dbl_mb(mass+k-1)
            enddo
         enddo 

         write(*,*) 'Acceleration:'
         do k=1,3*natom
            write(*,*) dbl_mb(acc+(i-1)*3*natom+k-1) 
         enddo
C ---- end of setting
      enddo             
      endif
      if(dbug) then
        do i=1,naccept
         write(*,*) 'Dbug Velocities:'
         do k=1,3*natom
            write(*,*) dbl_mb(vel+(i-1)*3*natom+k-1)
         enddo
         write(*,*) 'Dbug Acceleration:'
         do k=1,3*natom
            write(*,*) dbl_mb(acc+(i-1)*3*natom+k-1)
         enddo
 
        enddo
      endif         
      if (.not. geom_destroy(geom))
     &   call errquit('mc_init:  geom_destroy failed',
     &   911,GEOM_ERR)

      end
C_CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
