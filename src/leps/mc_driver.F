      function mc_driver(rtdb)
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
      integer natom,mcstep !number of Monte Carlo step
      integer h_in,i_in
      integer h_c,i_c
      integer h_q,i_q
      integer h_t,i_t
      integer h_s,i_s
      double precision E0,E1,R,ediff,temp
      integer i,j,step  
      natom=3
      mcstep=5  
      temp=300
C will need an array of coordinaates      
C allocate memory for the init_points array
      write(*,*) 'Entering mc_driver'
      if (.not. ma_alloc_get(mt_dbl, (mcstep+1)*3*natom, 'coordinate' ,
     & h_in,i_in))
     & call errquit('mc_driver: array of coords not allocated',
     & 3*natom,MA_ERR)
      write(*,*) 'iAllocated big array'
C allocate  memory for the property array
      if (.not. ma_alloc_get(mt_int, mcstep+1, 'mcstep' ,
     & h_s,i_s))
     & call errquit('mc_driver: array of properties not allocated',
     & 1,MA_ERR)

C  allocate memeory for the next point to be stored
      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
     & h_c,i_c))
     & call errquit('mc_driver: array of coords not allocated',
     & 3*natom,MA_ERR)
      write(*,*) 'Allocated coord'
      if (.not. ma_alloc_get(mt_dbl, natom, 'charge', h_q,i_q))
     & call errquit('leps_getgeom: charge not allocated',
     & natom,MA_ERR)
      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
     & call errquit('leps_getgeom: tags not allocated',
     & natom,MA_ERR)


C copy the initial geometry in the array:
c     create geometry object
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('mc_driver_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
      write(*,*) 'Created Geom object'
C     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('mc_driver_getgeom: no geometry',
     & 0,RTDB_ERR)
      write(*,*) 'Loaded Geom object'

      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_in), dbl_mb(i_q)))
     & call errquit('mc_driver_getgeom: coord error',0,GEOM_ERR)
      write(*,*) 'Filled array'


C  ! get ihe energy for the initial point
      R=util_random(1234) ! init random umber generator
      do step=1,mcstep
      write(*,*) 'Starting step',step

      call leps_energy(rtdb)
      
      if (.not. rtdb_get(rtdb,'leps:energy',mt_dbl,1,E0))
     & call errquit('leps_energy:  no energy ', 0,
     & RTDB_ERR)

      write(*,*) 'In MC-driver E0 is: ',E0 
      write(*,*) 'In MC-driver stored in the big array:'
      do i=1, natom
         write(*,*) dbl_mb(i_in+(step-1)*3*natom+3*i-3),
     &              dbl_mb(i_in+(step-1)*3*natom+3*i-2),
     &              dbl_mb(i_in+(step-1)*3*natom+3*i-1)
      enddo
C   store the next step in the temp array
        do i=1, natom
         do j=1,3
          if (3*i-j.eq.2 .or. 3*i-j.eq.8) then 
C this part is to be replaced with the real displacment from the eigenvectors
           dbl_mb(i_c+3*i-j)=dbl_mb(i_in+(step-1)*3*natom+3*i-j)+0.01
          else
           dbl_mb(i_c+3*i-j)=dbl_mb(i_in+(step-1)*3*natom+3*i-j)  
          endif
          dbl_mb(i_in+step*3*natom+3*i-j)=dbl_mb(i_c+3*i-j)
         enddo
        enddo
      write(*,*) 'In MC-driver the next state is:'
      do i=1, natom
         write(*,*) dbl_mb(i_c+3*i-3),dbl_mb(i_c+3*i-2),
     &              dbl_mb(i_c+3*i-1)
      enddo

C store the new coords in rtdb
        if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &   dbl_mb(i_c), dbl_mb(i_q)))
     &   call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

        if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &   call errquit('mc_driver: Failed to store geom to rtdb',0, 
     &               RTDB_ERR)

C call the energy calculation 
        call leps_energy(rtdb)

        if (.not. rtdb_get(rtdb,'leps:energy',mt_dbl,1,E1))
     &   call errquit('mc_driver:  no energy ', 0,
     &                 RTDB_ERR)

        write(*,*) 'In MC-driver second energy is: ',E1
        if (E1.lt.E0) then
          int_mb(i_s+step)=1
        else
          R=util_random(0)
          ediff=exp(-(E1-E0)/(temp*boltzmann*JToHartree))
          if(ediff.le.R) then
             int_mb(i_s+step)=1
          else
             int_mb(i_s+step)=0
C     ---and return the old coordinates in the rtdb  
             do i=1,natom
               do j=1,3 
                dbl_mb(i_c+3*i-j)=dbl_mb(i_in+(step-1)*3*natom+3*i-j)
               enddo
             enddo 
             if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &            dbl_mb(i_c), dbl_mb(i_q)))
     &       call errquit('mc_driver_putgeom: coord error',0,GEOM_ERR)

             if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     &       call errquit('mc_driver: Failed to store geom to rtdb',0,
     &               RTDB_ERR)

         endif 
        endif
        
      enddo
      do i=1,mcstep
          write(*,*) 'property is ',int_mb(i_s+i) 
      enddo 

      mc_driver=.true.
      return 
      END
