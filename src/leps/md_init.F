      subroutine md_init (rtdb)
C      subroutine md_init (rtdb,n,mass,vel,acc,frc)
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
        integer     natom
c        double precision vel(3*n), acc(3,n),frc(3,n),mass(n) 
        integer seed
        double precision kB,vau,temp 
        double precision gks_gasdev, md_data_get_temp
        external gks_gasdev
        integer i,j,ivib
        integer mc_data_get_i_eval,mc_data_get_i_evec,i_eval,i_evec 
        integer md_data_get_i_m,md_data_get_i_f,mass,force
        integer md_data_get_i_v,md_data_get_i_a,vel,acc
        integer mc_data_get_natom

c        call md_data_inp_set(rtdb)

        seed = -1234  
        kB = boltzmann 
        temp = md_data_get_temp()   
        i_eval = mc_data_get_i_eval()
        i_evec = mc_data_get_i_evec()
        mass = md_data_get_i_m()
        force = md_data_get_i_f()
        vel = md_data_get_i_v()
        acc = md_data_get_i_a()
c        natom=n
        natom = mc_data_get_natom()
C--- get handles for some things we need, they were allocatet in md_driver
       write(*,*) 'In MD_init:'
      
       write(*,*) 'natom =', natom
       write(*,*) 'Temp=',temp

       call printeg(3*natom,dbl_mb(i_eval),dbl_mb(i_evec))
       write(*,*) 'Mass'
       do i=1,natom
         write(*,*) dbl_mb(mass+i-1)
       enddo
       write(*,*) 'Force'
       do i=1,3*natom
         dbl_mb(force+i-1)=0.0
c         vel(i)=0.0
       enddo

       do i=1,3*natom
         dbl_mb(vel+i-1)=0.0         
c         vel(i)=0.0
       enddo
       do ivib=1,3*natom
        if (abs(dbl_mb(i_eval+ivib-1)).lt.0.1) then
           vau=0.0
        else 
           vau=gks_gasdev(seed)*sqrt(kb*temp*1.2582847e+14)
        endif
        write(*,*) 'For ivib:',ivib,' vau is: ', vau
        do i=1,3*natom
           dbl_mb(vel+i-1) = dbl_mb(vel+i-1) +
     &     dbl_mb(i_evec+(ivib-1)*3*natom+i-1)*vau  
c              vel(i)=vel(i)
c     &        +dbl_mb(i_evec+(ivib-1)*3*natom+i-1)*vau
        enddo
       enddo
       do i=1,natom
           do j=1,3
              dbl_mb(acc+(i-1)*natom+j-1) = 
     &        -dbl_mb(force+(i-1)*natom+j-1)/dbl_mb(mass+i-1)
c              acc(i,j)=-frc(i,j)/mass(i)
           enddo
       enddo

      do i=1,natom
         do j=1,3
           write (*,*) 'mode ',(i-1)*natom+j, 
     &     ' vel ',dbl_mb(vel+(i-1)*natom+j-1)
         enddo 
      enddo

      return
      end 


