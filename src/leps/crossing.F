       subroutine crossing(natom,c,trajnum,side_prev)
c      subroutine crossing(c,c_tst,v_tst,side,side_prev,nxing,forward,backward,forside,backside)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
      logical forward,backward,status,ifirc
c  forside0 and backside0 are determined from the irc (at low T=10 K)
c  c is the coordintate of the current position and the ts
C  v_tst is the imaginary mode eigenvector 
      integer natom,trajnum
      double precision dot,c(3*natom),c_tst(3*natom),v_tst(3*natom)
      double precision mass(3*natom)
      integer j,i,nxing,side,side_prev
      integer forside,backside
      integer i_c,h_c,i_m,h_m,i_sc 
      integer i_ctst,i_evec,i_eval
      integer mc_data_get_i_ctst, mc_data_get_i_evec, mc_data_get_i_eval
      integer mc_data_get_i_m
      logical mc_data_set_forside, mc_data_set_backside
      logical mc_data_set_forside0, mc_data_set_backside0 
      integer mc_data_get_forside0, mc_data_get_backside0
      integer mc_data_get_forside, mc_data_get_backside
      integer mc_data_get_i_sc
      logical mc_data_get_ifirc,mc_data_get_forward,mc_data_get_backward
      double precision small
      small =1E-8
      write(*,*) 'natom =', natom
      forward=mc_data_get_forward()
      backward=mc_data_get_backward()
 
      write(*,*) 'In crossing: TST coordinates' 
      write(*,*) ' Traj num: ',trajnum 
c      if (trajnum.eq.0) then
c         trajnum=1
c      endif
c      write(*,*) ' Traj num: ',trajnum 
      i_ctst=mc_data_get_i_ctst()
      i_evec=mc_data_get_i_evec()
      i_eval=mc_data_get_i_eval()
      i_m= mc_data_get_i_m()
      i_sc=mc_data_get_i_sc()
      write(*,*) 'In crossing: i_sc: ', i_sc 
      do i=1,3*natom
         write(*,*) dbl_mb(i_ctst+i-1)
      enddo
      do i=1,natom
        do j=1,3
          mass((i-1)*natom+j)=dbl_mb(i_m+i-1)
       enddo
      enddo
C      goto 7
      write(*,*) 'In crossing: C coordinates'
      do i=1,3*natom
         write(*,*) c(i)
      enddo
      write(*,*) 'In crossing: mass array'
      do i=1,3*natom
         write(*,*) mass(i)
      enddo
C ...c is the current structure, c_tst is the structure of the TST
C ...v_tst is the imaginary mode vector, but in not mass-weighted coordinates
C .. the dot product of the Vtst-mw and C-mw is a measurement of the dividing syurface
C therefore we need to multiply them both by sqrt(mass), so we get C*M*V_TST
C ... you need to know on which side were you at the previous(side_prev)point, so that you can determine if 
C ... you have crossed accross the dividing surface at this point (side)
      do j=1,3*natom
        write(*,*) 'In crossing: mode',dbl_mb(i_eval+j-1)
        if (dbl_mb(i_eval+j-1).lt.-1.0) then
          do i=1,3*natom
            v_tst(i)=dbl_mb(i_evec+(j-1)*3*natom+i-1)
          enddo
        endif
      enddo    
c      goto 7
      write(*,*) 'In crossing: TST_V coordinates'
      do i=1,3*natom
         write(*,*) v_tst(i)
      enddo

      dot =0.0
      do i=1,3*natom
         dot=dot+(c(i)-dbl_mb(i_ctst+i-1))*mass(i)*v_tst(i)
      enddo
      write(*,*) 'Dot is: ',dot  
C
      if (abs(dot).lt.small) then
         side=0
      else
         if(dot.lt.0.0) then
            side=-1
         else
            side=1
         endif
      endif
c      goto 7
      if (forward) then
         int_mb(i_sc+2*trajnum+1)=side
      write(*,*)'In crossing: (i_sc+2*trajnum+1) ',trajnum ,
     $  ' forw ',forward,' back ', backward,int_mb(i_sc+2*trajnum+1)
      else
         int_mb(i_sc+2*trajnum+2)=side
      write(*,*)'In crossing: (i_sc+2*trajnum+2) ',trajnum ,
     $  ' forw ',forward,' back ', backward,int_mb(i_sc+2*trajnum+2)
      endif

C
C     ----- Side ? Xing ? ... -----
C
c      if ((side.ne.side_prev).and.(side_prev.ne.0)) then
c         nxing=nxing+1
c      endif
      write(*,*) 'In crossing: side_prev ',side_prev
c      goto 7
      side_prev=side
      
      if (mc_data_get_ifirc()) then
        if (forward) then
           write(*,*) 'setting forside0'
           status= mc_data_set_forside0(side)
        endif
        if (backward) then 
           write(*,*) 'setting backside0'
           status= mc_data_set_backside0(side)
        endif
      else 
        if (forward) then
           write(*,*) 'setting forside'
           status= mc_data_set_forside(side)
        endif
        if (backward) then
           write(*,*) 'setting backside'
           status= mc_data_set_backside(side)
        endif
      endif
      
C forward backward are logical telling which part of the trajectory is running  
C their values are set in mc_main

C forside0 backside0 are the directions of the IRC trajectory
C they are set in this subroutine, when the ifirc flag is true
 
C forside backside  are the directions of the current trajectory
C are set in this subroutine 
      write(*,*) 'ifirc forward backward forside0 backside0 forside back
     &side ', 
     &   mc_data_get_ifirc(),  mc_data_get_forward(), 
     &   mc_data_get_backward(), mc_data_get_forside0(),
     &   mc_data_get_backside0(), mc_data_get_forside(), 
     &   mc_data_get_backside()
c      write(*,*) 'In crossing: nxing ',nxing
c    7 continue

      end
