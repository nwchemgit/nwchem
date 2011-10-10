       subroutine crossing(natom,c,mdstep)
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
      logical forward,backward,status
c  forside0 and backside0 are determined from the irc (at low T=10 K)
c  c is the coordintate of the current position and the ts
C  v_tst is the imaginary mode eigenvector 
      integer natom,trajnum,mdstep
      double precision dot,c(3*natom),v_tst(3*natom)
      double precision mass(3*natom)
      integer j,i,side,side_prev,nxing
      integer i_m,i_sc 
      integer i_ctst,i_evec,i_eval
      integer mc_data_get_i_ctst, mc_data_get_i_evec, mc_data_get_i_eval
      integer mc_data_get_i_m
      logical mc_data_set_forside, mc_data_set_backside
      logical mc_data_set_forside0, mc_data_set_backside0 
      integer mc_data_get_forside0, mc_data_get_backside0
      integer mc_data_get_forside, mc_data_get_backside
      integer mc_data_get_side_prev, mc_data_get_nxing
      integer mc_data_get_i_sc, mc_data_get_trajnum
      integer md_data_get_nsteps
      logical mc_data_get_ifirc,mc_data_get_forward,mc_data_get_backward
      logical mc_data_set_side_prev, mc_data_set_nxing
      double precision small
      integer iw, nsteps
      logical debug

      debug=.false.
      small =1E-8
      iw=6 
c      write(*,*) 'natom =', natom
      forward=mc_data_get_forward()
      backward=mc_data_get_backward()
      trajnum = mc_data_get_trajnum() 
c      write(*,*) 'In crossing: TST coordinates' 
c      write(*,*) ' Traj num: ',trajnum 
c      if (trajnum.eq.0) then
c         trajnum=1
c      endif
c      write(*,*) ' Traj num: ',trajnum 
      i_ctst=mc_data_get_i_ctst()
      i_evec=mc_data_get_i_evec()
      i_eval=mc_data_get_i_eval()
      i_m= mc_data_get_i_m()
      i_sc=mc_data_get_i_sc()
      nxing=mc_data_get_nxing()
      side_prev=mc_data_get_side_prev()
       nsteps=md_data_get_nsteps()
C      write(*,*) 'In crossing: i_sc: ', i_sc 
c      do i=1,3*natom
c         write(*,*) dbl_mb(i_ctst+i-1)
c      enddo
      do i=1,natom
        do j=1,3
          mass((i-1)*3+j)=dbl_mb(i_m+i-1)
       enddo
      enddo
C      goto 7
c      write(*,*) 'In crossing: C coordinates'
c      do i=1,3*natom
c         write(*,*) c(i)
c      enddo
c      write(*,*) 'In crossing: mass array'
c      do i=1,3*natom
c         write(*,*) mass(i)
c      enddo
C ...c is the current structure, c_tst is the structure of the TST
C ...v_tst is the imaginary mode vector, but in not mass-weighted coordinates
C .. the dot product of the Vtst-mw and C-mw is a measurement of the dividing syurface
C therefore we need to multiply them both by sqrt(mass), so we get C*M*V_TST
C ... you need to know on which side were you at the previous(side_prev)point, so that you can determine if 
C ... you have crossed accross the dividing surface at this point (side)
      do j=1,3*natom
c        write(*,*) 'In crossing: mode',dbl_mb(i_eval+j-1)
        if (dbl_mb(i_eval+j-1).lt.-1.0) then
          do i=1,3*natom
            v_tst(i)=dbl_mb(i_evec+(j-1)*3*natom+i-1)
          enddo
        endif
      enddo    
c      goto 7
c      write(*,*) 'In crossing: TST_V coordinates'
c      do i=1,3*natom
c         write(*,*) v_tst(i)
c      enddo

      dot =0.0
      do i=1,3*natom
         dot=dot+(c(i)-dbl_mb(i_ctst+i-1))*mass(i)*v_tst(i)
      enddo
c      write(*,*) 'Dot is: ',dot  
C
      if (abs(dot).lt.small) then
         side=0
      else
         if(dot.lt.0.0) then
            side=2
         else
            side=1
         endif
      endif
c      goto 7
c      write(*,*)'In crossing:trajnum ',trajnum 
c      write(iw,*)'In crossing: i_sc=',i_sc
C Do not cpunt trajectory as it is IRC
      if(.not.mc_data_get_ifirc()) then
       if (forward) then
c         call nga_put(i_sc,2*trajnum+1,1,side, ld)
c         int_mb(i_sc+(2*trajnum-2))=side
         int_mb(i_sc+(trajnum-1)*2*nsteps+mdstep-1) = side
c         write(iw,*)'In crossing: (i_sc+2*trajnum) ',trajnum ,
c     $  ' forw ',forward,' back ', backward,int_mb(i_sc+2*trajnum)
       else

c         int_mb(i_sc+(2*trajnum-1))=side
         int_mb(i_sc+(trajnum-1)*2*nsteps+nsteps+mdstep-1) = side
c      write(*,*)'In crossing: (i_sc+2*trajnum+1) ',trajnum ,
c     $  ' forw ',forward,' back ', backward,int_mb(i_sc+2*trajnum+1)
       endif
      endif

C
C     ----- Side ? Xing ? ... -----
C
      if ((side.ne.side_prev).and.(side_prev.ne.0)) then
         nxing=nxing+1
        write(iw,*) 'Recrossing: side_prev, side ',side_prev,side
        write(iw,*) 'Total nxing; ', nxing
      endif
c      goto 7
      status=mc_data_set_side_prev(side)
      status=mc_data_set_nxing(nxing)
  
      if (mc_data_get_ifirc()) then
        if (forward) then
           status= mc_data_set_forside0(side)
        endif
        if (backward) then 
           status= mc_data_set_backside0(side)
        endif
      else 
        if (forward) then
           status= mc_data_set_forside(side)
        endif
        if (backward) then
           status= mc_data_set_backside(side)
        endif
      endif
      
C forward backward are logical telling which part of the trajectory is running  
C their values are set in mc_main

C forside0 backside0 are the directions of the IRC trajectory
C they are set in this subroutine, when the ifirc flag is true
 
C forside backside  are the directions of the current trajectory
C are set in this subroutine 
      if (debug) then
      write(iw,9000) 'ifirc for back for0 back0 forside backside ',
     &   mc_data_get_ifirc(),  mc_data_get_forward(), 
     &   mc_data_get_backward(), mc_data_get_forside0(),
     &   mc_data_get_backside0(), mc_data_get_forside(), 
     &   mc_data_get_backside()
      endif
c      write(*,*) 'In crossing: nxing ',nxing
c    7 continue
 9000 format(A,3(L4),4(I4))
      return
      end
c $Id$
