      subroutine tce_print_x1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c $Id$
c
c Prints the X1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
         write(LuOut,9000) 'Singles contributions'
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1          .eq.irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
                if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1            call errquit('tce_print_x1: MA problem',0,MA_ERR)
                call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1            int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h2 = 1,int_mb(k_range+h2b-1)
                    i = i + 1
                    if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,
     1                                   irrep_p)
                      call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,
     1                                   irrep_h)
                      spin_p = spinnames(int_mb(k_spin+p1b-1))
                      spin_h = spinnames(int_mb(k_spin+h2b-1))
                      if (nodezero) write(LuOut,9010) 
c     1                int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
c     2                int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
c     3                dbl_mb(k_r1+i-1)
     1                (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1                 2-int_mb(k_spin+p1b-1))/2,
     1                irrep_p,spin_p,
     2                (int_mb(k_mo_index+int_mb(k_offset+h2b-1)+h2-1)+
     1                 2-int_mb(k_spin+h2b-1))/2,
     2                irrep_h,spin_h,
     2                dbl_mb(k_r1+i-1)
                    endif
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r1))
     1            call errquit('tce_print_x1: MA problem',1,MA_ERR)
              endif
            endif
          endif
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end

      subroutine tce_print_y1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c Prints the Y1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
         write(LuOut,9000) 'Singles contributions'
      do h2b = 1,noab
        do p1b = noab+1,noab+nvab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_y1: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((h2b-1)*nvab+p1b-noab-1))
            i = 0
            do p1 = 1,int_mb(k_range+p1b-1)
              do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
                if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                  call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p)
                  call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,irrep_h)
                  spin_p = spinnames(int_mb(k_spin+p1b-1))
                  spin_h = spinnames(int_mb(k_spin+h2b-1))
                  if (nodezero) write(LuOut,9010) 
     1              int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
     2              int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
     3              dbl_mb(k_r1+i-1)
                endif
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_y1: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
c
c
      subroutine tce_print_g1(rtdb,d_r1,k_r1_offset,r1_thresh,irrep_r1)
c
c Prints the G1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb
      integer d_r1
      integer g1b
      integer g2b
      integer g1
      integer g2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep_r1
      integer eiformat
      logical nodezero
      double precision r1_thresh
      character*4 irrep_g1,irrep_g2
      character*5 spinnames(2),spin_g1,spin_g2
      data spinnames/'alpha','beta '/
c
      if (.not.rtdb_get(rtdb,'tce:eiformat',mt_int,1,eiformat)) then
        eiformat = 1
        if (.not.rtdb_put(rtdb,'tce:eiformat',mt_int,1,eiformat))
     1  call errquit('tce_init: failed to write eiformat',0,RTDB_ERR)
      endif
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'Singles contributions'
      if(nodezero) then
      do g2b = 1,noab+nvab
        do g1b = 1,noab+nvab
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq.irrep_r1) then
            size = int_mb(k_range+g1b-1) * int_mb(k_range+g2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_y1: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),(g2b-1)*(noab+nvab)+g1b-1)
            i = 0
            do g1 = 1,int_mb(k_range+g1b-1)
             do g2 = 1,int_mb(k_range+g2b-1)
               i = i + 1
               if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                 call sym_irrepname(geom,int_mb(k_sym+g1b-1)+1,irrep_g1)
                 call sym_irrepname(geom,int_mb(k_sym+g2b-1)+1,irrep_g2)
                 spin_g1 = spinnames(int_mb(k_spin+g1b-1))
                 spin_g2 = spinnames(int_mb(k_spin+g2b-1))
                 if ((nodezero).and.(eiformat.eq.1)) then
                   write(LuOut,9011)
     1             int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2             int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     3             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.2)) then
                   write(LuOut,9012)
     1             int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2             int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     3             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.3)) then
                   write(LuOut,9013)
     1             int_mb(k_offset+g1b-1)+g1,int_mb(k_sym+g1b-1)+1,
     2             int_mb(k_spin+g1b-1),
     3             int_mb(k_offset+g2b-1)+g2,int_mb(k_sym+g2b-1)+1,
     4             int_mb(k_spin+g2b-1),
     5             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.4)) then
                   write(LuOut,9014)
     1             int_mb(k_offset+g1b-1)+g1,
     2             int_mb(k_offset+g2b-1)+g2,
     3             dbl_mb(k_r1+i-1)          
                 endif
               endif
             enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_y1: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9011 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.14)
 9012 format(1x,i3,',  ',a4,',',a5,',',i3,',  ',a4,',',a5,',',f20.14)
 9013 format(1x,i3,',',i3,',',i3,',',i3,',',i3,',',i3,',',f20.14)
 9014 format(1x,i3,', ',i3,', ',f20.14)
      return
      end
c
      subroutine tce_string_x1(d_r1,k_r1_offset,r1_thresh,irrep,xnorm)
c
c $Id: tce_print_x1.F 26707 2015-01-17 01:06:43Z edo $
c
c Prints the X1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep
c string
       integer spinsum,j,iposh,iposp
       integer nosa,nosb,nvsa,nvsb
c ------
      logical nodezero
      double precision r1_thresh,xnorm
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
c string 
      character*6 refstring(202)
      character*1 bra,num1
      character*2 num2
      character*3 keta,ketb,kea,keb
      character*4 cre4
      character*5 cre5
      character*6 cre6
c ------
      data spinnames/'alpha','beta '/
c
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
cc        write(6,*)'check 2: xnorm =',xnorm
cc        call util_flush(6)
c
c string
c
        bra='('
        keta='a)+'
        ketb='b)+'
        kea='a)'
        keb='b)'
c
        nosa=0
        nosb=0
        nvsa=0
        nvsb=0
        do i=1,noa
         nosa=nosa+int_mb(k_range+i-1)
        enddo
        do i=1,nob
         nosb=nosb+int_mb(k_range+noa+i-1)
        enddo
        do i=1,nva
         nvsa=nvsa+int_mb(k_range+noab+i-1)
        enddo
        do i=1,nvb
         nvsb=nvsb+int_mb(k_range+noab+nva+i-1)
        enddo
c reference string 
        do i=1,nosa
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//keta
           refstring(i+2)=cre5           
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//keta
           refstring(i+2)=cre6
          endif
        enddo
c
        do i=1,nosb
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//ketb
           refstring(i+nosa+2)=cre5
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//ketb
           refstring(i+nosa+2)=cre6
          endif
        enddo
c
c ------
         write(LuOut,9000) 'Singles strings'
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            spinsum=int_mb(k_spin+p1b-1)+int_mb(k_spin+h2b-1)
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1          .eq.irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
                if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1            call errquit('tce_print_x1: MA problem',0,MA_ERR)
                call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1            int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h2 = 1,int_mb(k_range+h2b-1)
                    i = i + 1
                    if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                     iposh=int_mb(k_offset+h2b-1)+h2
                     iposp=int_mb(k_offset+p1b-1)+p1-nosb 
c
                 if(iposh.le.9) then
                  write(num1,301) iposh
                  cre4=bra//num1//kea
                  refstring(2)=cre4     
                 else 
                  write(num2,302) iposh
                  cre5=bra//num2//kea
                  refstring(2)=cre5
                 endif     
c
                 if(iposp.le.9) then
                  write(num1,301) iposp
                  cre5=bra//num1//keta
                  refstring(1)=cre5
                 else
                  write(num2,302) iposp
                  cre6=bra//num2//keta
                  refstring(1)=cre6
                 endif
c yaml initial_state_suggestion state 
c yaml "Singles contribution"
        write(6,*) xnorm*dbl_mb(k_r1+i-1),':',
     &             (refstring(j),j=1,nosa+nosb+2),
     &             '|0>'
c spinsum.eq.4 - mirroring beta->beta excitation
             if(spinsum.eq.2) then  
                 if(iposh.le.9) then
                  write(num1,301) iposh
                  cre4=bra//num1//keb
                  refstring(2)=cre4
                 else
                  write(num2,302) iposh
                  cre5=bra//num2//keb
                  refstring(2)=cre5
                 endif
c
                 if(iposp.le.9) then
                  write(num1,301) iposp
                  cre5=bra//num1//ketb
                  refstring(1)=cre5
                 else
                  write(num2,302) iposp
                  cre6=bra//num2//ketb
                  refstring(1)=cre6
                 endif
c yaml initial_state_suggestion state
c yaml "Singles contribution"
        write(6,*) xnorm*dbl_mb(k_r1+i-1),':',
     &             (refstring(j),j=1,nosa+nosb+2),
     &             '|0>'
             endif  !spinsum
c ------
c
c
                    endif
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r1))
     1            call errquit('tce_print_x1: MA problem',1,MA_ERR)
              endif
            endif
          endif
        enddo
      enddo
      endif
      call ga_sync()
  301 format(i1)
  302 format(i2)
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
c
c
      subroutine tce_string_x1_norm(d_r1,k_r1_offset,r1_thresh,irrep,
     &           xnorm1)
c
c $Id: tce_print_x1.F 26707 2015-01-17 01:06:43Z edo $
c
c Prints the X1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep
c string
       integer spinsum,j
c ------
      logical nodezero
      double precision r1_thresh,xnorm1
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
c ------
      data spinnames/'alpha','beta '/
c xnorm defined only on nodezero
      xnorm1=0.0d0
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
c
         write(LuOut,9000) 'singles norm'
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            spinsum=int_mb(k_spin+p1b-1)+int_mb(k_spin+h2b-1)
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1          .eq.irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
                if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1            call errquit('tce_print_x1: MA problem',0,MA_ERR)
                call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1            int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h2 = 1,int_mb(k_range+h2b-1)
                    i = i + 1
                    if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
c
c string for RHF reference only  a->a 
                 xnorm1=xnorm1+dbl_mb(k_r1+i-1)*dbl_mb(k_r1+i-1)
c spinsum.eq.4 - mirroring beta->beta excitation
             if(spinsum.eq.2) then  
                 xnorm1=xnorm1+dbl_mb(k_r1+i-1)*dbl_mb(k_r1+i-1)
             endif  !spinsum
c ------
c
c
                    endif
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r1))
     1            call errquit('tce_print_x1: MA problem',1,MA_ERR)
              endif
            endif
          endif
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
c
c
c
      subroutine tce_string_ref(xnorm)
c
c $Id: tce_print_x1.F 26707 2015-01-17 01:06:43Z edo $
c
c Prints the X1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer i
c string
       integer spinsum,j,iposh,iposp
       integer nosa,nosb,nvsa,nvsb
c ------
      logical nodezero
      double precision xnorm
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
c string 
      character*6 refstring(202)
      character*1 bra,num1
      character*2 num2
      character*3 keta,ketb,kea,keb
      character*4 cre4
      character*5 cre5
      character*6 cre6
c ------
      data spinnames/'alpha','beta '/
c
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
c        write(6,*)'check 2: xnorm =',xnorm
c        call util_flush(6)
c
c string
c
        bra='('
        keta='a)+'
        ketb='b)+'
        kea='a)'
        keb='b)'
c
        nosa=0
        nosb=0
        nvsa=0
        nvsb=0
        do i=1,noa
         nosa=nosa+int_mb(k_range+i-1)
        enddo
        do i=1,nob
         nosb=nosb+int_mb(k_range+noa+i-1)
        enddo
        do i=1,nva
         nvsa=nvsa+int_mb(k_range+noab+i-1)
        enddo
        do i=1,nvb
         nvsb=nvsb+int_mb(k_range+noab+nva+i-1)
        enddo
c reference string 
        do i=1,nosa
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//keta
           refstring(i+2)=cre5           
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//keta
           refstring(i+2)=cre6
          endif
        enddo
c
        do i=1,nosb
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//ketb
           refstring(i+nosa+2)=cre5
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//ketb
           refstring(i+nosa+2)=cre6
          endif
        enddo
c
         write(LuOut,9000) 'Reference string'
c yaml initial_state_suggestion state
         write(6,*) xnorm,':',
     &             (refstring(j),j=3,nosa+nosb+2),
     &             '|0>'
      endif !nodezero
      call ga_sync()
  301 format(i1)
  302 format(i2)
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
c
c
