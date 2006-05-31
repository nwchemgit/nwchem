      subroutine get_block_ind(d_file,array,size,indexc)
c      subroutine get_block_ind(d_file,array,size,offset,indexc)
c
c $Id: get_block_ind.F,v 1.1 2006-05-31 23:15:11 kowalski Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sf.fh"
#include "eaf.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
      integer d_file
cccx      integer d_f
      integer size
cccx      integer offset
      double precision array(size)
      integer request
      integer l_a,k_a
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      integer length
      integer offset_alpha,key_alpha
      logical used
      logical parallel
      character*255 filename
c --- controls the logic of orbital integrals --
      logical l31s,l42s,lp31p42
      logical l32s,l41s,lp32p41
      integer ig1b,ig2b,ig3b,ig4b
      integer irow,icol,ires
      integer index_pair
      integer check_key
      logical uaadaa, ubbdbb
      logical uabdab,ubadba,ubadab,uabdba 
c ----------------------------------------------
c
      cpusecs(5) = cpusecs(5) - util_cpusec()
      cpusecs(55) = cpusecs(55) - util_wallsec()
      parallel = (ga_nnodes().gt.1)
cccx      if (util_print('get_block',print_debug)) then
cccx         write(LuOut,9000) ga_nodeid(),d_file,size,offset
cccx      endif
ccccx      if (offset .eq. -1)
ccccx     1  call errquit('get_block: illegal offset',0,UNKNOWN_ERR)
c
      if (ioalg.eq.2) then
c v^{ g3b < g4b}_{ g1b < g2b } => ( g3b g1b | g4b g2b ) - ( g3b g2b | g4b g1b )
        uaadaa=.false.
        ubbdbb=.false.
        uabdab=.false.
        ubadba=.false.
        ubadab=.false.
        uabdba=.false. 
        length = int_mb(k_o2ot_offset)
        g3b = int_mb(k_o2ot_offset+indexc)
        g4b = int_mb(k_o2ot_offset+length+indexc)
        g1b = int_mb(k_o2ot_offset+2*length+indexc)
        g2b = int_mb(k_o2ot_offset+3*length+indexc)
        ispin=int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1)+
     1        int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)
        if(ispin.eq.4) uaadaa=.true.
        if(ispin.eq.8) ubbdbb=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     1     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2) 
     2     uabdab=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     1     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     2     ubadba=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     1     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2)
     2     ubadab=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     1     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     2     uabdba=.true.
c        if ((ispin.eq.6).or.(ispin.eq.4).or.(ispin.eq.8)) then 
        if (uaadaa.or.ubbdbb.or.uabdab.or.ubadba) then
c first half
         key_alpha = int_mb(k_o2o_offset+indexc) 
c defining the order
         ig2b=int_mb(k_b2am+g2b-1)
         ig1b=int_mb(k_b2am+g1b-1)
         ig3b=int_mb(k_b2am+g3b-1)
         ig4b=int_mb(k_b2am+g4b-1)
         if(ig3b.ge.ig1b) then
            l31s=.false.
            irow=index_pair(ig3b,ig1b)
         else
            l31s=.true.
            irow=index_pair(ig1b,ig3b) 
         end if
         if(ig4b.ge.ig2b) then
            l42s=.false.
            icol=index_pair(ig4b,ig2b)
         else
            l42s=.true.
            icol=index_pair(ig2b,ig4b)
         end if
         if(irow.ge.icol) then
            lp31p42=.false.
            check_key=index_pair(irow,icol)
         else
            lp31p42=.true.
            check_key=index_pair(icol,irow)   
         end if
         call tce_hash_n(int_mb(k_v2_alpha_offset),
     1                 key_alpha,offset_alpha)
c check_key test
         if(check_key.ne.key_alpha) 
     1   call errquit('get_block_ind:check_key1',0,UNKNOWN_ERR)
c
c
c
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a)) 
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
           array(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,offset_alpha+1,offset_alpha+size,1,1,
     1         dbl_mb(k_a),1)
c
         size3 = int_mb(k_range+g3b-1)
         size1 = int_mb(k_range+g1b-1)
         size4 = int_mb(k_range+g4b-1)
         size2 = int_mb(k_range+g2b-1)
c
c
c --- ( g3 g1 | g4 g2 ) 
        if((.not.lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) + 
     1    dbl_mb(k_a+(((g2-1)*size4+g4-1)*size1+g1-1)*size3+g3-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g3 g1 | g2 g4 )
        if((.not.lp31p42).and.(.not.l31s).and.(l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g4-1)*size2+g2-1)*size1+g1-1)*size3+g3-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g1 g3 | g4 g2 )
        if((.not.lp31p42).and.(l31s).and.(.not.l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g2-1)*size4+g4-1)*size3+g3-1)*size1+g1-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g1 g3 | g2 g4 )
        if((.not.lp31p42).and.(l31s).and.(l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g4-1)*size2+g2-1)*size3+g3-1)*size1+g1-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g4 g2 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g1-1)*size3+g3-1)*size2+g2-1)*size4+g4-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g4 g2 | g1 g3 )
        if((lp31p42).and.(l31s).and.(.not.l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g3-1)*size1+g1-1)*size2+g2-1)*size4+g4-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c --- ( g2 g4 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g1-1)*size3+g3-1)*size4+g4-1)*size2+g2-1)
         enddo
         enddo
         enddo
         enddo
         go to 200 
        end if
c --- ( g2 g4 | g1 g3 )
        if((lp31p42).and.(l31s).and.(l42s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) +
     1    dbl_mb(k_a+(((g3-1)*size1+g1-1)*size4+g4-1)*size2+g2-1)
         enddo
         enddo
         enddo
         enddo
         go to 200
        end if
c
 200     continue
c
         IF(.not.MA_POP_STACK(l_a)) 
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
         end if !spin cases
c
c
c
c
c
c
        if (uaadaa.or.ubbdbb.or.uabdba.or.ubadab) then
c second half
         key_alpha = int_mb(k_o2o_offset+indexc+length) 
c defining the order
         ig2b=int_mb(k_b2am+g2b-1)
         ig1b=int_mb(k_b2am+g1b-1)
         ig3b=int_mb(k_b2am+g3b-1)
         ig4b=int_mb(k_b2am+g4b-1)
         if(ig3b.ge.ig2b) then
            l32s=.false.
            irow=index_pair(ig3b,ig2b)
         else
            l32s=.true.
            irow=index_pair(ig2b,ig3b)
         end if
         if(ig4b.ge.ig1b) then
            l41s=.false.
            icol=index_pair(ig4b,ig1b)
         else
            l41s=.true.
            icol=index_pair(ig1b,ig4b)
         end if
         if(irow.ge.icol) then
            lp32p41=.false.
            check_key=index_pair(irow,icol)
         else
            lp32p41=.true.
            check_key=index_pair(icol,irow)
         end if
c
         if(uabdba.or.ubadab) then
          do i = 1, size
            dbl_mb(k_a+i-1) = 0.0d0
            array(i) = 0.0d0
          enddo
         end if
c
         call tce_hash_n(int_mb(k_v2_alpha_offset),
     1                 key_alpha,offset_alpha)
c check_key test
         if(check_key.ne.key_alpha)
     1   call errquit('get_block_ind:check_key2',0,UNKNOWN_ERR)
c -------------
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a)) 
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
         call ga_get(d_v2orb,offset_alpha+1,offset_alpha+size,1,1,
     1         dbl_mb(k_a),1)
c
c
         size3 = int_mb(k_range+g3b-1)
         size1 = int_mb(k_range+g1b-1)
         size4 = int_mb(k_range+g4b-1)
         size2 = int_mb(k_range+g2b-1)
c
c
c --- ( g3 g2 | g4 g1 )
        if((.not.lp32p41).and.(.not.l32s).and.(.not.l41s)) then  
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g1-1)*size4+g4-1)*size2+g2-1)*size3+g3-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g3 g2 | g1 g4 )
        if((.not.lp32p41).and.(.not.l32s).and.(l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g4-1)*size1+g1-1)*size2+g2-1)*size3+g3-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g2 g3 | g4 g1 )
        if((.not.lp32p41).and.(l32s).and.(.not.l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g1-1)*size4+g4-1)*size3+g3-1)*size2+g2-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g2 g3 | g1 g4 )
        if((.not.lp32p41).and.(l32s).and.(l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g4-1)*size1+g1-1)*size3+g3-1)*size2+g2-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g4 g1 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(.not.l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g2-1)*size3+g3-1)*size1+g1-1)*size4+g4-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g4 g1 | g2 g3 )
        if((lp32p41).and.(l32s).and.(.not.l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g3-1)*size2+g2-1)*size1+g1-1)*size4+g4-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g1 g4 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g2-1)*size3+g3-1)*size4+g4-1)*size1+g1-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c --- ( g1 g4 | g2 g3 )
        if((lp32p41).and.(l32s).and.(l41s)) then
         i = 0
         do g3 = 1,int_mb(k_range+g3b-1)
         do g4 = 1,int_mb(k_range+g4b-1)
         do g1 = 1,int_mb(k_range+g1b-1)
         do g2 = 1,int_mb(k_range+g2b-1)
         i=i+1
          array(i) = array(i) -
     1    dbl_mb(k_a+(((g3-1)*size2+g2-1)*size4+g4-1)*size1+g1-1)
         enddo
         enddo
         enddo
         enddo
         go to 300
        end if
c
c
c
 300    continue
c
         IF(.not.MA_POP_STACK(l_a)) 
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
         end if !spin cases
c
c
c
      end if    !ioalg
c
      cpusecs(5) = cpusecs(5) + util_cpusec()
      cpusecs(55) = cpusecs(55) + util_wallsec()
      return
      end


