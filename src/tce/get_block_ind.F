      subroutine get_block_ind(d_file,array,size,key_ix,indexc)
c      subroutine get_block_ind(d_file,array,size,offset,indexc)
c
c $Id: get_block_ind.F,v 1.8 2008-09-17 03:46:20 jhammond Exp $
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
#include "tce_main.fh"
#include "sym.fh"
c      integer dummy
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
c--      logical parallel
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
c --- keys to orbital part
      integer key_ix,nkey,first_h,second_h
      integer i1,j1,k1,l1
      integer inx1,inx2,inx3,inx4
c
c         print*,'tce_sort version get_block_ind'
c
      cpusecs(5) = cpusecs(5) - util_cpusec()
      cpusecs(55) = cpusecs(55) - util_wallsec()
c--      parallel = (ga_nnodes().gt.1)
cccx      if (util_print('get_block',print_debug)) then
cccx         write(LuOut,9000) ga_nodeid(),d_file,size,offset
cccx      endif
ccccx      if (offset .eq. -1)
ccccx     1  call errquit('get_block: illegal offset',0,UNKNOWN_ERR)
c
c
c
c
c keys to orbital integrals!
c
c
      DO g3b = 1,noab+nvab
      DO g4b = g3b,noab+nvab
      DO g1b = 1,noab+nvab
      DO g2b = g1b,noab+nvab
      IF (int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1) .eq. int_mb(k_spin+g
     &1b-1)+int_mb(k_spin+g2b-1)) THEN
      IF (ieor(int_mb(k_sym+g3b-1),ieor(int_mb(k_sym+g4b-1),ieor(int_mb(
     &k_sym+g1b-1),int_mb(k_sym+g2b-1)))) .eq. irrep_v) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1
     &)+int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1).ne.8)) THEN
c
      nkey=g2b-1+(noab+nvab)*(g1b-1+
     &(noab+nvab)*(g4b-1+(noab+nvab)*(g3b-1)))
c
      if(key_ix.eq.nkey) then !----
       inx3=g3b
       inx4=g4b
       inx1=g1b
       inx2=g2b
c first half
        j1 = int_mb(k_b2am+g3b-1)
        i1 = int_mb(k_b2am+g1b-1)
        l1 = int_mb(k_b2am+g4b-1)
        k1 = int_mb(k_b2am+g2b-1)
        if(i1.ge.j1) then
         i = i1
         j = j1
        else
         i = j1
         j = i1
        end if
        if(k1.ge.l1) then
         k = k1
         l = l1
        else
         k = l1
         l = k1
        end if
         irow = INDEX_PAIR(i,j)
         icol = INDEX_PAIR(k,l)
         if(irow.ge.icol) then
          ires = INDEX_PAIR(irow,icol)
         else
          ires = INDEX_PAIR(icol,irow)
         end if
         first_h = IRES
        j1 = int_mb(k_b2am+g3b-1)
        i1 = int_mb(k_b2am+g2b-1)
        l1 = int_mb(k_b2am+g4b-1)
        k1 = int_mb(k_b2am+g1b-1)
        if(i1.ge.j1) then
         i = i1
         j = j1
        else
         i = j1
         j = i1
        end if
        if(k1.ge.l1) then
         k = k1
         l = l1
        else
         k = l1
         l = k1
        end if
         irow = INDEX_PAIR(i,j)
         icol = INDEX_PAIR(k,l)
         if(irow.ge.icol) then
          ires = INDEX_PAIR(irow,icol)
         else
          ires = INDEX_PAIR(icol,irow)
         end if
         second_h = ires
c ---------
         go to 1000
      end if ! -----
c
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
c
 1000 CONTINUE
c
c
         irow = 0
         icol = 0
         ires = 0
         g1b  = 0
         g2b  = 0
         g3b  = 0
         g4b  = 0
c
c
c
      if (ioalg.eq.2) then
c v^{ g3b < g4b}_{ g1b < g2b } => ( g3b g1b | g4b g2b ) - ( g3b g2b | g4b g1b )
        uaadaa=.false.
        ubbdbb=.false.
        uabdab=.false.
        ubadba=.false.
        ubadab=.false.
        uabdba=.false.
c--        length = int_mb(k_o2ot_offset)
c--        g3b = int_mb(k_o2ot_offset+indexc)
c--        g4b = int_mb(k_o2ot_offset+length+indexc)
c--        g1b = int_mb(k_o2ot_offset+2*length+indexc)
c--        g2b = int_mb(k_o2ot_offset+3*length+indexc)
        g3b=inx3
        g4b=inx4
        g1b=inx1
        g2b=inx2
c CHECK POINT
        if((g3b.ne.inx3).or.(g4b.ne.inx4).or.(g1b.ne.inx1).or.
     &     (g2b.ne.inx2)) then
           call errquit('get_block_ind:g-ind failed',0,UNKNOWN_ERR)
        end if
c
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
c--         key_alpha = int_mb(k_o2o_offset+indexc)
        key_alpha = first_h
c CHECK POINT
         if(key_alpha.ne.first_h) then
          write(6,*)'key_alpha first_h',key_alpha,first_h
          call util_flush(6)
          call errquit('get_block_ind:first_h',0,UNKNOWN_ERR)
         end if
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
c           array(i) = 0.0d0
         enddo
         do i = 1, size
c           dbl_mb(k_a+i-1) = 0.0d0
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
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g2-1)*size4+g4-1)*size1+g1-1)*size3+g3-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 3142
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
         go to 200
        end if
c --- ( g3 g1 | g2 g4 )
        if((.not.lp31p42).and.(.not.l31s).and.(l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g4-1)*size2+g2-1)*size1+g1-1)*size3+g3-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 3124
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
         go to 200
        end if
c --- ( g1 g3 | g4 g2 )
        if((.not.lp31p42).and.(l31s).and.(.not.l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g2-1)*size4+g4-1)*size3+g3-1)*size1+g1-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 1342
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
         go to 200
        end if
c --- ( g1 g3 | g2 g4 )
        if((.not.lp31p42).and.(l31s).and.(l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g4-1)*size2+g2-1)*size3+g3-1)*size1+g1-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 1324
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
         go to 200
        end if
c --- ( g4 g2 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(.not.l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g1-1)*size3+g3-1)*size2+g2-1)*size4+g4-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 4231
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
         go to 200
        end if
c --- ( g4 g2 | g1 g3 )
        if((lp31p42).and.(l31s).and.(.not.l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g3-1)*size1+g1-1)*size2+g2-1)*size4+g4-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 4213
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
         go to 200
        end if
c --- ( g2 g4 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g1-1)*size3+g3-1)*size4+g4-1)*size2+g2-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 2431
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
         go to 200
        end if
c --- ( g2 g4 | g1 g3 )
        if((lp31p42).and.(l31s).and.(l42s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) +
c     1    dbl_mb(k_a+(((g3-1)*size1+g1-1)*size4+g4-1)*size2+g2-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 2413
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
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
c--         key_alpha = int_mb(k_o2o_offset+indexc+length)
        key_alpha = second_h
c CHECK POINT
         if(key_alpha.ne.second_h) then
          write(6,*)'key_alpha second_h',key_alpha,second_h
          call util_flush(6)
          call errquit('get_block_ind:second_h',0,UNKNOWN_ERR)
         end if
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
c--- not opened yet            dbl_mb(k_a+i-1) = 0.0d0
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
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g1-1)*size4+g4-1)*size2+g2-1)*size3+g3-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 3241
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
         go to 300
        end if
c --- ( g3 g2 | g1 g4 )
        if((.not.lp32p41).and.(.not.l32s).and.(l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g4-1)*size1+g1-1)*size2+g2-1)*size3+g3-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 3214
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
         go to 300
        end if
c --- ( g2 g3 | g4 g1 )
        if((.not.lp32p41).and.(l32s).and.(.not.l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g1-1)*size4+g4-1)*size3+g3-1)*size2+g2-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 2341
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
         go to 300
        end if
c --- ( g2 g3 | g1 g4 )
        if((.not.lp32p41).and.(l32s).and.(l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g4-1)*size1+g1-1)*size3+g3-1)*size2+g2-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 2314
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
         go to 300
        end if
c --- ( g4 g1 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(.not.l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g2-1)*size3+g3-1)*size1+g1-1)*size4+g4-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 4132
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
         go to 300
        end if
c --- ( g4 g1 | g2 g3 )
        if((lp32p41).and.(l32s).and.(.not.l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g3-1)*size2+g2-1)*size1+g1-1)*size4+g4-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 4123
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
         go to 300
        end if
c --- ( g1 g4 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g2-1)*size3+g3-1)*size4+g4-1)*size1+g1-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 1432
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
         go to 300
        end if
c --- ( g1 g4 | g2 g3 )
        if((lp32p41).and.(l32s).and.(l41s)) then
c         i = 0
c         do g3 = 1,int_mb(k_range+g3b-1)
c         do g4 = 1,int_mb(k_range+g4b-1)
c         do g1 = 1,int_mb(k_range+g1b-1)
c         do g2 = 1,int_mb(k_range+g2b-1)
c         i=i+1
c          array(i) = array(i) -
c     1    dbl_mb(k_a+(((g3-1)*size2+g2-1)*size4+g4-1)*size1+g1-1)
c         enddo
c         enddo
c         enddo
c         enddo
c         dummy = 1423
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
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
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
      subroutine get_block_ind_i(d_file,array,size,key_ix,indexc,
     &                         w2b,w1b,w4b,w3b)
c      subroutine get_block_ind(d_file,array,size,offset,indexc)
c
c $Id: get_block_ind.F,v 1.8 2008-09-17 03:46:20 jhammond Exp $
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
#include "tce_main.fh"
#include "sym.fh"
c      integer dummy
      integer d_file
cccx      integer d_f
      integer size
cccx      integer offset
      double precision array(size)
c#ifdef TCESORT
c      double precision array2(size)
c#endif
      integer request
      integer l_a,k_a
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer w1b,w2b,w3b,w4b
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      integer length
      integer offset_alpha,key_alpha
      logical used
c--      logical parallel
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
c --- keys to orbital part
      integer key_ix,nkey,first_h,second_h
      integer i1,j1,k1,l1
      integer inx1,inx2,inx3,inx4
      external index_pair
c
c         print*,'tce_sort version get_block_ind_i'
c
c
c
c      cpusecs(5) = cpusecs(5) - util_cpusec()
c      cpusecs(55) = cpusecs(55) - util_wallsec()
c
c
c
c
c keys to orbital integrals!
c
c *** debug ***
c       write(6,*)'from gbbi g-s',w2b,w1b,w4b,w3b
c       write(6,*)'from gbbi d_file',d_file
c       call util_flush(6)
c ***************
c
       inx3=w3b
       inx4=w4b
       inx1=w1b
       inx2=w2b
        j1 = int_mb(k_b2am+w3b-1)
        i1 = int_mb(k_b2am+w1b-1)
        l1 = int_mb(k_b2am+w4b-1)
        k1 = int_mb(k_b2am+w2b-1)
        if(i1.ge.j1) then
         i = i1
         j = j1
        else
         i = j1
         j = i1
        end if
        if(k1.ge.l1) then
         k = k1
         l = l1
        else
         k = l1
         l = k1
        end if
         irow = INDEX_PAIR(i,j)
         icol = INDEX_PAIR(k,l)
         if(irow.ge.icol) then
          ires = INDEX_PAIR(irow,icol)
         else
          ires = INDEX_PAIR(icol,irow)
         end if
         first_h = IRES
        j1 = int_mb(k_b2am+w3b-1)
        i1 = int_mb(k_b2am+w2b-1)
        l1 = int_mb(k_b2am+w4b-1)
        k1 = int_mb(k_b2am+w1b-1)
        if(i1.ge.j1) then
         i = i1
         j = j1
        else
         i = j1
         j = i1
        end if
        if(k1.ge.l1) then
         k = k1
         l = l1
        else
         k = l1
         l = k1
        end if
         irow = INDEX_PAIR(i,j)
         icol = INDEX_PAIR(k,l)
         if(irow.ge.icol) then
          ires = INDEX_PAIR(irow,icol)
         else
          ires = INDEX_PAIR(icol,irow)
         end if
         second_h = ires
c
c
         irow = 0
         icol = 0
         ires = 0
         g1b  = 0
         g2b  = 0
         g3b  = 0
         g4b  = 0
c
c
c
      if (ioalg.eq.2) then
c v^{ g3b < g4b}_{ g1b < g2b } => ( g3b g1b | g4b g2b ) - ( g3b g2b | g4b g1b )
        uaadaa=.false.
        ubbdbb=.false.
        uabdab=.false.
        ubadba=.false.
        ubadab=.false.
        uabdba=.false.
c--        length = int_mb(k_o2ot_offset)
c--        g3b = int_mb(k_o2ot_offset+indexc)
c--        g4b = int_mb(k_o2ot_offset+length+indexc)
c--        g1b = int_mb(k_o2ot_offset+2*length+indexc)
c--        g2b = int_mb(k_o2ot_offset+3*length+indexc)
        g3b=inx3
        g4b=inx4
        g1b=inx1
        g2b=inx2
c CHECK POINT
        if((g3b.ne.inx3).or.(g4b.ne.inx4).or.(g1b.ne.inx1).or.
     &     (g2b.ne.inx2)) then
           call errquit('get_block_ind:g-ind failed',0,UNKNOWN_ERR)
        end if
c
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
c--         key_alpha = int_mb(k_o2o_offset+indexc)
        key_alpha = first_h
c CHECK POINT
         if(key_alpha.ne.first_h) then
          write(6,*)'key_alpha first_h',key_alpha,first_h
          call util_flush(6)
          call errquit('get_block_ind:first_h',0,UNKNOWN_ERR)
         end if
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
c *** debug ***
c       write(6,*)'from gbbi before ma push get'
c       call util_flush(6)
c ***************
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
c           array(i) = 0.0d0
         enddo
         do i = 1, size
c           dbl_mb(k_a+i-1) = 0.0d0
           array(i) = 0.0d0
         enddo
c *** debug ***
c       write(6,*)'from gbbi d_v2orb 1',d_v2orb
c       call util_flush(6)
c ***************
         call ga_get(d_v2orb,offset_alpha+1,offset_alpha+size,1,1,
     1         dbl_mb(k_a),1)
c *** debug ***
c       write(6,*)'from gbbi d_v2orb 1',d_v2orb
c       call util_flush(6)
c ***************
c
         size1 = int_mb(k_range+g1b-1)
         size2 = int_mb(k_range+g2b-1)
         size3 = int_mb(k_range+g3b-1)
         size4 = int_mb(k_range+g4b-1)
c
c
c --- ( g3 g1 | g4 g2 )
        if((.not.lp31p42).and.(.not.l31s).and.(.not.l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 3142000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
c#endif
         go to 200
        end if
c --- ( g3 g1 | g2 g4 )
        if((.not.lp31p42).and.(.not.l31s).and.(l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 3124000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g1 g3 | g4 g2 )
        if((.not.lp31p42).and.(l31s).and.(.not.l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 1342000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g1 g3 | g2 g4 )
        if((.not.lp31p42).and.(l31s).and.(l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 1324000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g4 g2 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(.not.l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 4231000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g4 g2 | g1 g3 )
        if((lp31p42).and.(l31s).and.(.not.l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 4213000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g2 g4 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 2431000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 200
        end if
c --- ( g2 g4 | g1 g3 )
        if((lp31p42).and.(l31s).and.(l42s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 2413000
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
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
c *** debug ***
c       write(6,*)'from gbbi second_half'
c       call util_flush(6)
c ***************
c
c
c
        if (uaadaa.or.ubbdbb.or.uabdba.or.ubadab) then
c *** debug ***
c       write(6,*)'from gbbi in sec_half'
c       call util_flush(6)
c ***************
c second half
c--         key_alpha = int_mb(k_o2o_offset+indexc+length)
        key_alpha = second_h
c CHECK POINT
         if(key_alpha.ne.second_h) then
          write(6,*)'key_alpha second_h',key_alpha,second_h
          call util_flush(6)
          call errquit('get_block_ind:second_h',0,UNKNOWN_ERR)
         end if
c *** debug ***
c       write(6,*)'from gbbi in sec_half s1'
c       call util_flush(6)
c ***************
c defining the order
         ig2b=int_mb(k_b2am+g2b-1)
         ig1b=int_mb(k_b2am+g1b-1)
         ig3b=int_mb(k_b2am+g3b-1)
         ig4b=int_mb(k_b2am+g4b-1)
c *** debug ***
c       write(6,*)'from gbbi in sec_half s2'
c       write(6,*)'g2b g1b g3b g4b',g2b,g1b,g3b,g4b
c       write(6,*)'ig2b ig1b ig3b ig4b',ig2b,ig1b,ig3b,ig4b
c       call util_flush(6)
c ***************
         if(ig3b.ge.ig2b) then
c *** debug ***
c         write(6,*)'sh2 in1'
c         call util_flush(6)
c *************
            l32s=.false.
            irow=index_pair(ig3b,ig2b)
c *** debug ***
c        write(6,*)'from gbbi sh 1 irow',ig3b,ig2b,irow
c *************
         else
c *** debug ***
c         write(6,*)'sh2 in2'
c         call util_flush(6)
c *************
            l32s=.true.
            irow=index_pair(ig2b,ig3b)
c *** debug ***
c        write(6,*)'from gbbi sh 2 irow',ig2b,ig3b,irow
c *************
         end if
         if(ig4b.ge.ig1b) then
c *** debug ***
c         write(6,*)'sh2 in3'
c         call util_flush(6)
c *************
            l41s=.false.
            icol=index_pair(ig4b,ig1b)
c *** debug ***
c        write(6,*)'from gbbi sh 3 irow',ig4b,ig1b,icol
c *************
         else
c *** debug ***
c         write(6,*)'sh2 in4'
c         call util_flush(6)
c *************
            l41s=.true.
            icol=index_pair(ig1b,ig4b)
c *** debug ***
c        write(6,*)'from gbbi sh 3 irow',ig1b,ig4b,icol
c *************
         end if
         if(irow.ge.icol) then
            lp32p41=.false.
            check_key=index_pair(irow,icol)
         else
            lp32p41=.true.
            check_key=index_pair(icol,irow)
         end if
c *** debug ***
c       write(6,*)'from gbbi in sec_half s2a'
c       write(6,*)'from gbbi in seci_half size',size
c       write(6,*)'from gbbi in seci_half k_a',k_a
c       call util_flush(6)
c ***************
c
         if(uabdba.or.ubadab) then
          do i = 1, size
c--- k_a is not opened          dbl_mb(k_a+i-1) = 0.0d0
            array(i) = 0.0d0
          enddo
         end if
c *** debug ***
c       write(6,*)'from gbbi in sec_half s3'
c       call util_flush(6)
c ***************
c
         call tce_hash_n(int_mb(k_v2_alpha_offset),
     1                 key_alpha,offset_alpha)
c *** debug ***
c       write(6,*)'from gbbi in sec_half s4'
c       call util_flush(6)
c ***************
c check_key test
         if(check_key.ne.key_alpha)
     1   call errquit('get_block_ind:check_key2',0,UNKNOWN_ERR)
c -------------
c *** debug ***
c       write(6,*)'from gbbi before second ma push get'
c       call util_flush(6)
c ***************
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
c *** debug ***
c       write(6,*)'from gbbi d_v2orb 2',d_v2orb
c       call util_flush(6)
c ***************
         call ga_get(d_v2orb,offset_alpha+1,offset_alpha+size,1,1,
     1         dbl_mb(k_a),1)
c *** debug ***
c       write(6,*)'after gbbi d_v2orb 2'
c       call util_flush(6)
c ***************
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
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 3241001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g3 g2 | g1 g4 )
        if((.not.lp32p41).and.(.not.l32s).and.(l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 3214001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g2 g3 | g4 g1 )
        if((.not.lp32p41).and.(l32s).and.(.not.l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 2341001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g2 g3 | g1 g4 )
        if((.not.lp32p41).and.(l32s).and.(l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
c         call zero_array(array,size,size1,size2,size3,size4)
c         call zero_array(array2,size,size1,size2,size3,size4)
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
         dummy = 2314001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g4 g1 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(.not.l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 4132001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g4 g1 | g2 g3 )
        if((lp32p41).and.(l32s).and.(.not.l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 4123001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g1 g4 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 1432001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
         go to 300
        end if
c --- ( g1 g4 | g2 g3 )
        if((lp32p41).and.(l32s).and.(l41s)) then
c         call copy_array(array2,array,size,size1,size2,size3,size4)
#ifdef LOOPS
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
         dummy = 1423001
#endif
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
c         call test_array(array2,array,size,size1,size2,size3,size4,
c          &                    dummy)
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
c *** debug ***
c       write(6,*)'leaving gbbi g-s'
c       call util_flush(6)
c ***************
c
c      cpusecs(5) = cpusecs(5) + util_cpusec()
c      cpusecs(55) = cpusecs(55) + util_wallsec()
      return
      end

#ifdef LOOPS

      subroutine copy_array(array2,array,size,size1,size2,size3,size4)
      implicit none
c#include "util.fh"
c#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
      integer size,i
      double precision array(size)
      double precision array2(size)
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4

         i = 0
         do g3 = 1,size3
         do g4 = 1,size4
         do g1 = 1,size1
         do g2 = 1,size2
         i=i+1
          array2(i) = array(i)
         enddo
         enddo
         enddo
         enddo

      return
      end

      subroutine test_array(array2,array,size,size1,size2,size3,size4,
     &                      dummy)
      implicit none
c#include "util.fh"
c#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
      integer size,dummy,i
      double precision array(size)
      double precision array2(size)
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      logical throw

         throw=.false.
         i = 0
         do g3 = 1,size3
         do g4 = 1,size4
         do g1 = 1,size1
         do g2 = 1,size2
         i=i+1
          if ( (abs(array(i)).gt.(1.0d-9)).or.
     &         (abs(array2(i)).gt.(1.0d-9)) ) then
          if (array(i).ne.array2(i)) then
            print*,'PROBLEM',i,array(i),array2(i)
            throw=.true.
          endif
          endif
         enddo
         enddo
         enddo
         enddo
         if (throw)
     &     call errquit('get_block_ind: test failed',dummy,UNKNOWN_ERR)

      return
      end

      subroutine zero_array(array,size,size1,size2,size3,size4)
      implicit none
c#include "util.fh"
c#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
      integer size,i
      double precision array(size)
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4

         i = 0
         do g3 = 1,size3
         do g4 = 1,size4
         do g1 = 1,size1
         do g2 = 1,size2
         i=i+1
          array(i) = 0.0d0
         enddo
         enddo
         enddo
         enddo

      return
      end

#endif
