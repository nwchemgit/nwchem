!
! History: At one point, all the sortacc4 routines were explicit loops.
!          These were carefully eliminated in favor of calls to the
!          potentially optimized subroutine.  The loops remain in the
!          comments since back-translating - if ever necessary - would
!          be tedious otherwise.
!                                   - Jeff Hammond, 2015
!
      subroutine get_block_ind(d_file,array,size,key_ix,indexc)
c
c $Id$
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
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f_a(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_a
#endif
      integer alloc_error
#else
      integer l_a,k_a
#endif
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      integer length
      integer off_a,key_alpha ! off_a was formerly offset_alpha
      logical used
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
ccx         call tce_hash_n(int_mb(k_v2_alpha_offset),key_alpha,off_a)
c - exascale ---
        if(model.eq.'ccsd_act') then
         call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        end if
c ---
c check_key test
ccx         if(check_key.ne.key_alpha)
ccx     1   call errquit('get_block_ind:check_key1',0,UNKNOWN_ERR)
c
c
c
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         do i = 1, size
           f_a(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
         do i = 1, size
           array(i) = 0.0d0
         enddo
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#endif
         go to 200
        end if
c
 200     continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
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
ccx         call tce_hash_n(int_mb(k_v2_alpha_offset),key_alpha,off_a)
c - exascale ---
        if(model.eq.'ccsd_act') then
         call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        end if
c ---
c check_key test
ccx         if(check_key.ne.key_alpha)
ccx     1   call errquit('get_block_ind:check_key2',0,UNKNOWN_ERR)
c -------------
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         do i = 1, size
           f_a(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#endif
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
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#endif
         go to 300
        end if
c
c
c
 300    continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
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
      subroutine get_block_ind_i(d_file,array,size,key_ix,indexc,
     &                         w2b,w1b,w4b,w3b)
c
c $Id$
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
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f_a(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_a
#endif
      integer alloc_error
#else
      integer l_a,k_a
#endif
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer w1b,w2b,w3b,w4b
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      integer length
      integer off_a,key_alpha ! off_a was formerly offset_alpha
      logical used
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
c - exascale ---
cc      LOGICAL is_active_1,is_active_2,is_active_3,is_active_4
c --- keys to orbital part
      integer key_ix,nkey,first_h,second_h
      integer i1,j1,k1,l1
      integer inx1,inx2,inx3,inx4
      external index_pair
c
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
           first_h= k - 1 + (noa+nva) * (l - 1 + (noa+nva) * (i-
     &               1 + (noa+nva) * (j - 1)))
         else
           first_h= i - 1 + (noa+nva) * (j - 1 + (noa+nva) * (k-
     &               1 + (noa+nva) * (l - 1)))
         end if
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
          second_h= k - 1 + (noa+nva) * (l - 1 + (noa+nva) * (i-
     &               1 + (noa+nva) * (j - 1)))
         else
          second_h= i - 1 + (noa+nva) * (j - 1 + (noa+nva) * (k-
     &               1 + (noa+nva) * (l - 1)))
         end if
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
        key_alpha = first_h
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
c - exascale ---
       if(model.eq.'ccsd_act') then
        call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       end if
c ---
c check_key test
cc         if(check_key.ne.key_alpha)
cc     1   call errquit('get_block_ind:check_key1',0,UNKNOWN_ERR)
c
c
c
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         do i = 1, size
           array(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
c --- ACHTUNG ATTENTION UWAGA BHUMAHE
ccx           dbl_mb(k_a+i-1) = 0.0d0
           array(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
c
         size1 = int_mb(k_range+g1b-1)
         size2 = int_mb(k_range+g2b-1)
         size3 = int_mb(k_range+g3b-1)
         size4 = int_mb(k_range+g4b-1)
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
c         dummy = 3142000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#endif
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
c         dummy = 3124000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#endif
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
c         dummy = 1342000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#endif
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
c         dummy = 1324000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#endif
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
c         dummy = 4231000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#endif
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
c         dummy = 4213000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#endif
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
c         dummy = 2431000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#endif
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
c         dummy = 2413000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#endif
         go to 200
        end if
c
 200     continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
         end if !spin cases
c
c
c
c
        if (uaadaa.or.ubbdbb.or.uabdba.or.ubadab) then
c second half
        key_alpha = second_h
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
c
         if(uabdba.or.ubadab) then
          do i = 1, size
c--- k_a is not opened          dbl_mb(k_a+i-1) = 0.0d0
            array(i) = 0.0d0
          enddo
         end if
c
c - exascale ---
       if(model.eq.'ccsd_act') then
         call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       end if
c ---
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
c --- ACHTUNG ATTENTION UWAGA BHUMAHE
c         do i = 1, size
c           dbl_mb(k_a+i-1) = 0.0d0
c         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
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
c         dummy = 3241001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#endif
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
c         dummy = 3214001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#endif
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
c         dummy = 2341001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#endif
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
c         dummy = 2314001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#endif
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
c         dummy = 4132001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#endif
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
c         dummy = 4123001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#endif
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
c         dummy = 1432001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#endif
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
c         dummy = 1423001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#endif
         go to 300
        end if
c
c
c
 300    continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
         end if !spin cases
c
c
c
      end if    !ioalg
c
c
c      cpusecs(5) = cpusecs(5) + util_cpusec()
c      cpusecs(55) = cpusecs(55) + util_wallsec()
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
      subroutine get_block_ind_i_przerobiony(
     &d_file,array,size,key_ix,indexc,
     &                         w2b,w1b,w4b,w3b)
c      subroutine get_block_ind(d_file,array,size,offset,indexc)
c
c $Id$
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
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f_a(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_a
#endif
      integer alloc_error
#else
      integer l_a,k_a
#endif
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer w1b,w2b,w3b,w4b
      integer g1,g2,g3,g4
      integer size1,size2,size3,size4
      integer length
      integer off_a,key_alpha
      logical used
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
ccx         call tce_hash_n(int_mb(k_v2_alpha_offset),key_alpha,off_a)
c - exascale ---
       if(model.eq.'ccsd_act') then
        call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       end if
c ---
c check_key test
ccx         if(check_key.ne.key_alpha)
ccx     1   call errquit('get_block_ind:check_key1',0,UNKNOWN_ERR)
c
c
c
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         do i = 1, size
           f_a(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
         do i = 1, size
           array(i) = 0.0d0
         enddo
c
         size1 = int_mb(k_range+g1b-1)
         size2 = int_mb(k_range+g2b-1)
         size3 = int_mb(k_range+g3b-1)
         size4 = int_mb(k_range+g4b-1)
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
c         dummy = 3142000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
#endif
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
c         dummy = 3124000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
#endif
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
c         dummy = 1342000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
#endif
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
c         dummy = 1324000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
#endif
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
c         dummy = 4231000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
#endif
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
c         dummy = 4213000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
#endif
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
c         dummy = 2431000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
#endif
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
c         dummy = 2413000
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
#endif
         go to 200
        end if
c
 200     continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
         end if !spin cases
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
c--- k_a is not opened          dbl_mb(k_a+i-1) = 0.0d0
            array(i) = 0.0d0
          enddo
         end if
c
ccx         call tce_hash_n(int_mb(k_v2_alpha_offset),key_alpha,off_a)
c - exascale ---
       if(model.eq.'ccsd_act') then
         call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
       end if
c ---
c check_key test
ccx         if(check_key.ne.key_alpha)
ccx     1   call errquit('get_block_ind:check_key2',0,UNKNOWN_ERR)
c -------------
#ifdef USE_F90_ALLOCATABLE
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind1',2,MA_ERR)
         endif
         do i = 1, size
           f_a(i) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
#else
         IF(.not.MA_PUSH_GET(mt_dbl,size,'noname',l_a,k_a))
     1           CALL ERRQUIT('get-block-ind1',2,MA_ERR)
         do i = 1, size
           dbl_mb(k_a+i-1) = 0.0d0
         enddo
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,dbl_mb(k_a),size)
#endif
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
c         dummy = 3241001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
#endif
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
c         dummy = 3214001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
#endif
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
c         dummy = 2341001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
#endif
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
c         dummy = 2314001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
#endif
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
c         dummy = 4132001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
#endif
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
c         dummy = 4123001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
#endif
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
c         dummy = 1432001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
#endif
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
c         dummy = 1423001
#ifdef USE_F90_ALLOCATABLE
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#else
         call tce_sortacc_4(dbl_mb(k_a),array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
#endif
         go to 300
        end if
c
c
c
 300    continue
c
#ifdef USE_F90_ALLOCATABLE
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get-block-ind2',3,MA_ERR)
         endif
#else
         IF(.not.MA_POP_STACK(l_a))
     1            CALL ERRQUIT('get-block-ind2',3,MA_ERR)
#endif
         end if !spin cases
c
c
c
      end if    !ioalg
c
c      cpusecs(5) = cpusecs(5) + util_cpusec()
c      cpusecs(55) = cpusecs(55) + util_wallsec()
      return
      end

#if 0

      subroutine copy_array(array2,array,size,size1,size2,size3,size4)
      implicit none
#include "errquit.fh"
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
#include "errquit.fh"
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
#include "errquit.fh"
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
