      subroutine add_hash_block_r(d_file,array,size,hash,key)
      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      double precision array(size)
      call tce_hash(hash,key,offset)
#ifdef TCE_TIME_COMMS
      cpusecs(4) = cpusecs(4) - util_cpusec()
      cpusecs(54) = cpusecs(54) - util_wallsec()
#endif
      if (util_print('add_block_r',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     &  call errquit('add_block_r: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.2) then
!$omp critical
        call ga_acc(d_file,offset+1,offset+size,1,1,array,1,1.0d0)
!$omp end critical
      else
        call errquit('add_block_r: ioalg=2 only!',ioalg,INPUT_ERR)
      endif
 9000 format(1x,'node',i3,' add_block_r request to file:',i10,
     &  ' size:',i10,' offset:',i10)
#ifdef TCE_TIME_COMMS
      cpusecs(4) = cpusecs(4) + util_cpusec()
      cpusecs(54) = cpusecs(54) + util_wallsec()
#endif
      return
      end

      subroutine add_hash_block_nb_r(d_file,array,size,hash,key,nbh)
      implicit none
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer nbh
      double precision array(size)
      call tce_hash(hash,key,offset)
#ifdef TCE_TIME_COMMS
      cpusecs(4) = cpusecs(4) - util_cpusec()
      cpusecs(54) = cpusecs(54) - util_wallsec()
#endif
      if (util_print('add_block_nb_r',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     &  call errquit('add_block_nb_r: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.2) then
!$omp critical
        call ga_nbacc(d_file,offset+1,offset+size,1,1,array,1,1.0d0,nbh)
!$omp end critical
      else
        call errquit('add_block_nb_r: ioalg=2 only!',ioalg,INPUT_ERR)
      endif
#ifdef TCE_TIME_COMMS
      cpusecs(4) = cpusecs(4) + util_cpusec()
      cpusecs(54) = cpusecs(54) + util_wallsec()
#endif
      return
 9000 format(1x,'node',i3,' add_block_nb_r request to file:',i10,
     &  ' size:',i10,' offset:',i10)
      end

      subroutine put_hash_block_r(d_file,array,size,hash,key)
      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      double precision array(size)
      call tce_hash(hash,key,offset)
#ifdef TCE_TIME_COMMS
      cpusecs(3) = cpusecs(3) - util_cpusec()
      cpusecs(53) = cpusecs(53) - util_wallsec()
#endif
      if (util_print('put_block_r',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     &  call errquit('put_block_r: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.2) then
!$omp critical
        call ga_put(d_file,offset+1,offset+size,1,1,array,size)
!$omp end critical
      else
        call errquit('put_block_r: ioalg=2 only!',ioalg,UNKNOWN_ERR)
      endif
 9000 format(1x,'node',i3,' put_block_r request to file:',i10,
     &  ' size:',i10,' offset:',i10)
#ifdef TCE_TIME_COMMS
      cpusecs(3) = cpusecs(3) + util_cpusec()
      cpusecs(53) = cpusecs(53) + util_wallsec()
#endif
      return
      end

      subroutine get_hash_block_r(d_file,array,size,hash,key)
      implicit none
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer indexc
      integer iauxo
      double precision array(size)
      if (intorb.and.(d_file.eq.d_v2orb)) then
        iauxo=k_v2spin_offset
        call tce_hash_ind(iauxo,key,indexc)
        call get_block_ind_r(d_file,array,size,key,indexc)
      else
        call tce_hash(hash,key,offset)
        call get_block_r(d_file,array,size,offset)
      endif
      return
      end

      subroutine get_hash_block_i_r(d_file,array,size,hash,key,
     &                            g2b,g1b,g4b,g3b)
      implicit none
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer indexc
      integer g2b,g1b,g4b,g3b
      integer qy2,qy1,qy4,qy3
      double precision array(size)
      if (intorb.and.(d_file.eq.d_v2orb)) then
        qy2=g2b
        qy1=g1b
        qy4=g4b
        qy3=g3b
        call get_block_ind_i_r(d_file,array,size,key,indexc,
     &                       qy2,qy1,qy4,qy3)
      else
        call tce_hash(hash,key,offset)
        call get_block_r(d_file,array,size,offset)
      end if
      return
      end

      subroutine get_block_ind_r(d_file,array,size,key_ix,indexc)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "sym.fh"
      integer d_file
      integer size
      double precision array(size)
      double precision, allocatable :: f_a(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_a
#endif
      integer alloc_error
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer size1,size2,size3,size4
      integer off_a,key_alpha ! off_a was formerly offset_alpha
! --- controls the logic of orbital integrals --
      logical l31s,l42s,lp31p42
      logical l32s,l41s,lp32p41
      integer ig1b,ig2b,ig3b,ig4b
      integer irow,icol,ires
      integer index_pair
      integer check_key
      logical uaadaa, ubbdbb
      logical uabdab,ubadba,ubadab,uabdba
      integer key_ix,nkey,first_h,second_h
      integer i1,j1,k1,l1
      integer inx1,inx2,inx3,inx4
#ifdef TCE_TIME_COMMS
      cpusecs(5) = cpusecs(5) - util_cpusec()
      cpusecs(55) = cpusecs(55) - util_wallsec()
#endif
      DO g3b = 1,noab+nvab
      DO g4b = g3b,noab+nvab
      DO g1b = 1,noab+nvab
      DO g2b = g1b,noab+nvab
      IF (int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1) .eq.
     &    int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)) THEN
      IF (ieor(int_mb(k_sym+g3b-1),ieor(int_mb(k_sym+g4b-1),
     &    ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))))
     &    .eq. irrep_v) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+g3b-1)+
     &    int_mb(k_spin+g4b-1)+int_mb(k_spin+g1b-1)+
     &    int_mb(k_spin+g2b-1).ne.8)) THEN
      nkey=g2b-1+(noab+nvab)*(g1b-1+(noab+nvab)*
     &     (g4b-1+(noab+nvab)*(g3b-1)))
      if(key_ix.eq.nkey) then
       inx3=g3b
       inx4=g4b
       inx1=g1b
       inx2=g2b
! first half
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
         go to 1000
      end if
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
 1000 CONTINUE
         irow = 0
         icol = 0
         ires = 0
         g1b  = 0
         g2b  = 0
         g3b  = 0
         g4b  = 0
      if (ioalg.eq.2) then
! v^{ g3b < g4b}_{ g1b < g2b } => ( g3b g1b | g4b g2b ) - ( g3b g2b | g4b g1b )
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
        if((g3b.ne.inx3).or.(g4b.ne.inx4).or.(g1b.ne.inx1).or.
     &     (g2b.ne.inx2)) then
           call errquit('get_block_ind_r: g-ind failed',0,UNKNOWN_ERR)
        end if
        ispin=int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1)+
     &        int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)
        if(ispin.eq.4) uaadaa=.true.
        if(ispin.eq.8) ubbdbb=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     &     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2)
     &     uabdab=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     &     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     &     ubadba=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     &     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2)
     &     ubadab=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     &     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     &     uabdba=.true.
        if (uaadaa.or.ubbdbb.or.uabdab.or.ubadba) then
         key_alpha = first_h
         if(key_alpha.ne.first_h) then
          write(6,*)'key_alpha first_h',key_alpha,first_h
          call util_flush(6)
          call errquit('get_block_ind_r: H1',0,UNKNOWN_ERR)
         end if
! defining the order
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
         if(model.eq.'ccsd_act') then
           call tce_hash_v2_act(int_mb(k_v2_alpha_offset),
     &                          key_alpha,off_a)
         else
           call tce_hash_v2(int_mb(k_v2_alpha_offset),
     &                      key_alpha,off_a)
         end if
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_r',2,MA_ERR)
         endif
         do i = 1, size
           f_a(i) = 0.0d0
         enddo
!$omp critical
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
!$omp end critical
         do i = 1, size
           array(i) = 0.0d0
         enddo
         size3 = int_mb(k_range+g3b-1)
         size1 = int_mb(k_range+g1b-1)
         size4 = int_mb(k_range+g4b-1)
         size2 = int_mb(k_range+g2b-1)
!
! --- ( g3 g1 | g4 g2 )
        if((.not.lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
         go to 200
        end if
! --- ( g3 g1 | g2 g4 )
        if((.not.lp31p42).and.(.not.l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
         go to 200
        end if
! --- ( g1 g3 | g4 g2 )
        if((.not.lp31p42).and.(l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
         go to 200
        end if
! --- ( g1 g3 | g2 g4 )
        if((.not.lp31p42).and.(l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
         go to 200
        end if
! --- ( g4 g2 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
         go to 200
        end if
! --- ( g4 g2 | g1 g3 )
        if((lp31p42).and.(l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
         go to 200
        end if
! --- ( g2 g4 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
         go to 200
        end if
! --- ( g2 g4 | g1 g3 )
        if((lp31p42).and.(l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
         go to 200
        end if
 200     continue
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_r',3,MA_ERR)
         endif
        end if !spin cases
        if (uaadaa.or.ubbdbb.or.uabdba.or.ubadab) then
         key_alpha = second_h
         if(key_alpha.ne.second_h) then
          write(6,*)'key_alpha second_h',key_alpha,second_h
          call util_flush(6)
          call errquit('get_block_ind_r: H2',0,UNKNOWN_ERR)
         end if
! defining the order
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
         if(uabdba.or.ubadab) then
          do i = 1, size
            array(i) = 0.0d0
          enddo
         end if
        if(model.eq.'ccsd_act') then
         call tce_hash_v2_act(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        else
         call tce_hash_v2(int_mb(k_v2_alpha_offset),key_alpha,off_a)
        end if
        allocate( f_a(1:size), stat=alloc_error)
        if (alloc_error.ne.0) then
          call errquit('get_block_ind_r',2,MA_ERR)
        endif
        do i = 1, size
          f_a(i) = 0.0d0
        enddo
!$omp critical
        call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
!$omp end critical
        size3 = int_mb(k_range+g3b-1)
        size1 = int_mb(k_range+g1b-1)
        size4 = int_mb(k_range+g4b-1)
        size2 = int_mb(k_range+g2b-1)
! --- ( g3 g2 | g4 g1 )
        if((.not.lp32p41).and.(.not.l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
         go to 300
        end if
! --- ( g3 g2 | g1 g4 )
        if((.not.lp32p41).and.(.not.l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
         go to 300
        end if
! --- ( g2 g3 | g4 g1 )
        if((.not.lp32p41).and.(l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
         go to 300
        end if
! --- ( g2 g3 | g1 g4 )
        if((.not.lp32p41).and.(l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
         go to 300
        end if
! --- ( g4 g1 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
         go to 300
        end if
! --- ( g4 g1 | g2 g3 )
        if((lp32p41).and.(l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
         go to 300
        end if
! --- ( g1 g4 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
         go to 300
        end if
! --- ( g1 g4 | g2 g3 )
        if((lp32p41).and.(l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
         go to 300
        end if
 300    continue
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_r',3,MA_ERR)
         endif
         end if !spin cases
      end if    !ioalg
#ifdef TCE_TIME_COMMS
      cpusecs(5) = cpusecs(5) + util_cpusec()
      cpusecs(55) = cpusecs(55) + util_wallsec()
#endif
      return
      end

      subroutine get_block_ind_i_r(d_file,array,size,key_ix,indexc,
     &                         w2b,w1b,w4b,w3b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "sym.fh"
      integer d_file
      integer size
      double precision array(size)
      double precision, allocatable :: f_a(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_a
#endif
      integer alloc_error
      integer i,j,k,l
      integer ispin
      integer indexc
      integer g1b,g2b,g3b,g4b
      integer w1b,w2b,w3b,w4b
      integer size1,size2,size3,size4
      integer off_a,key_alpha
      logical l31s,l42s,lp31p42
      logical l32s,l41s,lp32p41
      integer ig1b,ig2b,ig3b,ig4b
      integer irow,icol,ires
      integer index_pair
      integer check_key
      logical uaadaa, ubbdbb
      logical uabdab,ubadba,ubadab,uabdba
      integer key_ix,first_h,second_h
      integer i1,j1,k1,l1
      integer inx1,inx2,inx3,inx4
      external index_pair
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
         irow = 0
         icol = 0
         ires = 0
         g1b  = 0
         g2b  = 0
         g3b  = 0
         g4b  = 0
      if (ioalg.eq.2) then
! v^{ g3b < g4b}_{ g1b < g2b } => ( g3b g1b | g4b g2b ) - ( g3b g2b | g4b g1b )
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
! CHECK POINT
        if((g3b.ne.inx3).or.(g4b.ne.inx4).or.(g1b.ne.inx1).or.
     &     (g2b.ne.inx2)) then
           call errquit('get_block_ind_i_r: g-ind failed',0,UNKNOWN_ERR)
        end if
        ispin=int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1)+
     &        int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)
        if(ispin.eq.4) uaadaa=.true.
        if(ispin.eq.8) ubbdbb=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     &     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2)
     &     uabdab=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     &     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     &     ubadba=.true.
        if(int_mb(k_spin+g3b-1).eq.2.and.int_mb(k_spin+g4b-1).eq.1.and.
     &     int_mb(k_spin+g1b-1).eq.1.and.int_mb(k_spin+g2b-1).eq.2)
     &     ubadab=.true.
        if(int_mb(k_spin+g3b-1).eq.1.and.int_mb(k_spin+g4b-1).eq.2.and.
     &     int_mb(k_spin+g1b-1).eq.2.and.int_mb(k_spin+g2b-1).eq.1)
     &     uabdba=.true.
        if (uaadaa.or.ubbdbb.or.uabdab.or.ubadba) then
! first half
        key_alpha = first_h
! defining the order
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
         if(model.eq.'ccsd_act') then
           call tce_hash_v2_act(int_mb(k_v2_alpha_offset),
     &                          key_alpha,off_a)
         else
           call tce_hash_v2(int_mb(k_v2_alpha_offset),
     &                      key_alpha,off_a)
         end if
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_i_r',2,MA_ERR)
         endif
         do i = 1, size
           array(i) = 0.0d0
         enddo
!$omp critical
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
!$omp end critical
         size1 = int_mb(k_range+g1b-1)
         size2 = int_mb(k_range+g2b-1)
         size3 = int_mb(k_range+g3b-1)
         size4 = int_mb(k_range+g4b-1)
! --- ( g3 g1 | g4 g2 )
        if((.not.lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,2,3,1,1.0d0)
         go to 200
        end if
! --- ( g3 g1 | g2 g4 )
        if((.not.lp31p42).and.(.not.l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size1,size3,4,1,3,2,1.0d0)
         go to 200
        end if
! --- ( g1 g3 | g4 g2 )
        if((.not.lp31p42).and.(l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size4,size3,size1,3,2,4,1,1.0d0)
         go to 200
        end if
! --- ( g1 g3 | g2 g4 )
        if((.not.lp31p42).and.(l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size4,size2,size3,size1,3,1,4,2,1.0d0)
         go to 200
        end if
! --- ( g4 g2 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size2,size4,2,4,1,3,1.0d0)
         go to 200
        end if
! --- ( g4 g2 | g1 g3 )
        if((lp31p42).and.(l31s).and.(.not.l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size2,size4,1,4,2,3,1.0d0)
         go to 200
        end if
! --- ( g2 g4 | g3 g1 )
        if((lp31p42).and.(.not.l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size3,size4,size2,2,3,1,4,1.0d0)
         go to 200
        end if
! --- ( g2 g4 | g1 g3 )
        if((lp31p42).and.(l31s).and.(l42s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size1,size4,size2,1,3,2,4,1.0d0)
         go to 200
        end if
 200     continue
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_i_r',3,MA_ERR)
         endif
         end if !spin cases
        if (uaadaa.or.ubbdbb.or.uabdba.or.ubadab) then
! second half
        key_alpha = second_h
! defining the order
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
         if(uabdba.or.ubadab) then
          do i = 1, size
            array(i) = 0.0d0
          enddo
         end if
         if(model.eq.'ccsd_act') then
           call tce_hash_v2_act(int_mb(k_v2_alpha_offset),
     &                          key_alpha,off_a)
         else
           call tce_hash_v2(int_mb(k_v2_alpha_offset),
     &                      key_alpha,off_a)
         end if
         allocate( f_a(1:size), stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_i_r',2,MA_ERR)
         endif
!$omp critical
         call ga_get(d_v2orb,off_a+1,off_a+size,1,1,f_a,size)
!$omp end critical
         size3 = int_mb(k_range+g3b-1)
         size1 = int_mb(k_range+g1b-1)
         size4 = int_mb(k_range+g4b-1)
         size2 = int_mb(k_range+g2b-1)
! --- ( g3 g2 | g4 g1 )
        if((.not.lp32p41).and.(.not.l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size2,size3,4,2,1,3,-1.0d0)
         go to 300
        end if
! --- ( g3 g2 | g1 g4 )
        if((.not.lp32p41).and.(.not.l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,4,1,2,3,-1.0d0)
         go to 300
        end if
! --- ( g2 g3 | g4 g1 )
        if((.not.lp32p41).and.(l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,2,1,4,-1.0d0)
         go to 300
        end if
! --- ( g2 g3 | g1 g4 )
        if((.not.lp32p41).and.(l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size1,size4,size3,size2,3,1,2,4,-1.0d0)
         go to 300
        end if
! --- ( g4 g1 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size3,size1,size4,2,4,3,1,-1.0d0)
         go to 300
        end if
! --- ( g4 g1 | g2 g3 )
        if((lp32p41).and.(l32s).and.(.not.l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,4,3,2,-1.0d0)
         go to 300
        end if
! --- ( g1 g4 | g3 g2 )
        if((lp32p41).and.(.not.l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size2,size1,size4,size3,2,3,4,1,-1.0d0)
         go to 300
        end if
! --- ( g1 g4 | g2 g3 )
        if((lp32p41).and.(l32s).and.(l41s)) then
         call tce_sortacc_4(f_a,array,
     &                      size3,size2,size1,size4,1,3,4,2,-1.0d0)
         go to 300
        end if
 300    continue
         deallocate( f_a, stat=alloc_error)
         if (alloc_error.ne.0) then
           call errquit('get_block_ind_i_r',3,MA_ERR)
         endif
         end if !spin cases
      end if !ioalg
      return
      end


      subroutine get_block_r(d_file,array,size,offset)
      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer size
      integer offset
      double precision array(size)
#ifdef TCE_TIME_COMMS
      cpusecs(5) = cpusecs(5) - util_cpusec()
      cpusecs(55) = cpusecs(55) - util_wallsec()
#endif
      if (util_print('get_block_r',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     &  call errquit('get_block_r: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.2) then
!$omp critical
        call ga_get(d_file,offset+1,offset+size,1,1,array,size)
!$omp end critical
      else
        call errquit('get_block_r: ioalg=2 only!',ioalg,INPUT_ERR)
      endif
 9000 format(1x,'node',i3,' get_block_r request to file:',i10,
     &  ' size:',i10,' offset:',i10)
#ifdef TCE_TIME_COMMS
      cpusecs(5) = cpusecs(5) + util_cpusec()
      cpusecs(55) = cpusecs(55) + util_wallsec()
#endif
      return
      end
