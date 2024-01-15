      subroutine egrad_trace(ibas, aodens, ninfo, info, nbft, g_force,
     &  frc_nuc, frc_kin, frc_eff, frc_2el, frc_sow, frc_sox, frc_soy,
     &  frc_soz, frc_tot, nat, rtdb, nadct, aodensa, infoa, trans )
      implicit none
#include "errquit.fh"
* #include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "msgids.fh"
* #include "inp.fh"
* #include "apiP.fh"
* #include "cscfps.fh"
* #include "sym.fh"
* #include "cdft.fh"
#include "stdio.fh"
#include "cdft.fh"
    
      
      integer ibas
      integer aodens
      integer ninfo
      integer info(ninfo)  ! not to be confused with NWChem info
      integer nbft
      integer g_force, g_d1, g_eff
      integer g_sod(3)
      integer lforce
      double precision frc_nuc(3,nat)
      double precision frc_kin(3,nat)
      double precision frc_eff(3,nat)
      double precision frc_2el(3,nat)
      double precision frc_sow(3,nat)
      double precision frc_sox(3,nat)
      double precision frc_soy(3,nat)
      double precision frc_soz(3,nat)
      double precision frc_tot(3,nat)
      integer nat
      integer rtdb
      logical nadct

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
!      integer geom
!      logical oskel
      integer nshell
      integer ishell, mem1, max1e, lrecal, n1mx
      integer max2e, mscratch_1e, mscratch_2e, lbuf1, lscratch,
     &        max_at_bf, lsqatom
      integer l_buf, l_scr, l_ilab, l_jlab, l_sifbuf, l_d1, l_eff, 
     &        l_sifval, l_slab, l_sm, l_sym, l_map, l_b2s, 
     &        l_smin, l_smax, l_sodx, l_sody, l_sodz, l_buf2, l_scr2
      integer k_buf, k_scr, k_ilab, k_jlab, k_sifbuf, k_d1, k_eff,
     &        k_sifval, k_slab, k_sm, k_sym, k_map, k_b2s, 
     &        k_smin, k_smax, k_sodx, k_sody, k_sodz, k_buf2, k_scr2
      integer sizeof_int, szlabs, numd1, sz_bf_list
      integer i, j, ij, k
      integer me, nproc
      integer max2e4c, mscr2e4c, lbuf2, lscr2
c  
c     sifs parameters
c
      integer l1rec, n1max, l2rec, n2max, ierr 
      integer nsym, nmap
      parameter(nsym=1, nmap=0)
*     integer ietype(mxenrgy)
      integer kntin(nsym), nmpsy(nsym)
*     real*8 energy(nenrgy)
      integer itypea, itypeb
      integer numd1s, nbfs
      integer aodens2
      logical status
      INTEGER :: g_ad1, g_hd, l_asym, k_asym, l_hd, k_hd, nadx
      DOUBLE PRECISION :: nadxfl(nat,3)
      INTEGER :: hdol ! Handle or half-drv ovlp GA
      INTEGER :: aodensa, infoa(ninfo)
      LOGICAL :: trans
      INTEGER :: dims(3),chunk(3),iend(2),jstrt(3),jend(3) 
      INTEGER, PARAMETER :: istrt(2) = (/1,1/)

      me = ga_nodeid()
      nproc = ga_nnodes()
c
cgk debug
*      write(*,*)'gk: entered egrad_trace'
*      write(*,*)'gk: nbft = ', nbft
*      write(*,*)'gk: nadct = ',nadct
cgk end
c     
c     Get info about the basis sets
c     
cgk debug
*      write(*,*)'gk: ibas = ', ibas
cgk end
      if (.not. bas_numcont(ibas, nshell)) call errquit
     $     ('egrad_trace: bas_numcont failed for ibas', 
     &       ibas, BASIS_ERR)
cgk debug
*     write(*,*)'gk: nshell = ', nshell
cgk end

      l1rec = info(2)
      n1max = info(3)

c   allocate necessary local temporary arrays on the stack
      
c     get memory for sifs 1e labels
      sizeof_int=MA_sizeof(MT_INT, 1, MT_BYTE)
      if(sizeof_int.eq.4) then
        szlabs=(n1max+1)/2
        sz_bf_list=(nbft+1)/2
      elseif(sizeof_int.eq.8) then
        szlabs=n1max
        sz_bf_list=nbft
      else
        call errquit('wrt_dft_aoints: unsupporte integer type size=',
     &    sizeof_int, MA_ERR)
      endif

c     sifs integral label buffer (1 electron)
      if (.not. MA_push_get(MT_Int,2*szlabs,'wrt_dft_aoints:ilab',
     &        l_slab, k_slab))
     $     call errquit('wrt_dft_aoints: ma failed for sif labs', 
     &        2*szlabs, MA_ERR)

      lrecal=l1rec

c     get memory for sifs work buffer space
      if (.not. MA_push_get(MT_DBL, lrecal,'egrad_trace:sifbuf',
     &     l_sifbuf, k_sifbuf))
     $     call errquit('egrad_trace: ma failed getting sifbuf', lrecal,
     &        MA_ERR)

      n1mx=n1max

c     sifs value buffer
      if (.not. MA_push_get(MT_DBL, n1mx,'egrad_trace:sifval',
     &     l_sifval, k_sifval))
     $     call errquit('egrad_trace: ma failed getting sifval', n1mx, 
     &        MA_ERR)

      numd1=nbft*(nbft+1)/2
      numd1s=numd1
cgk debug
*      write(*,*)'gk: numd1=',numd1,' nbft=',nbft
cgk end

c     sifs total effective density
      if (.not. MA_push_get(MT_DBL, numd1,'egrad_trace:d1',
     &     l_sm, k_sm))
     $     call errquit('egrad_trace: ma failed getting d1', numd1, 
     &        MA_ERR)
c     memory for symmetry labels (all one)
      if (.not. MA_push_get(MT_Int, sz_bf_list,'egrad_trace:sym',
     &     l_sym, k_sym))
     $     call errquit('egrad_trace: ma failed getting sym',  
     &        sz_bf_list, MA_ERR)
c     memory for maps 
      if (.not. MA_push_get(MT_Int, sz_bf_list,'egrad_trace:map',
     &     l_map, k_map))
     $     call errquit('egrad_trace: ma failed getting map',  
     &        sz_bf_list, MA_ERR)
      
      nmpsy(1)=nbft
      nbfs=nbft

      if (.not. bas_geom(ibas, geom)) call errquit
     $   ('egrad_trace: basis corrupt?', 0, BASIS_ERR)

c     global 1e columbus matrices
      
c     effective total non spin orbit density
      g_d1 = ga_create_atom_blocked(geom, ibas, 
     &  'total eff density matrix')
c     effective non spin orbit fock matrix
      g_eff = ga_create_atom_blocked(geom, ibas,
     &  'effective fock matrix')
C LB
c     antisymmetric transition density
      g_ad1 = ga_create_atom_blocked(geom, ibas, 
     &  'antisymmteric density matrix')
c     half-derivative overlap, single derivative scratch
      g_hd = ga_create_atom_blocked(geom, ibas,
     &  'half-derivative overlap matrix')
c     half-derivative overlap, all derivatives (3*# atoms)
      dims(1) = 3*nat
      dims(2) = nbft
      dims(3) = nbft
      chunk(1)= dims(1)
      chunk(2)= -1
      chunk(3)= -1
      IF(.NOT. nga_create(MT_DBL,3,dims,
     &  'half-derivative overlap GA',chunk,hdol))     
     $  call errquit('egrad_trace: failed to create hdol GA',
     &       3*nat*nbft*nbft,GA_ERR)
      CALL ga_zero(hdol)
      CALL ga_print(hdol)
        
C LB
     
c   Do this for node 0 only to avoid conflics reading aodens
      if(me.eq.0) then
 
c       read one electron density
 
        itypea=0
        itypeb=7
cgk debug
!      write(*,*)'gk:  getting effective density'
!      write(*,*)'gk:  calling rd1mat from egrad_trace'
!      WRITE(*,*)"LB egrad_trace, l1rec=",l1rec,"n1max=",n1max,
!     &   "numd1s=",numd1s,"nsym=",nsym,"nmpsy=",nmpsy,"kntin=",kntin,
!     &   "nbfs=",nbfs
cgkd end
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 
cgk debug
!       write(*,*)'gk:  back from rd1mat'
!       write(*,*)'gk: total effective density matrix:'
!       ij=0
!       do i=1,nbft
!         do j=1,i
!           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
!           ij=ij+1
!         enddo
!       enddo
cgk end
c       put effective density in global arrays
        call sif2ga(g_d1,dbl_mb(k_sm),ibas,nbft)
 
c       get sifs effective fock matrix
        itypea=0
        itypeb=8
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 

cgk debug
!       write(*,*)'gk:  back from rd1mat'
!       write(*,*)'gk: total effective fock matrix:'
!       ij=0
!       do i=1,nbft
!         do j=1,i
!           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
!           ij=ij+1
!         enddo
!       enddo
cgk end

c       put effective density in global arrays
        call sif2ga(g_eff,dbl_mb(k_sm),ibas,nbft)

C LB

c       get antisymmetric density  matrix
!       USING MEMORY ALLOCATIONS FROM AODENS -- will this always work?
        IF(trans) THEN
          itypea=2
          itypeb=9
          call rd1mat(aodensa,infoa(2),infoa(3), info, dbl_mb(k_sifbuf),
     &      int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &      int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &      itypeb, nbfs) 

!        WRITE(*,*)"LB egrad_trace, back from rd1mat for asym dens"

c         put  antisymmetric density in global arrays
          call asif2ga(g_ad1,dbl_mb(k_sm),ibas,nbft)

*          write(*,*)'gk: antisymmteric density matrix global array'
          call ga_print(g_ad1)
        ENDIF

*        write(*,*)'gk: effective densitiy matrix global array'
        call ga_print(g_d1)
*        write(*,*)'gk: effective fock matrix global array'
        call ga_print(g_eff)
C LB
      end if
      call ga_sync()
c
c using ~/src/gradients/grad_force.F as guide.
c
      call schwarz_tidy()
      call intd_terminate()
      call int_init(rtdb, 1, ibas)
      call schwarz_init(geom, ibas)
      call int_terminate()
      call intd_init(rtdb, 1, ibas)

      max1e=0
      max2e=0
      mscratch_1e=0
      mscratch_2e=0
cgk debug
      call int_mem_1e(max1e,mscratch_1e)
*      write(*,*)'gk: max1e=',max1e
*      write(*,*)'gk: mscratch_1e=',mscratch_1e
cgk end
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
cgk debug
*      write(*,*)'gk: max1e=',max1e,' max2e=',max2e
*      write(*,*)'gk: mscratch_1e=',mscratch_1e,' mscratch_2e=',
*     & mscratch_2e
cgk end
*     max2e = max(max2e,1296*100)          ! 100 D quartets 
      lbuf1 = max(max1e, max2e)
      lbuf1 = 3*2*lbuf1
      lscratch = max(mscratch_1e,mscratch_2e)
cgk debug
*     write(*,*)'gk: in egrad_trace before meme alloc'
*     write(*,*)'gk:  max2e = ', max2e
*     write(*,*)'gk:  lbuf1 = ', lbuf1
*     write(*,*)'gk:  lscratch = ', lscratch
cgk end

c     buffers for one electron integral derivatives
      if (.not. ma_push_get(mt_dbl,lbuf1,'deriv buffer',l_buf,k_buf))
     $     call errquit('egrad_trace:could not allocate buffer',lbuf1,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()
cgk debug
*     write(*,*)'gk: l_scr=',l_scr,' k_scr=',k_scr,' lscratch=',
*    & lscratch
cgk end
      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch',
     $     l_scr, k_scr))
     &       call errquit('egrad_trace: scratch alloc failed',
     $         lscratch, MA_ERR)
      
c     local density matrix block 
      if (.not. bas_nbf_ce_max(ibas,max_at_bf))
     $     call errquit('egrad_trace: could not get max_at_bf',0, 
     $       BASIS_ERR)
      lsqatom = max_at_bf * max_at_bf
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_d1,k_d1))
     $     call errquit('egrad_trace:could not allocate l_d1',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_eff_fock',l_eff,
     $     k_eff)) call errquit
     &      ('egrad_trace;could not allocate l_eff',lsqatom, MA_ERR)

      oskel=.false.
c     g_force does not appear to be used
      call grad1(dbl_mb(k_buf), lbuf1, dbl_mb(k_scr), lscratch,
     $     dbl_mb(k_d1), dbl_mb(k_eff), frc_nuc,
     $     frc_kin, frc_eff, g_force,
     $     g_d1, g_eff, ibas, geom, nproc, nat, max_at_bf,
     $     rtdb, oskel)


      if(me.eq.0) then
        write(6,2200) 'kinetic energy gradient',' ',' ',
     $    ((frc_kin(i,j),i=1,3),j=1,nat)
        if(.not.nadct)
     $    write(6,2200) 'nuclear repulsion gradient',' ',' ',
     $      ((frc_nuc(i,j),i=1,3),j=1,nat)
        write(6,2200) 'weighted density gradient',' ',' ',
     $    ((frc_eff(i,j),i=1,3),j=1,nat)
 2200   format(A,A,A/,1000(3(1x,F16.10),/))
      endif

*     goto 1000

C LB  calculate half-derivative overlap matrices
      IF (trans) THEN
!        ENDIF
       if (.not. bas_numbf(AO_bas_han, nbf_ao)) then
           call errquit('Exiting from nadct_trace',1, BASIS_ERR)
        endif 
        CALL hdoverlap(dbl_mb(k_buf), lbuf1, lscratch, AO_bas_han, geom,
     &                  nproc, nat, max1e, mem1, hdol, nbft)

        nadxfl=0D0
        iend(1) =nbft
        iend(2) =nbft
        jstrt(2)=1
        jstrt(3)=1
        jend(2) =nbft
        jend(3) =nbft


        !CALL ga_print(hdol)
        DO i = 0, nat-1 !atom number
          DO j = 0, 2   !derivative direction
            jstrt(1)=i*3+j+1
            jend(1) =i*3+j+1
            !WRITE(*,*)"In egrad_trace, jstrt=",jstrt(1)
            CALL nga_copy_patch('N',hdol,jstrt,jend,g_hd,istrt,iend)
            CALL ga_print(g_hd)
            nadxfl(i+1,j+1)=-2*ga_ddot(g_hd,g_ad1)
          ENDDO
        ENDDO
        !CALL ga_print(hdol)
        WRITE(LUOUT,*)"nadxfl :"
        CALL output(nadxfl,1,nat,1,3,nat,3,1)
        OPEN(nadx,file='nadxfl',form='formatted')
        WRITE(nadx,1010) (nadxfl(i,:), i=1,nat)
        
1010    FORMAT (3d15.6)

      ENDIF
C LB

c     spin orbit contributions
      do i=1, 3
         status = ga_create(mt_dbl,nbft,nbft,'so density matrix', 0, 0,
     &        g_sod(i))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' egrad_trace: ', 'so density matrix'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('egrad_trace: ga_create ', 0, GA_ERR)
         endif
         call ga_zero(g_sod(i))
      end do

      call ga_zero(g_eff) ! reuse effective density matrix

      if (.not. ma_push_get(mt_dbl,lsqatom,'local sox den',l_sodx,
     $     k_sodx)) call errquit
     &      ('egrad_trace;could not allocate l_sodx',lsqatom, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local sox den',l_sody,
     $     k_sody)) call errquit
     &      ('egrad_trace;could not allocate l_sody',lsqatom, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local soz den',l_sodz,
     $     k_sodz)) call errquit
     &      ('egrad_trace;could not allocate l_sodz',lsqatom, MA_ERR)
      if(me.eq.0) then
 
        itypea=2
        itypeb=11
cgk debug
*     write(*,*)'gk:  x spin orbit density'
cgkd end
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 
cgk debug
*       write(*,*)'gk: x spin orbit density matrix matrix:'
*       ij=0
*       do i=1,nbft
*         do j=1,i
*           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
*           ij=ij+1
*         enddo
*       enddo
cgk end
c       put x so effective density in global arrays
        call asif2ga(g_sod(1),dbl_mb(k_sm),ibas,nbft)
        write(*,*)'x so effective density global array'
        call ga_print(g_sod(1))

        itypea=2
        itypeb=12
cgk debug
*     write(*,*)'gk:  y spin orbit density'
cgkd end
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 
cgk debug
*       write(*,*)'y spin orbit density matrix matrix:'
*       ij=0
*       do i=1,nbft
*         do j=1,i
*           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
*           ij=ij+1
*         enddo
*       enddo
cgk end
c       put y so effective density in global arrays
        call asif2ga(g_sod(2),dbl_mb(k_sm),ibas,nbft)
        write(*,*)'y so effective density global array'
        call ga_print(g_sod(2))

        itypea=2
        itypeb=13
cgk debug
*     write(*,*)'gk:  z spin orbit density'
cgkd end
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 
cgk debug
*       write(*,*)'gk: z spin orbit dinsity matrix matrix:'
*       ij=0
*       do i=1,nbft
*         do j=1,i
*           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
*           ij=ij+1
*         enddo
*       enddo
cgk end
c       put z so effective density in global arrays
        call asif2ga(g_sod(3),dbl_mb(k_sm),ibas,nbft)
        write(*,*)'z so effective density global array'
        call ga_print(g_sod(3))
 
c       get spin orbit effective fock matrix
        itypea=1
        itypeb=35
        call rd1mat(aodens, l1rec, n1max, info, dbl_mb(k_sifbuf), 
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_sm), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea, 
     &    itypeb, nbfs) 

cgk debug
*       write(*,*)'gk: total spin orbit effective fock matrix:'
*       ij=0
*       do i=1,nbft
*         do j=1,i
*           write(*,'(2i4,f20.12)') i,j,dbl_mb(k_sm+ij)
*           ij=ij+1
*         enddo
*       enddo
cgk end

c       put effective density in global arrays
        call sif2ga(g_eff,dbl_mb(k_sm),ibas,nbft)

cgk debug
*        write(*,*)'spin orbit effective fock matrix global array'
*        call ga_print(g_eff)
cgk end
      end if

cgk debug
*      write(*,*)'gk:calling sopgrdtrc'
cgk end
      call sopgrdtrc(dbl_mb(k_buf),lbuf1,dbl_mb(k_scr),lscratch,
     & dbl_mb(k_eff),frc_sow,g_sod,g_eff,ibas,geom,nproc,nat,max_at_bf,
     & oskel,frc_sox,frc_soy,frc_soz, dbl_mb(k_sodx),dbl_mb(k_sody),
     & dbl_mb(k_sodz))
cgk debug
*      write(*,*)'gk:  back from sopgrdtrc'
cgk end

      if(me.eq.0) then
        write(6,2201) 'SOx gradient',' ',' ',
     $    ((frc_sox(i,j),i=1,3),j=1,nat)
        write(6,2201) 'SOy gradient',' ',' ',
     $    ((frc_soy(i,j),i=1,3),j=1,nat)
        write(6,2201) 'SOz gradient',' ',' ',
     $    ((frc_soz(i,j),i=1,3),j=1,nat)
        write(6,2201) 'SO F.dS gradient',' ',' ',
     $    ((frc_sow(i,j),i=1,3),j=1,nat)
      endif

 2201   format(A,A,A/,1000(3(1x,E18.6),/))

      do j=1, nat
       do i=1,3
        frc_sow(i,j)=-frc_sow(i,j)+frc_sox(i,j)+frc_soy(i,j)+
     &                frc_soz(i,j)
       enddo
      enddo

      if(me.eq.0) then
        write(6,2201) 'total SO gradient',' ',' ',
     $    ((frc_sow(i,j),i=1,3),j=1,nat)
      endif

      
c     chop stack at first item allocated
      
      if (.not. MA_pop_stack(l_sodz)) call errquit
     $     ('egrad_trace: pop failed at l_sodz', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sody)) call errquit
     $     ('egrad_trace: pop failed at l_sody', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sodx)) call errquit
     $     ('egrad_trace: pop failed at l_sodx', 0, GA_ERR)

c     jump to here of no spin orbit

 1000 continue

      call schwarz_tidy()
      call intd_terminate()
      call int_init(rtdb, 1, ibas)
      call schwarz_init(geom, ibas)
      call int_terminate()
      call intd_init(rtdb, 1, ibas)

      if(me.eq.0) then
        call sifo2f(aodens, aodens, 'aodens', info, aodens2, ierr)
      endif

      l2rec = info(4)
      n2max = info(5)

c     memory for basis function to shell mapping array
      if (.not. ma_push_get(mt_dbl,nbft,'bf2shell_map',l_b2s,
     $     k_b2s)) call errquit
     &      ('egrad_trace;could not allocate l_b2s',nbft, MA_ERR)
      if (.not. ma_push_get(mt_dbl,nshell,'shell_min',l_smin,
     $     k_smin)) call errquit
     &      ('egrad_trace;could not allocate l_smin',nshell, MA_ERR)
      if (.not. ma_push_get(mt_dbl,nshell,'shell_min',l_smax,
     $     k_smax)) call errquit
     &      ('egrad_trace;could not allocate l_smax',nshell, MA_ERR)

      call int_mem_2e4c(max2e4c, mscr2e4c) ! blocking algorithm
*      write(*,*)'gk:  max2e4c = ', max2e4c
*      write(*,*)'gk:  mscr2e4c = ', mscr2e4c
      lbuf2=3*4*max2e4c
      lscr2=mscr2e4c

      if (.not. ma_push_get(mt_dbl,lbuf2,'buf2',l_buf2,
     $     k_buf2)) call errquit
     &      ('egrad_trace;could not allocate l_b2s',lbuf2, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lscr2,'buf2',l_scr2,
     $     k_scr2)) call errquit
     &      ('egrad_trace;could not allocate l_b2s',lscr2, MA_ERR)

c     two electron trace
      call sif2da_2e_sort(aodens, l2rec, n2max, info, ibas, nbft, 
     &  nshell, dbl_mb(k_b2s), dbl_mb(k_smin), dbl_mb(k_smax), frc_2el,
     &  nat, lbuf2, dbl_mb(k_buf2), lscr2, dbl_mb(k_scr2) )

      call schwarz_tidy()
      call intd_terminate()

      if(nadct) then
        do j=1, nat
         do i=1,3
          frc_tot(i,j)=frc_kin(i,j)+frc_eff(i,j)+
     &                 frc_sow(i,j)+frc_2el(i,j)
         enddo
        enddo
      else
        do j=1, nat
         do i=1,3
          frc_tot(i,j)=frc_kin(i,j)+frc_nuc(i,j)+frc_eff(i,j)+
     &                 frc_sow(i,j)+frc_2el(i,j)
         enddo
        enddo
      endif

      if(me.eq.0) then
        write(6,2200) '2 electron energy gradient',' ',' ',
     $    ((frc_2el(i,j),i=1,3),j=1,nat)
        write(6,2200) 'total gradient',' ',' ',
     $    ((frc_tot(i,j),i=1,3),j=1,nat)
        write(6,2201) 'sci total gradient',' ',' ',
     $    ((frc_tot(i,j),i=1,3),j=1,nat)
      endif

      if (.not. MA_pop_stack(l_scr2)) call errquit
     $     ('egrad_trace: pop failed at l_scr2', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf2)) call errquit
     $     ('egrad_trace: pop failed at l_buf2', 0, GA_ERR)
      if (.not. MA_pop_stack(l_smax)) call errquit
     $     ('egrad_trace: pop failed at l_smax', 0, GA_ERR)
      if (.not. MA_pop_stack(l_smin)) call errquit
     $     ('egrad_trace: pop failed at l_smin', 0, GA_ERR)
      if (.not. MA_pop_stack(l_b2s)) call errquit
     $     ('egrad_trace: pop failed at l_b2s', 0, GA_ERR)
      if (.not. MA_pop_stack(l_eff)) call errquit
     $     ('egrad_trace: pop failed at l_eff', 0, GA_ERR)
      if (.not. MA_pop_stack(l_d1)) call errquit
     $     ('egrad_trace: pop failed at l_d1', 0, GA_ERR)
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('egrad_trace: pop failed at l_scr', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('egrad_trace: pop failed at l_buf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_map)) call errquit
     $     ('egrad_trace: pop failed at l_map', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sym)) call errquit
     $     ('egrad_trace: pop failed at l_sym', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sm)) call errquit
     $     ('egrad_trace: pop failed at l_sm', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifval)) call errquit
     $     ('egrad_trace: pop failed at l_sifval', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifbuf)) call errquit
     $     ('egrad_trace: pop failed at l_sifbuf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_slab)) call errquit
     $     ('egrad_trace: pop failed at l_slab', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_bfnlab)) call errquit
*     if (.not. MA_pop_stack(l_bfnlab)) call errquit
*    $     ('egrad_trace: pop failed at l_bfnlab', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_info)) call errquit
*    $     ('egrad_trace: pop failed at l_info', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_jlab)) call errquit
*    $     ('egrad_trace: pop failed at l_jlab', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_ilab)) call errquit
*    $     ('egrad_trace: pop failed at l_ilab', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_scr)) call errquit
*    $     ('egrad_trace: pop failed at l_scr', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_buf)) call errquit
*    $     ('egrad_trace: pop failed at l_buf', 0, GA_ERR)
c
c     ----- restore ocache -----
c
      call ga_sync()            ! So that no nasty races can result
c     
cgk debug
*     write(*,*)'gk: leaving egrad_trace at end of routine'
cgk end
      end
