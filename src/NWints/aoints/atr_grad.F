      subroutine atr_grad(ibas,aodensa,ninfo,info,nbft,nat,rtdb,frc_atr)
c     perform trace between 1 electron antisymmetric transition density
c     matrix and antisymmetric AO Fock matrix
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer ibas
      integer*4 aodensa
      integer*4 ninfo
      integer*4 info(ninfo)
      integer nbft
      integer nat
      integer rtdb
      double precision frc_atr(3,nat)

      integer*4 l1rec, n1max
      integer*4 nsym
      parameter(nsym=1) 
      integer*4 kntin(nsym), nmpsy(nsym)
      integer*4 numd1s, nbfs
      integer*4 itypea, itypeb

      integer geom
      integer g_ad1
      integer me, nproc
      integer sizeof_int, szlabs, sz_bf_list, lrecal, numd1
      integer l_slab, l_sifbuf, l_sifval, l_am, l_map, l_sym, l_buf,
     & l_scr, l_lad
      integer k_slab, k_sifbuf, k_sifval, k_am, k_map, k_sym, k_buf,
     & k_scr, k_lad
      integer n1mx
      integer i, j, ij
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer max1e, max2e, mscratch_1e, mscratch_2e, max_at_bf,
     & lspatom, lbuf, lscr
      
cgk debug
*     write(*,*)'gk: in atr_grad'
*     write(*,*)'gk;  nbft = ', nbft
cgk end

      me = ga_nodeid()
      nproc = ga_nnodes()
      l1rec = info(2)
      n1max = info(3)
cgk debug
*     write(*,*)'gk: l1rec = ', l1rec
*     write(*,*)'gk: n1max = ', n1max
cgk end

      call int_init(rtdb, 1, ibas)
      call schwarz_init(geom,ibas)
      call int_terminate()
      call intd_init(rtdb, 1, ibas)
    

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
cgk debug
*      write(*,*)'gk: sz_bf_list = ', sz_bf_list
*      write(*,*)'gk: nbft = ', nbft
cgk end

c     sifs integral label buffer (1 electron)
      if (.not. MA_push_get(MT_Int,2*szlabs,'wrt_dft_aoints:ilab',
     &        l_slab, k_slab))
     $     call errquit('wrt_dft_aoints: ma failed for sif labs',
     &        2*szlabs, MA_ERR)

      lrecal=l1rec

c     get memory for sifs work buffer space
      if (.not. MA_push_get(MT_DBL, lrecal,'atr_trace:sifbuf',
     &     l_sifbuf, k_sifbuf))
     $     call errquit('atr_trace: ma failed getting sifbuf', lrecal,
     &        MA_ERR)

      n1mx=n1max

c     sifs value buffer
      if (.not. MA_push_get(MT_DBL, n1mx,'atr_trace:sifval',
     &     l_sifval, k_sifval))
     $     call errquit('atr_trace: ma failed getting sifval', n1mx,
     &        MA_ERR)

      numd1=nbft*(nbft+1)/2
      numd1s=numd1

cgk debug
*     write(*,*)'gk: numd1 = ', numd1
cgk end

c     memory for symmetry labels (all one)
      if (.not. MA_push_get(MT_Int, sz_bf_list,'atr_grad:sym',
     &     l_sym, k_sym))
     $     call errquit('atr_grad: ma failed getting sym',
     &        sz_bf_list, MA_ERR)
c     memory for maps
      if (.not. MA_push_get(MT_Int, sz_bf_list,'atr_grad:map',
     &     l_map, k_map))
     $     call errquit('atr_grad: ma failed getting map',
     &        sz_bf_list, MA_ERR)
c     sifs total effective density
      if (.not. MA_push_get(MT_DBL, numd1,'atr_grad:d1',
     &     l_am, k_am))
     $     call errquit('atr_grad: ma failed getting d1', numd1,
     &        MA_ERR)

      nmpsy(1)=nbft
      nbfs=nbft

      if (.not. bas_geom(ibas, geom)) call errquit
     $   ('atr_grad: basis corrupt?', 0, BASIS_ERR)

c     global 1e columbus matrices

c     effective total non spin orbit density
      g_ad1 = ga_create_atom_blocked(geom, ibas,
     &  'total eff density matrix')

      if (me.eq.0) then
c       read one electron density

        itypea=2
        itypeb=9
cgk debug
*     write(*,*)'gk:  getting effective density'
*     write(*,*)'gk:  calling rd1mat from atr_grad'
cgkd end
        call rd1mat(aodensa, l1rec, n1max, info, dbl_mb(k_sifbuf),
     &    int_mb(k_slab), dbl_mb(k_sifval), dbl_mb(k_am), numd1s,
     &    int_mb(k_sym), int_mb(k_map), nsym, nmpsy, kntin, itypea,
     &    itypeb, nbfs)
cgk debug
        write(*,*)'gk:  back from rd1mat'
        write(*,*)'gk: total antisymmetric transition density matrix:'
        ij=0
        do i=1,nbft
          do j=1,i
            write(*,'(2i4,f20.12)') i,j,dbl_mb(k_am+ij)
            ij=ij+1
          enddo
        enddo
cgk end
c       put effective density in global arrays
        call asif2ga(g_ad1,dbl_mb(k_am),ibas,nbft)
      endif

cgk debug
*     write(*,*)'gk: antisymmetric transistion density global array'
*     call ga_print(g_ad1)
cgk end

c     we no longer need the sifs matrix
      if (.not. MA_pop_stack(l_am)) call errquit
     $     ('atr_grad: pop failed at l_am', 0, GA_ERR)

      call int_mem_1e(max1e, mscratch_1e)
cgk debug
*     write(*,*)'gk: max1e=',max1e
*     write(*,*)'gk: mscratch_1e=',mscratch_1e
cgk end
      lbuf = 3*nat*max1e
      lscr = 3*nat*mscratch_1e

c     buffers for one electron integral derivatives
      if (.not. ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit('atr_grad:could not allocate buffer',lbuf,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()
      if (.not. ma_push_get(mt_dbl,lscr,'deriv scratch',
     $     l_scr, k_scr))
     &       call errquit('atr_grad: scratch alloc failed',
     $         lscr, MA_ERR)
c     local density matrix block 
      if (.not. bas_nbf_ce_max(ibas,max_at_bf))
     $     call errquit('atr_grad: could not get max_at_bf',0,
     $       BASIS_ERR)
      lspatom = max_at_bf * max_at_bf
      if (.not. ma_push_get(mt_dbl,lspatom,'local_density',
     $     l_lad,k_lad))
     &       call errquit('atr_grad: local_density alloc failed',
     $         lspatom, MA_ERR)

   
      call adaovl(dbl_mb(k_buf),lbuf,dbl_mb(k_scr),lscr,dbl_mb(k_lad),
     &  frc_atr, g_ad1, ibas, geom, nproc, nat, max_at_bf,.false.)

      call schwarz_tidy()
      call intd_terminate()


      if (.not. MA_pop_stack(l_lad)) call errquit
     $     ('atr_grad: pop failed at l_lad', 0, GA_ERR)
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('atr_grad: pop failed at l_scr', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('atr_grad: pop failed at l_buf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_map)) call errquit
     $     ('atr_grad: pop failed at l_map', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sym)) call errquit
     $     ('atr_grad: pop failed at l_sym', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifval)) call errquit
     $     ('atr_grad: pop failed at l_sifval', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifbuf)) call errquit
     $     ('atr_grad: pop failed at l_sifbuf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_slab)) call errquit
     $     ('atr_grad: pop failed at l_slab', 0, GA_ERR)

      return
      end
