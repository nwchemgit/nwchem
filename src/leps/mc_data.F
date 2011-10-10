      subroutine mc_data_init(rtdb)
      implicit none
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
        if(.not. rtdb_get(rtdb,'mc_data:eigenvalues',mt_dbl,3*natom,
     &      dbl_mb(i_eval)))
     &      call errquit('mc_data:rtdb_get of eigenvalues failed',
     &      0, RTDB_ERR)

        if(.not. rtdb_get(rtdb,'mc_data:eigenvectors',mt_dbl,9*natom**2,
     &      dbl_mb(i_evec)))
     &      call errquit('mc_data:rtdb_get of eigenvectors failed',
     &      0, RTDB_ERR)

c        if(.not. rtdb_get(rtdb,'tst:coords',mt_dbl,3*natom,
c     &      dbl_mb(i_ctst)))
c     &      call errquit('mc_data:rtdb_get tst_coord',0, RTDB_ERR)

      end

      subroutine mc_data_geom_copy(rtdb,geom)
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb,geom


      if (.not. geom_masses_get(geom,natom,dbl_mb(i_m) ))
     & call errquit('mc_data_getgeom:  mass error',0,GEOM_ERR)

      if (.not. geom_cart_get(geom,natom, byte_mb(i_t),
     &                         dbl_mb(i_c),dbl_mb(i_q) ))
     & call errquit('mc_data_getgeom: coord error',0,GEOM_ERR)
c       write(*,*) 'mc_data geom_copy:Filled cartesians and masses'


      end

      function mc_data_get_i_ctst()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_ctst

      mc_data_get_i_ctst = i_ctst

      end

      subroutine mc_data_allocate(n,rtdb)
      implicit none
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer md_data_get_nsteps,rtdb
      integer i,n,nsteps
   
      call md_data_inp_set(rtdb)
      nsteps=md_data_get_nsteps() 
      write(*,*) 'In mc  data allocate:n,mcsteps,nsteps: ', n, 
     &  mcsteps,nsteps 
      if (.not. ma_alloc_get(mt_int,(mcsteps+1)*2*nsteps, 'sidecount',
     & h_sc,i_sc))
     & call errquit('mc_data: array of sidecount not allocated',
     & 1,MA_ERR)
       do i=0,(mcsteps+1)*2*nsteps-1
          int_mb(i_sc+i)=0
c          write(*,*) 'In mcdata i_sc,val=',i_sc+i,int_mb(i_sc+i)
       enddo
      if (.not. ma_alloc_get(mt_dbl, (mcsteps+1)*3*n, 'coordinates' ,
     & h_in,i_in))
     & call errquit('mc_data: array of coords not allocated',
     & 3*n,MA_ERR)

      if (.not. ma_alloc_get(mt_int, mcsteps+1, 'mcsteps' ,
     & h_s,i_s))
     & call errquit('mc_data: array of properties not allocated',
     & 1,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, n, 'mass' ,
     & h_m,i_m))
     & call errquit('mc_data: masses not allocated',
     & n,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, n, 'charge' ,
     & h_q,i_q))
     & call errquit('mc_data: charge not allocated',
     & n,MA_ERR)

      if (.not.ma_alloc_get(mt_byte, 16*n, 'tags', h_t,i_t))
     & call errquit('mc_data: tags not allocated', n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'force' ,
     & h_f,i_f))
     & call errquit('mc_data: force not allocated',
     & 3*n,MA_ERR)
c allocate memory for the eigenvalues and eigenvectors
      if (.not. ma_alloc_get(mt_dbl, 3*n, 'eigenvalues',
     & h_eval,i_eval))
     & call errquit('mc_data_geteigenvalue: not allocated',
     & n,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, 9*n*n, 'eigenvectors',
     & h_evec,i_evec))
     & call errquit('mc_data_geteigenvectors: not allocated',
     & n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'coords',
     & h_c,i_c))
     & call errquit('crossing tst allocate: not allocated',
     & n,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*n, 'tst_coords',
     & h_ctst,i_ctst))
     & call errquit('crossing tst allocate: not allocated',
     & n,MA_ERR)
c      write(*,*) 'In mc_data_allocatte: i_sc: ',i_sc
      end

      subroutine mc_data_allocate_vel_acc()
      implicit none
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

      if (.not. ma_alloc_get(mt_dbl, (naccept+1), 'property' ,
     & h_prp,i_prp))
     & call errquit('mc_data: array of property not allocated',
     & 0,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, (naccept+1), 'property_tmp' ,
     & h_prp_tmp,i_prp_tmp))
     & call errquit('mc_data: array of property_tmp not allocated',
     & 0,MA_ERR)
       
      if (.not. ma_alloc_get(mt_dbl, (naccept+1)*3*natom, 'velocity' ,
     & h_in_vel,i_in_vel))
     & call errquit('mc_data: array of init vels not allocated',
     & 3*natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, (naccept+1)*3*natom,'acceleration',
     & h_in_acc,i_in_acc))
     & call errquit('mc_data: array of init acc not allocated',
     & 3*natom,MA_ERR)
      end


      function mc_data_get_side_prev()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_side_prev

      mc_data_get_side_prev = side_prev

      end

      function mc_data_set_side_prev(sp)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_side_prev
      integer sp

      side_prev=sp
      mc_data_set_side_prev=.true. 
      end

      function mc_data_get_nxing()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_nxing

      mc_data_get_nxing = nxing

      end

      function mc_data_set_nxing(nx)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_nxing
      integer nx

      nxing=nx
      mc_data_set_nxing=.true.
      end

      function mc_data_get_trajnum()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_trajnum

      mc_data_get_trajnum = trajnum 

      end

      function mc_data_set_trajnum(tn)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_trajnum
      integer tn

      trajnum=tn
      mc_data_set_trajnum=.true.
      end

      function mc_data_get_i_prp()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_prp

      mc_data_get_i_prp = i_prp

      end

      function mc_data_get_i_prp_tmp()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_prp_tmp

      mc_data_get_i_prp_tmp = i_prp_tmp

      end

      function mc_data_get_i_sc()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_sc

      mc_data_get_i_sc = i_sc

      end

      function mc_data_get_i_m()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_m

      mc_data_get_i_m = i_m

      end

      function mc_data_get_i_f()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_f

      mc_data_get_i_f = i_f

      end

      function mc_data_get_i_in()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_in

      mc_data_get_i_in = i_in

      end

      function mc_data_get_i_s()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_s

      mc_data_get_i_s = i_s

      end

      function mc_data_get_i_c()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_c

      mc_data_get_i_c = i_c

      end

      function mc_data_get_i_t()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_t

      mc_data_get_i_t = i_t

      end

      function mc_data_get_i_q()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_q

      mc_data_get_i_q = i_q

      end

      function mc_data_get_i_in_vel()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_in_vel

      mc_data_get_i_in_vel = i_in_vel

      end

      function mc_data_get_i_in_acc()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_in_acc

      mc_data_get_i_in_acc = i_in_acc

      end


      subroutine mc_data_inp_set(rtdb)
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb
      if (.not.rtdb_get(rtdb,"mymc:temp",mt_dbl,1,temp))
     >      call errquit('failed to get temp',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymc:stepsize",mt_dbl,1,stepsize))
     >      call errquit('failed to get stepsize',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymc:mcsteps",mt_int,1,mcsteps))
     >      call errquit('failed to get mcsteps',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymc:ratio",mt_dbl,1,ratio))
     >      call errquit('failed to get ratio',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymc:mode",mt_int,1,mode))
     >      call errquit('failed to get mode',0,
     >       RTDB_ERR)
       if (.not.rtdb_get(rtdb,"mymc:seed",mt_dbl,1,seed))
     >      call errquit('failed to get seed',0,
     >       RTDB_ERR)
       if (.not.rtdb_get(rtdb,"mymc:gseed",mt_int,1,gseed))
     >      call errquit('failed to get gseed',0,
     >       RTDB_ERR)
      end

      subroutine mc_data_set(rtdb,geom)
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb
      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_c), dbl_mb(i_q)))
     & call errquit('mc_driver_getgeom: coord error',0,GEOM_ERR)
      if (.not. geom_masses_get(geom,natom,dbl_mb(i_m) ))
     & call errquit('mc_data_getgeom:  mass error',0,GEOM_ERR)

      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_in), dbl_mb(i_q)))
     & call errquit('mc_driver_getgeom: coord error',0,GEOM_ERR)

      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_ctst), dbl_mb(i_q)))
     & call errquit('mc_driver_getgeom: coord error',0,GEOM_ERR)

      end

      function mc_data_get_i_c_array()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_c_array
        mc_data_get_i_c_array = i_in    
      end 

      function mc_data_get_stepsize()
      implicit none
#include "mc_data.fh"
      double precision mc_data_get_stepsize

      mc_data_get_stepsize = stepsize

      end

      function mc_data_get_ratio()
      implicit none
#include "mc_data.fh"
      double precision mc_data_get_ratio

      mc_data_get_ratio = ratio

      end

      function mc_data_get_temp()
      implicit none
#include "mc_data.fh"
      double precision mc_data_get_temp

      mc_data_get_temp = temp

      end
      function mc_data_get_natom()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_natom

      mc_data_get_natom = natom

      end

      function mc_data_set_natom(n)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_natom 
      integer n
      
      natom=n
      mc_data_set_natom=.true. 
      end

      function mc_data_get_mcsteps()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_mcsteps

      mc_data_get_mcsteps = mcsteps

      end

      function mc_data_get_mode()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_mode

      mc_data_get_mode = mode

      end
      function mc_data_get_seed()
      implicit none
#include "mc_data.fh"
      double precision mc_data_get_seed
 
      mc_data_get_seed = seed 
 
      end
 
       function mc_data_get_gseed()
       implicit none
#include "mc_data.fh"
       integer mc_data_get_gseed
 
       mc_data_get_gseed = gseed
 
       end

      function mc_data_get_i_evec()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_evec

      mc_data_get_i_evec = i_evec

      end


      function mc_data_get_i_eval()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_i_eval

      mc_data_get_i_eval = i_eval

      end

      function mc_data_get_forside0()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_forside0

      mc_data_get_forside0 = forside0

      end

      function mc_data_get_forside()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_forside

      mc_data_get_forside = forside

      end

      function mc_data_get_backside0()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_backside0

      mc_data_get_backside0= backside0

      end

      function mc_data_get_backside()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_backside

      mc_data_get_backside= backside

      end

      function mc_data_set_forside0(f0)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_forside0
      integer f0
      forside0=f0
      mc_data_set_forside0=.true.
      end

      function mc_data_set_forside(f)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_forside
      integer f
      forside=f
      mc_data_set_forside=.true.
      end

      function mc_data_set_backside0(b0)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_backside0
      integer b0
      backside0=b0
      mc_data_set_backside0=.true. 
      end

      function mc_data_set_backside(b)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_backside
      integer b
      backside=b
      mc_data_set_backside=.true.
      end

      function mc_data_get_forward()
      implicit none
#include "mc_data.fh"
      logical mc_data_get_forward

      mc_data_get_forward = forward

      end

      function mc_data_get_backward()
      implicit none
#include "mc_data.fh"
      logical mc_data_get_backward

      mc_data_get_backward= backward

      end


      function mc_data_set_forward(f)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_forward
      logical f
      forward=f
      mc_data_set_forward=.true.
      end

      function mc_data_set_backward(b)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_backward
      logical b
      backward=b
      mc_data_set_backward=.true.
      end

      function mc_data_set_ifirc(irc)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_ifirc
      logical irc
      ifirc=irc
      mc_data_set_ifirc=.true.
      end

      function mc_data_get_ifirc()
      implicit none
#include "mc_data.fh"
      logical mc_data_get_ifirc
      mc_data_get_ifirc=ifirc

      end



      function mc_data_set_naccept(n)
      implicit none
#include "mc_data.fh"
      logical mc_data_set_naccept
      integer n 
      naccept=n
      mc_data_set_naccept=.true.
      end

      function mc_data_get_naccept()
      implicit none
#include "mc_data.fh"
      integer mc_data_get_naccept
      mc_data_get_naccept=naccept
      end

      subroutine mc_data_free_all()
      implicit none
#include "mc_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

      if(.not.ma_free_heap(h_c))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_q))
     &     call errquit('memory deallocation:charge',0,RTDB_ERR)
      if(.not.ma_free_heap(h_t))
     &     call errquit('memory deallocation:tag',0,RTDB_ERR)
      if(.not.ma_free_heap(h_m))
     &     call errquit('memory deallocation:mass',0,RTDB_ERR)
      if(.not.ma_free_heap(h_f))
     &     call errquit('memory deallocation:force',0,RTDB_ERR)
      if(.not.ma_free_heap(h_eval))
     &     call errquit('memory deallocation:eigenvalues',0,RTDB_ERR)
      if(.not.ma_free_heap(h_evec))
     &     call errquit('memory deallocation:eigenvectors',0,RTDB_ERR)
      if(.not.ma_free_heap(h_ctst))
     &     call errquit('memory deallocation:ctst',0,RTDB_ERR)
      if(.not.ma_free_heap(h_in))
     &     call errquit('memory deallocation:h_in',0,RTDB_ERR)
      if(.not.ma_free_heap(h_s))
     &     call errquit('memory deallocation:s',0,RTDB_ERR)
      if(.not.ma_free_heap(h_sc))
     &     call errquit('memory deallocation:sc',0,RTDB_ERR)
      if(.not.ma_free_heap(h_prp))
     &     call errquit('memory deallocation:prp',0,RTDB_ERR)
      if(.not.ma_free_heap(h_prp_tmp))
     &     call errquit('memory deallocation:prp_tmp',0,RTDB_ERR)
      if(.not.ma_free_heap(h_in_vel))
     &     call errquit('memory deallocation:in_vel',0,RTDB_ERR)
      if(.not.ma_free_heap(h_in_acc))
     &     call errquit('memory deallocation:in_acc',0,RTDB_ERR)

      end
c $Id$
