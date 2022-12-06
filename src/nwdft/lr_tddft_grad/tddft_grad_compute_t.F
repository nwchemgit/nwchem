c
c     Compute unrelaxed difference density matrix T
c     HvD 9/2006, NG 11/2012   
c
      subroutine tddft_grad_compute_t(tda,ipol,naoc,nav,nroot,
     +      g_xpy,g_xmy,g_tp)
c
c     Eq. 19: Furche & Ahlrich JCP (2002)
c     Tijs = -1/2 Sum_a [(X+Y)ias(X+Y)jas+(X-Y)ias(X-Y)jas]
c     Tabs =  1/2 Sum_i [(X+Y)ias(X+Y)ibs+(X-Y)ias(X-Y)ibs]
c     Tias = Tais = 0 (but we won't actually initialise these)
c
      implicit none
c
#include "global.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
c
c     Input:
c
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nroot     ! the number of states to consider
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
c     Output:
c
      integer g_tp(2)   ! global array handle for T
c
c     Local:
c
      integer alo(3)    ! lower chunk limits on A
      integer ahi(3)    ! upper chunk limits on A
      integer blo(3)    ! lower chunk limits on B
      integer bhi(3)    ! upper chunk limits on B
      integer clo(3)    ! lower chunk limits on C
      integer chi(3)    ! upper chunk limits on C
      integer ip        ! counter on polarization
      integer ir        ! counter on roots
      logical oroot
c
      character*32 pname
      pname="tddft_grad_compute_t: "
c
c     Loop over spins...
c
      do ip = 1, ipol
c
        call ga_zero(g_tp(ip))
c
c       Do the occupied-occupied block
c
        do ir = 1, nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = nav(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = nav(ip)
          blo(3) = 1
          bhi(3) = naoc(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = 1
          chi(2) = naoc(ip)
          clo(3) = 1
          chi(3) = naoc(ip)
c Daniel (12-3-12): Rather than calling nga_matmul_patch followed by
c nga_scale_patch, it is probably better to check whether the
c calculation is CIS first.
          if (tda) then
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)jas = Xjas
c (X-Y)ias = Xias , (X-Y)jas = Xjas
c Tij = - sum_a Xias*Xjas 
            call tddft_patch3mxm('n','t',-1d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          else
c (X+Y)ias(X+Y)jas 
            call tddft_patch3mxm('n','t',-0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
c (X-Y)ias(X-Y)jas
            call tddft_patch3mxm('n','t',-0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          endif ! tda
        enddo  ! ir = 1, nroot
c
c       Do the virtual-virtual block
c
        do ir = 1, nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = nav(ip)
          alo(3) = 1
          ahi(3) = naoc(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = naoc(ip)
          blo(3) = 1
          bhi(3) = nav(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = naoc(ip)+1
          chi(2) = naoc(ip)+nav(ip)
          clo(3) = naoc(ip)+1
          chi(3) = naoc(ip)+nav(ip)
c Daniel (12-3-12): Rather than calling nga_matmul_patch followed by
c nga_scale_patch, it is probably better to check whether the
c calculation is CIS first.
          if (tda) then
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)ibs = Xibs
c (X-Y)ias = Xias , (X-Y)ibs = Xibs
c Tab = sum_i Xias*Xibs 
            call tddft_patch3mxm('t','n',1d0,1.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          else
c (X+Y)ias(X+Y)ibs 
            call tddft_patch3mxm('t','n',0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
c (X-Y)ias(X-Y)ibs 
            call tddft_patch3mxm('t','n',0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          endif  ! tda
        enddo  ! ir = 1,nroot
c
      enddo ! ip = 1, ipol
c
      if (util_print('tddft grad t',print_debug)) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'T'
        call tddft_grad_print_array(ipol,nroot,g_tp,dble(ipol))
      endif
c
      end
c $Id$
      subroutine tddft_patch3mxm(transa,transb,alpha,beta,
     a     g_a, alo, ahi,
     b     g_b, blo, bhi,
     c     g_c, clo, chi)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      character*1 transa,transb
      integer	g_a, alo(3), ahi(3)
      integer	g_b, blo(3), bhi(3)
      integer	g_c, clo(3), chi(3)
      double precision	alpha, beta
c
      integer g_in1,g_in2
      g_in1=g_a
      g_in2=g_b
      if(transa.eq.'t'.or.transa.eq.'T') then
         call tddft_transpatch(g_a,g_in1,alo,ahi)
      endif
      if(transb.eq.'t'.or.transa.eq.'T') then
         call tddft_transpatch(g_b,g_in2,blo,bhi)
      endif
      call nga_matmul_patch('n','n',alpha,beta,
     +                            g_in1,alo,ahi,
     +                            g_in2,blo,bhi,
     +                            g_c,clo,chi)
      
      if(transa.eq.'t'.or.transa.eq.'T') then
         if (.not. ga_destroy(g_in1)) call errquit(
     e        ' tddft_patch3mxm_tn: ga_destroy',0,GA_ERR)
      endif
      if(transb.eq.'t'.or.transb.eq.'T') then
         if (.not. ga_destroy(g_in2)) call errquit(
     e        ' tddft_patch3mxm_tn: ga_destroy',0,GA_ERR)
      endif
     
      return
      end
      subroutine tddft_transpatch(g_a,g_t,alo,ahi)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer g_a
      integer g_t
      integer alo(3),ahi(3)
c      
      integer a_type,a_ndims,a_dims(3),t_dims(3)
      integer ichunk(3)
      integer mylo(3),myhi(3)
c      
      call  nga_inquire(g_a, a_type, a_ndims, a_dims)
      if(a_ndims.ne.3) call errquit(
     e     ' tddft_patch3mxm ready only for ndim=3',0,GA_ERR)
      
c     swap 2- and 3- for transpose
      t_dims(1)=a_dims(1)
      t_dims(2)=a_dims(3)
      t_dims(3)=a_dims(2)
      ichunk(1) = t_dims(1)
      ichunk(2) = -1
      ichunk(3) = t_dims(3)
      if (.not.nga_create(mt_dbl,a_ndims,t_dims,
     T     'transp g_b',ichunk,g_t)) call errquit(
     e     ' tddft_patch3mxm_tn: ga_destroy',0,GA_ERR)
#if 0      
      call nga_distribution(g_a,ga_nodeid(),mylo,myhi)
      write(6,'(i3,a,6i3)') ga_nodeid(), ' a distr ',mylo,myhi
      call nga_distribution(g_t,ga_nodeid(),mylo,myhi)
      write(6,'(i3,a,6i3)') ga_nodeid(), ' t distr ',mylo,myhi
      write(6,*) ' call copy_pa'
      write(6,'(i3,a,6i3)') ga_nodeid(), ' a ',alo,ahi

      write(6,'(i3,a,6i3)') ga_nodeid(), ' b ',blo,bhi
      write(6,'(i3,a,6i3)') ga_nodeid(), ' t ',alo,a_dims
#endif      
      call nga_copy_patch('t',
     I     g_a,alo,a_dims,
     O     g_t,alo,t_dims)
      return
      end
