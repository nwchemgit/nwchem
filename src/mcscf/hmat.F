c     
c     =========================================================
c     Set of routines to generate CASSCF orbital Hessian
c     explicitly in core. For debugging.
c     =========================================================
c     
c     
      subroutine hmat( nbf, nclosed, nact, orlen, dm1, dm2,
     $                 g_ifock, g_afock, g_gfock,
     $                 g_coul, g_exch, h )
      implicit none
      integer nbf, nclosed, nact, orlen
      double precision dm1(nact,nact)
      double precision dm2(nact,nact,nact,nact)
      integer g_coul, g_exch
      integer g_ifock, g_afock, g_gfock
      double precision h(orlen,orlen)
      double precision xx(1000),yy(1000),zz(1000)
      double precision hlp(1000000)
      integer i,j,nocc,nnocc
      logical odofock, odoao, odomo
      data odofock, odoao, odomo/.false.,.true.,.true./
c     
      call dfill((orlen*orlen),0.d0,h,1)
      nocc = nclosed + nact
      nnocc = (nocc*(nocc+1))/2
c
c  1e fock section:
c
      if (odofock) then
        call ga_get(g_ifock,1,nbf,1,nbf,xx,nbf)
        call ga_get(g_afock,1,nbf,1,nbf,yy,nbf)
        call ga_get(g_gfock,1,nbf,1,nbf,zz,nbf)
        call hmat_fock( nbf, nclosed, nact, orlen, dm1,
     $                  xx, yy, zz, h )
      endif
c     
c  AO section: P and K supermatrix contracted against 1-pdm
c     

      if (odoao) then
        call getaoints( nbf, nclosed, nact, nnocc, hlp )
        call hmat_ao_iaia( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
        call hmat_ao_iaiv( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
        call hmat_ao_iaav( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
        call hmat_ao_iviv( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
        call hmat_ao_ivav( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
      endif
c     
c  MO section: MO integrals contracted against 2-pdm
c     
      if (odomo) then
        call hmat_mo_iaia( nbf, nclosed, nact, orlen, dm2,
     $                     g_coul, g_exch, h, xx, yy )
        call hmat_mo_iaav( nbf, nclosed, nact, orlen, dm2,
     $                     g_coul, g_exch, h, xx, yy )
        call hmat_mo_avav( nbf, nclosed, nact, orlen, dm2,
     $                     g_coul, g_exch, h, xx, yy )
      endif
c     
c  Scale
c     
      call dscal((orlen*orlen),2.d0,h,1)
c     
c  Write to file
c     
      open(unit=55,file='hessX',form='unformatted',
     $     status='unknown')
      do i=1,orlen
        write(55) (h(j,i),j=1,orlen)
      enddo
      close(55)
c$$$      print*,'Exact Hessian diagonal from HMAT'
c$$$      write(6,222) (h(i,i),i=1,orlen)
c$$$ 222  format(10f12.6)
      
      
      return
      end
      





      subroutine hmat_mo_iaia( nbf, nclosed, nact, orlen, dm2,
     $                      g_coul, g_exch, h, xj, xk )
      implicit none
      integer nbf, nclosed, nact, orlen
      double precision dm2(nact,nact,nact,nact)
      integer g_coul, g_exch
      double precision h(orlen,orlen)
      double precision xj(nbf,nbf),xk(nbf,nbf)
c
      integer nvir
      integer i,j,t,u,v,x,vx,it,ju,nn,iaoff
      double precision xx,yy
c
      nvir = nbf - nclosed - nact
      nn = nbf*nbf
      iaoff = (nclosed+nact)*nvir
      do v=1,nact
        do x=1,nact
          vx = (max(v,x)*(max(v,x)-1))/2+min(v,x)
          call ga_get(g_coul,1,nn,vx,vx,xj,1)
          call ga_get(g_exch,1,nn,vx,vx,xk,1)
          if (x.gt.v) then
            call mattrp(nbf,xj)
            call mattrp(nbf,xk)
          endif
          do i=1,nclosed
            do t=1,nact
              it = iaoff + (i-1)*nact + t
              do j=1,nclosed
                do u=1,nact
                  ju = iaoff + (j-1)*nact + u
                  xx = 2.d0*xj(j,i)*dm2(u,t,v,x)
                  yy = 2.d0*xk(j,i)*(dm2(u,x,v,t)+dm2(u,x,t,v))
                  h(it,ju) = h(it,ju) + xx + yy
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c$$$      call gmatprint(orlen,orlen,iaoff+1,orlen,iaoff+1,orlen,h)
c$$$      write(6,*)
      return
      end




      subroutine hmat_mo_iaav( nbf, nclosed, nact, orlen, dm2,
     $                      g_coul, g_exch, h, xj, xk )
      implicit none
      integer nbf, nclosed, nact, orlen
      double precision dm2(nact,nact,nact,nact)
      integer g_coul, g_exch
      double precision h(orlen,orlen)
      double precision xj(nbf,nbf),xk(nbf,nbf)
c
      integer nvir
      integer i,t,u,a,v,x,vx,it,ua,nn,iaoff,avoff,voff,aa
      double precision xx,yy
c
      nvir = nbf - nclosed - nact
      nn = nbf*nbf
      iaoff = (nclosed+nact)*nvir
      avoff = nclosed*nvir
      voff = nclosed+nact
      do v=1,nact
        do x=1,nact
          vx = (max(v,x)*(max(v,x)-1))/2+min(v,x)
          call ga_get(g_coul,1,nn,vx,vx,xj,1)
          call ga_get(g_exch,1,nn,vx,vx,xk,1)
          if (x.gt.v) then
            call mattrp(nbf,xj)
            call mattrp(nbf,xk)
          endif
          do i=1,nclosed
            do t=1,nact
              it = iaoff + (i-1)*nact + t
              do u=1,nact
                do a=1,nvir
                  ua = avoff + (u-1)*nvir + a
                  aa = voff + a
                  xx = -2.d0*xj(i,aa)*dm2(t,u,v,x)
                  yy = -2.d0*xk(aa,i)*(dm2(t,v,u,x)+dm2(t,v,x,u))
                  h(it,ua) = h(it,ua) + xx + yy
                  h(ua,it) = h(ua,it) + xx + yy
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c$$$      call gmatprint(orlen,orlen,iaoff+1,orlen,avoff+1,iaoff,h)
c$$$      write(6,*)
c$$$      call gmatprint(orlen,orlen,avoff+1,iaoff,iaoff+1,orlen,h)
c$$$      write(6,*)
      return
      end






      subroutine hmat_mo_avav( nbf, nclosed, nact, orlen, dm2,
     $                      g_coul, g_exch, h, xj, xk )
      implicit none
      integer nbf, nclosed, nact, orlen
      double precision dm2(nact,nact,nact,nact)
      integer g_coul, g_exch
      double precision h(orlen,orlen)
      double precision xj(nbf,nbf),xk(nbf,nbf)
c
      integer nvir
      integer t,u,a,b,v,x,vx,ta,ub,nn,iaoff,avoff,voff,aa,bb
      double precision xx,yy
c
      nvir = nbf - nclosed - nact
      nn = nbf*nbf
      iaoff = (nclosed+nact)*nvir
      avoff = nclosed*nvir
      voff = nclosed+nact
      do v=1,nact
        do x=1,nact
          vx = (max(v,x)*(max(v,x)-1))/2+min(v,x)
          call ga_get(g_coul,1,nn,vx,vx,xj,1)
          call ga_get(g_exch,1,nn,vx,vx,xk,1)
          if (x.gt.v) then
            call mattrp(nbf,xj)
            call mattrp(nbf,xk)
          endif
          do t=1,nact
            do a=1,nvir
              ta = avoff + (t-1)*nvir + a
              aa = voff + a
              do u=1,nact
                do b=1,nvir
                  ub = avoff + (u-1)*nvir + b
                  bb = voff + b
                  xx = 2.d0*xj(bb,aa)*dm2(t,u,v,x)
                  yy = 2.d0*xk(aa,bb)*(dm2(t,x,v,u)+dm2(t,x,u,v))
                  h(ta,ub) = h(ta,ub) + xx + yy
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c$$$      call gmatprint(orlen,orlen,avoff+1,iaoff,avoff+1,iaoff,h)
c$$$      write(6,*)
      return
      end





      

      subroutine mattrp( n, x )
      integer n
      double precision x(n,n)
      double precision y(1000)
      integer i,j,ij

      do i=1,n
        do j=1,n
          ij =  (i-1)*n + j
          y(ij) = x(i,j)
        enddo
      enddo
      call dcopy((n*n),y,1,x,1)
      return
      end

        
      



      subroutine gmatprint(n,m,ilo,ihi,jlo,jhi,x)
      integer n,m,ilo,ihi,jlo,jhi
      double precision x(n,m)
      integer i,j,jtop

      jtop = min((jlo+9),jhi)
      do i=ilo,ihi
        write(6,881) (x(i,j),j=jlo,jtop)
 881    format(10f12.5)
      enddo
      return
      end











       subroutine mcscf_0hessv( geom, basis, nbf, nclosed, nact,
     $                         oskel, pflg, lshift, dm1, dm2,
     $                         g_movecs, g_ifock, g_afock, g_gfock,
     $                         g_coul, g_exch, g_x, g_ax )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       logical oskel                                ! [input] Symmetry toggle
       integer pflg                                 ! [input] Hessian-type
       double precision lshift                      ! [input] Level-shift
       double precision dm1(nact,nact)              ! [input] 1PDM
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_ifock                              ! [input] Inactive Fock matrix
       integer g_afock                              ! [input] Active Fock matrix
       integer g_gfock                              ! [input] Generalized Fock matrix
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
c
c
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
c
c
c
       if (.not.ga_duplicate(g_afock, g_tmp1, 'temp1'))
     $      call errquit('mcscf_hessv: cannot duplicate',0)
       if (.not.ga_duplicate(g_afock, g_tmp2, 'temp2'))
     $      call errquit('mcscf_hessv: cannot duplicate',0)
       if (.not.ga_duplicate(g_afock, g_tmp3, 'temp3'))
     $      call errquit('mcscf_hessv: cannot duplicate',0)
c
c
c
       call ga_zero(g_ax)
       call mcscf_hessv_1e00( geom, basis, nbf, nclosed, nact,
     $                        lshift, dm1, g_ifock, g_afock, g_gfock,
     $                        g_tmp1, g_tmp2, g_tmp3,
     $                        g_x, g_ax )

       call mcscf_hessv_2e00( geom, basis, nbf, nclosed, nact,
     $                      oskel,
     $                      dm1, dm2, g_movecs, g_coul, g_exch, 
     $                      g_tmp1, g_tmp2, g_tmp3, g_x, g_ax )
c
c
c

       if (.not.ga_destroy(g_tmp1))
     $      call errquit('mcscf_hessv: cannot destroy temp1',0)
       if (.not.ga_destroy(g_tmp2))
     $      call errquit('mcscf_hessv: cannot destroy temp2',0)
       if (.not.ga_destroy(g_tmp3))
     $      call errquit('mcscf_hessv: cannot destroy temp3',0)
c
c
c
       return
       end









c
c Return the MCSCF 1e Hessian-vector product
c
c NB: some analogous 1e-terms from ROHF expression
c     cannot be computed from the Fock matrices
c     and are actually included in the 2e component
c

       subroutine mcscf_hessv_1e00( geom, basis, nbf, nclosed, nact,
     $                            lshift, dm1, g_ifock, g_afock,
     $                            g_gfock, g_tmp1, g_tmp2, g_tmp3,
     $                            g_x, g_ax )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision lshift                      ! [input] Level shift
       double precision dm1(nact,nact)              ! [input] 1PDM
       integer g_ifock                              ! [input] Inactive Fock matrix
       integer g_afock                              ! [input] Active Fock matrix
       integer g_gfock                              ! [input] Generalized Fock matrix
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
c
       integer g_tmp4
       integer nvir, vlen, voff, aoff, aend
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ga_duplicate(g_tmp3, g_tmp4, 'temp4'))
     $      call errquit('mcscf_hessv_1e: cannot duplicate',0)
       call ga_zero(g_tmp4)
       call ga_zero(g_tmp3)
c
c Reshape argument vector by copying into matrix
c
       call ga_zero(g_tmp1)
       call mcscf_vec2mat( nbf, nclosed, nact, 1.d0, g_x, 0.d0, g_tmp1)
c
c 1-PDM density matrix is hidden in unused 
c active-active subblock
c
       call ga_put(g_tmp1, aoff, aend, aoff, aend, dm1, nact )
c
c                                          cv
c Useful intermediate (equivalent to ROHF F  )
c
c    Z     I     A
c   F  =  F  +  F
c
c
       call ga_dadd(1.d0, g_ifock, 1.d0, g_afock, g_tmp2 )
c
c Closed-virtual bit
c
       call ga_matmul_patch( 'n', 'n', 4.d0, 0.d0,
     $                       g_tmp2, voff, nbf, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
       call ga_matmul_patch( 'n', 'n', -4.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp2, 1, nclosed, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
       call ga_matmul_patch( 'n', 'n', -1.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_tmp1, voff, nbf, aoff, aend,     ! Careful here ...        
     $                       g_gfock, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
       call ga_matmul_patch( 'n', 'n', 4.d0, 1.d0,
     $                       g_tmp2, voff, nbf, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
       call ga_matmul_patch( 't', 'n', -1.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_gfock, voff, nbf, aoff, aend,    ! Careful here ...        
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c
c Active-virtual bit
c       
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp2, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 't', -1.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_gfock, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,     
     $                       g_tmp1, aoff, aend, aoff, aend,    !                                     
     $                       g_tmp4, voff, nbf, aoff, aend )    ! X' = X.D    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,              !                                     
     $                       g_ifock, voff, nbf, voff, nbf,     ! Z = F.X' = F.X.D                    
     $                       g_tmp4, voff, nbf, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_tmp1, voff, nbf, aoff, aend,     ! Careful here ...        
     $                       g_gfock, aoff, aend, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_tmp1, 1, nclosed, aoff, aend,
     $                       g_tmp1, aoff, aend, aoff, aend,    !                                         
     $                       g_tmp4, 1, nclosed, aoff, aend )   ! X' = tr(X).D    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,             !                                         
     $                       g_ifock, voff, nbf, 1, nclosed,    ! Z = F.X' = F.tr(X).D                    
     $                       g_tmp4, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 't', 2.d0, 1.d0,
     $                       g_tmp2, voff, nbf, 1, nclosed,
     $                       g_tmp1, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c
c Closed-active bit
c
c      ia      iv         av
c    Y1   = 4.X  (iF + aF)                    (eqn A2)
c
       call ga_matmul_patch( 'n', 'n', 4.d0, 1.d0,
     $                       g_tmp2, aoff, aend, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia    iv av
c    Y2   = X  F                              (eqn A2)
c
       call ga_matmul_patch( 'n', 'n', -1.d0, 1.d0,
     $                       g_gfock, aoff, aend, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia         av   iv
c    Y1   = -2.d.X  .iF                       (eqn A3)
c               
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, voff, nbf,    !                                         
     $                       g_tmp4, aoff, aend, voff, nbf )   ! X' = D.tr(X)    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,            !                                         
     $                       g_tmp4, aoff, aend, voff, nbf,    ! Z = X'.F = D.tr(X).F                  
     $                       g_ifock, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia      av          iv
c    Y2   = 2.X  .(iF + aF)                  (eqn A3)
c
       call ga_matmul_patch( 't', 'n', 2.d0, 1.d0,
     $                       g_tmp1, aoff, aend, voff, nbf,
     $                       g_tmp2, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia        ia   ii
c    Y3   = 2.d.X  .iF                       (eqn A1)
c
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,   !                                       
     $                       g_tmp4, aoff, aend, 1, nclosed )  ! X' = D.X    (multiply arg by 1PDM)
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,             !                                       
     $                       g_tmp4, aoff, aend, 1, nclosed,   ! Z = X'.F = D.X.F                  
     $                       g_ifock, 1, nclosed, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia      ia               aa
c    Y4   = 2.X  (2iF + 2aF - F)             (eqn A1)
c             
                                                               !  x       z         Should really use patch routine
       call ga_dadd(2.d0, g_tmp2, -1.d0, g_gfock, g_tmp4 )     ! F  =  2.F  - F     here but only active-active
                                                               !                    block of g_tmp4 used.
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,
     $                       g_tmp4, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia       ia         ii
c    Y5   = -4.X  (iF + aF)                  (eqn A1)
c
       call ga_matmul_patch( 'n', 'n', -4.d0, 1.d0,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp2, 1, nclosed, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c Copy result from matrix format to vector format
c
       call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
       call mcscf_mat2vec( nbf, nclosed, nact, 1.d0, g_tmp3, 1.d0, g_ax)
c
c
c
       if (.not.ga_destroy(g_tmp4))
     $      call errquit('mcscf_hessv_1e: cannot destroy',0)
       return
       end

       















       subroutine mcscf_hessv_2e00( geom, basis, nbf, nclosed, nact,
     $                            oskel, dm1, dm2, g_movecs, g_coul,
     $                            g_exch, g_tmp1, g_tmp2, g_tmp3,
     $                            g_x, g_ax )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       logical oskel                                ! [input] Symmetry toggle
       double precision dm1(nact,nact)              ! [input] 1PDM
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
       integer nvir, vlen, voff, aoff, aend
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
c
       call ga_zero( g_tmp3 )
       call ga_zero( g_tmp1 )
       call mcscf_vec2mat( nbf, nclosed, nact, 1.d0, g_x, 0.d0, g_tmp1 )
       call mcscf_hessv_2e_ao00( geom, basis, nbf, nclosed, nact, oskel, 
     $                         dm1, g_movecs, g_tmp1, g_tmp2, g_tmp3 )
c
c
c
       call mcscf_hessv_2e_mo00( nbf, nclosed, nact, dm2,
     $                         g_coul, g_exch, g_tmp1, g_tmp3 )
c
c
c
       call mcscf_mat2vec(nbf, nclosed, nact, 1.d0, g_tmp3, 1.d0, g_ax )
c
c
c
       return
       end
















       subroutine mcscf_hessv_2e_mo00( nbf, nclosed, nact,
     $                               dm2, g_coul, g_exch,
     $                               g_x, g_ax )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_x                                  ! [input] Argument matrix
       integer g_ax                                 ! [output] Hessian product (in matrix format)
c
c
       integer nvir, vlen, voff, aoff, aend
       integer l_xit, k_xit, l_xua, k_xua, l_xju, k_xju
       logical status
c
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact

c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'xit',l_xit,k_xit))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0)
       call ga_get(g_x,aoff,aend,1,nclosed,dbl_mb(k_xit),nact)
       call mcscf_hessv_2e_mo_avia( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xit), g_ax )
       if (.not.ma_pop_stack(l_xit))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0)
c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nvir),'xua',l_xua,k_xua))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0)
       call ga_get(g_x,voff,nbf,aoff,aend,dbl_mb(k_xua),nvir)
       call mcscf_hessv_2e_mo_avav( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xua), g_ax )
       call mcscf_hessv_2e_mo_iaav( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xua), g_ax)
       if (.not.ma_pop_stack(l_xua))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0)
c
c
c
       if (.not.ma_push_get(MT_DBL,(nact*nclosed),'xju',l_xju,k_xju))
     $      call errquit('mcscf_hessv_2e_mo: cannot allocate',0)
       call ga_get(g_x,aoff,aend,1,nclosed,dbl_mb(k_xju),nact)
       call mcscf_hessv_2e_mo_iaia( nbf, nclosed, nact, dm2, g_coul,
     $                              g_exch, dbl_mb(k_xju), g_ax )
       if (.not.ma_pop_stack(l_xju))
     $      call errquit('mcscf_hessv_2e_mo: failed pop stack',0)
c
c
c
       return
       end




c
c
c
c
       subroutine mcscf_hessv_2e_ao00( geom, basis, nbf, nclosed, nact,
     $                               oskel, dm1, g_movecs, g_x, g_tmp,
     $                               g_ax )
#include "global.fh"
#include "mafdecls.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       logical oskel                                ! [input] Symmetry selection
       double precision dm1(nact,nact)              ! [input] 1PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_tmp                                ! [input] Temporary  (nbf * nbf)
       integer g_x                                  ! [input] Argument parameter matrix 
       integer g_ax                                 ! [output] Hessian product (in matrix rep)
c
c
c
       integer g_f1, g_f2, g_f3, g_f4, g_f5, g_f6   ! Atom-blocked Fock matrices
       integer g_d1, g_d2, g_d3, g_d4, g_d5, g_d6   ! Atom-blocked densities
       integer nvir, vlen, voff, aoff, aend, nn
       integer g_tmp1, g_tmp2
c
c
       integer nsets
       parameter(nsets=6)
       integer iv_dens(nsets), iv_fock(nsets)
       double precision jfac(nsets),kfac(nsets)
       data jfac/6*1.d0/
       data kfac/6*-0.5d0/
       double precision tol2e
       data tol2e/1.d-12/
c
c
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
c
c
c
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
       if (.not.ga_duplicate(g_tmp,g_tmp1,'temp1'))
     $      call errquit('mcscf_hessv_2e_ao: cannot duplicate',0)
       if (.not.ga_duplicate(g_tmp,g_tmp2,'temp2'))
     $      call errquit('mcscf_hessv_2e_ao: cannot duplicate',0)
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_zero(g_tmp)
       if (ga_nodeid().eq.0)
     $      call ga_put(g_tmp,aoff,aend,aoff,aend,dm1,nact)
       call ga_sync()
       g_d1 = ga_create_atom_blocked( geom, basis, 'Density 1')
       g_d2 = ga_create_atom_blocked( geom, basis, 'Density 2')
       g_d3 = ga_create_atom_blocked( geom, basis, 'Density 3')
       g_d4 = ga_create_atom_blocked( geom, basis, 'Density 4')
       g_d5 = ga_create_atom_blocked( geom, basis, 'Density 5')
       g_d6 = ga_create_atom_blocked( geom, basis, 'Density 6')
       g_f1 = ga_create_atom_blocked( geom, basis, 'Fock 1')
       g_f2 = ga_create_atom_blocked( geom, basis, 'Fock 2')
       g_f3 = ga_create_atom_blocked( geom, basis, 'Fock 3')
       g_f4 = ga_create_atom_blocked( geom, basis, 'Fock 4')
       g_f5 = ga_create_atom_blocked( geom, basis, 'Fock 5')
       g_f6 = ga_create_atom_blocked( geom, basis, 'Fock 6')
c
c Inactive-Virtual "density"
c                                               iv   t iv
c                                              D  = C.X .C
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_x, voff, nbf, 1, nclosed,
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy( g_tmp2, g_d1 )
       call ga_copy( g_d1, g_d2 )
       call ga_copy( g_d1, g_d3 )
c
c Inactive-Virtual, Active-Virtual "density"    
c                                               ,av   t av
c                                              D   = C.X .d.C
c
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_x, voff, nbf, aoff, aend,
     $                       g_tmp, aoff, aend, aoff, aend,
     $                       g_tmp1, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_tmp2, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, voff, nbf,
     $                       g_tmp2, voff, nbf, 1, nbf,
     $                       g_tmp1, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp1)
       call ga_copy( g_tmp1, g_d4 )
       call ga_dadd( 1.d0, g_d4, 2.d0, g_d1, g_d1 )                  
c
c Inactive-Virtual, Inactive-Active "density"
c                                               ,ia   t      ia
c                                              D   = C.(2-d)X .C
c
       call ga_zero(g_tmp2)
       call diagfill_patch(g_tmp2, 2.d0, aoff, aend )
       call ga_dadd(1.d0, g_tmp2, -1.d0, g_tmp, g_tmp2 )
       call ga_zero(g_tmp1)
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                       g_tmp2, aoff, aend, aoff, aend,
     $                       g_x, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_matmul_patch('n', 't', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp2, aoff, aend, 1, nbf )
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nbf,
     $                       g_tmp1, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp1)
       call ga_copy( g_tmp1, g_d5 )
       call ga_dadd( 1.d0, g_d5, 1.d0, g_d1, g_d1 )
c
c Inactive-Active density
c                                               ia    t ia
c                                              D   = C.X .C
c
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_matmul_patch('n', 't', 1.d0, 0.d0,
     $                      g_x, aoff, aend, 1, nclosed,
     $                      g_movecs, 1, nclosed, 1, nbf,
     $                      g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                      g_movecs, 1, nbf, aoff, aend,
     $                      g_tmp1, aoff, aend, 1, nbf,
     $                      g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy( g_tmp2, g_d5 )
       call ga_dadd( 1.d0, g_d5, 1.d0, g_d2, g_d2 )
c
c
c CA, CA density
c                                               ,,ia    t       ia
c                                              D    = -C.(1-d).X .C
c
       call ga_zero(g_tmp1)
       call diagfill_patch(g_tmp1, 1.d0, aoff, aend )                  
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1)               
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_x, aoff, aend, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp2, aoff, aend, 1, nclosed, 
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy(g_tmp2, g_d6)
c 
c
c
c  Density summary
c
c              iv   ,av   ,ia                 ,av
c       d1 : 2D  + D   + D              d4 : D
c
c             iv   ia                         ia
c       d2 : D  + D                     d5 : D
c
c             iv                              ,,ia
c       d3 : D                          d6 : D
c
c
       iv_dens(1) = g_d1
       iv_dens(2) = g_d2
       iv_dens(3) = g_d3
       iv_dens(4) = g_d4
       iv_dens(5) = g_d5
       iv_dens(6) = g_d6
c
c Fock build
c
       nn = 6
       call ga_zero(g_f1)
       call ga_zero(g_f2)
       call ga_zero(g_f3)
       call ga_zero(g_f4)
       call ga_zero(g_f5)
       call ga_zero(g_f6)
       iv_fock(1) = g_f1
       iv_fock(2) = g_f2
       iv_fock(3) = g_f3
       iv_fock(4) = g_f4
       iv_fock(5) = g_f5
       iv_fock(6) = g_f6
       call fock_2e( geom, basis, nn, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock )
c
c  Symmetrize AO Fock matrices
c
       if (oskel) then
         do i=1,nn
           call sym_symmetrize(geom, basis, .false., iv_fock(i))
         enddo
       endif
c
c Inactive-Virtual contribution
c 
c
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, voff, nbf, 1, nbf,
     $                       g_f1, 1, nbf, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,                      ! Where does this factor come from?
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, voff, nbf, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp2, voff, nbf, 1, nclosed,
     $                     1.d0, g_ax, voff, nbf, 1, nclosed,
     $                           g_ax, voff, nbf, 1, nclosed )
c
c Active-Virtual contribution
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, voff, nbf, 1, nbf,
     $                       g_f2, 1, nbf, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp2, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp2, voff, nbf, aoff, aend,
     $                       g_tmp, aoff, aend, aoff, aend,
     $                       g_tmp1, voff, nbf, aoff, aend )
       call ga_dadd_patch( 1.d0, g_tmp1, voff, nbf, aoff, aend,
     $                     1.d0, g_ax, voff, nbf, aoff, aend,
     $                           g_ax, voff, nbf, aoff, aend )
c
c Inactive-Active contributions
c                                          t   iv
c                                 (2 - d) C.F[D  ].C
c
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f3, 1, nbf, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_zero(g_tmp1)
       call diagfill_patch( g_tmp1, 2.d0, aoff, aend )
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1 )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
c
c                                     t   ,av  
c                                 -d C.F[D   ].C
c
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f4, 1, nbf, 1, nbf,
     $                       g_tmp2, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 1.d0,
     $                       g_tmp2, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp1, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c Inactive-Active, Inactive-Active contribution
c    (note this last section has zero contribution
c    in ROHF theory and has not been debugged against the ROHF
c    Hessian product)
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,                         
     $                       g_movecs, aoff, aend, 1, nbf,                
     $                       g_f6, 1, nbf, 1, nbf,                        
     $                       g_tmp1, aoff, aend, 1, nbf )                 
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,                        
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp2, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f5, 1, nbf, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call diagfill_patch( g_tmp1, 1.d0, aoff, aend )
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1 )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp1, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c
       if (.not.ga_destroy(g_tmp1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_tmp2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d3))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d4))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d5))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_d6))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f3))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f4))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f5))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
       if (.not.ga_destroy(g_f6))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0)
c
c
c
       return
       end




















       subroutine getaoints( nbf, nclosed, nact, nnocc, hlp )
       implicit none
       integer nbf, nclosed, nact, nnocc
       double precision hlp(nbf*nbf,nnocc*2)
       integer i

       open(unit=88,file='JKints',status='old',form='unformatted')
       do i=1,2*nnocc
         call hv_readv(88,(nbf*nbf),hlp(1,i))
C         CALL MOINTS_MATPRINT(NBF,NBF,HLP(1,I))
       enddo
       close(88)
       return
       end

       





c
c  Contribution from eqn A1 (Siegbahn, Almlof, Heiberg and Roos)
c
c
       subroutine hmat_ao_iaia( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)

       integer nvir
       integer i,j,t,u,v,it,ju,tt,uu,vv
       integer iaoff,nn,ij
       double precision dd, xx, ee, yy

       nvir = nbf - nclosed - nact
       iaoff = (nclosed+nact)*nvir
       nn = nbf*nbf
       do i=1,nclosed
         do j=1,nclosed
           ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
           call dcopy(nn,xeri(1,ij),1,xj,1)
           call dcopy(nn,xeri(1,ij+nnocc),1,xk,1)
           if (j.gt.i) call mattrp(nbf,xk)
           do t=1,nact
             it = iaoff + (i-1)*nact + t
             tt = t + nclosed
             do u=1,nact
               ju = iaoff + (j-1)*nact + u
               uu = u + nclosed
               do v=1,nact
                 vv = v + nclosed
                 dd = -dm1(t,v)
                 if (v.eq.t) dd = dd + 1.d0
                 xx = 4.d0*xk(uu,vv) - xk(vv,uu) - xj(vv,uu)
                 
                 ee = -dm1(u,v)
                 if (v.eq.u) ee = ee + 1.d0
                 yy = 4.d0*xk(vv,tt) - xk(tt,vv) - xj(vv,tt)
                 
                 h(it,ju) = h(it,ju) + xx*dd + yy*ee
               enddo
             enddo
           enddo
         enddo
       enddo
                 
       return
       end



c
c  Contribution from eqn A2 (Siegbahn, Almlof, Heiberg and Roos)
c
c
       subroutine hmat_ao_iaiv( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)
c
       integer nvir,i,j,t,v,a,it,ja,vv,aa,vj
       integer iaoff,avoff,nn
       double precision dd, xx

       nvir = nbf - nclosed - nact
       iaoff = (nclosed+nact)*nvir
       avoff = nclosed*nvir
       nn = nbf*nbf
       do j=1,nclosed
         do v=1,nact
           vv = v + nclosed
           vj = (vv*(vv-1))/2 + j
           call dcopy(nn,xeri(1,vj),1,xj,1)
           call dcopy(nn,xeri(1,vj+nnocc),1,xk,1)
           do i=1,nclosed
             do t=1,nact
               it = iaoff + (i-1)*nact + t
               do a=1,nvir
                 aa = a + nclosed + nact
                 ja = (j-1)*nvir + a
                 
                 dd = -dm1(t,v)
                 if (v.eq.t) dd = 2.d0 + dd
                 xx = 4.d0*xk(aa,i) - xk(i,aa) - xj(i,aa)

                 h(it,ja) = h(it,ja) + xx*dd
                 h(ja,it) = h(ja,it) + xx*dd
               enddo
             enddo
           enddo
         enddo
       enddo
                 
       return
       end





c
c  Contribution from eqn A3 (Siegbahn, Almlof, Heiberg and Roos)
c
c
       subroutine hmat_ao_iaav( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)

       integer nvir
       integer i,t,u,v,a,it,ua,tt,uu,vv,aa,vi
       integer iaoff,avoff,nn
       double precision dd, xx

       nvir = nbf - nclosed - nact
       iaoff = (nclosed+nact)*nvir
       avoff = nclosed*nvir
       nn = nbf*nbf
       do i=1,nclosed
         do v=1,nact
           vv = v + nclosed
           vi = (vv*(vv-1))/2 + i
           call dcopy(nn,xeri(1,vi),1,xj,1)
           call dcopy(nn,xeri(1,vi+nnocc),1,xk,1)
           do t=1,nact
             it = iaoff + (i-1)*nact + t
             tt = t + nclosed
             do u=1,nact
               uu = u + nclosed
               do a=1,nvir
                 aa = a + nclosed + nact
                 ua = avoff + (u-1)*nvir + a
                 
                 dd = dm1(u,v)
                 xx = 4.d0*xk(tt,aa) - xk(aa,tt) - xj(tt,aa)

                 h(it,ua) = h(it,ua) + xx*dd
                 h(ua,it) = h(ua,it) + xx*dd
               enddo
             enddo
           enddo
         enddo
       enddo
                 
       return
       end







c
c  Contribution from eqn A4 (Siegbahn, Almlof, Heiberg and Roos)
c
c
       subroutine hmat_ao_iviv( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)
c
       integer nvir,i,j,a,b,ia,jb,bb,aa,nn,ij
       double precision xx

       nvir = nbf - nclosed - nact
       nn = nbf*nbf
       do i=1,nclosed
         do j=1,nclosed
           ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
           call dcopy(nn,xeri(1,ij),1,xj,1)
           call dcopy(nn,xeri(1,ij+nnocc),1,xk,1)
           if (j.gt.i) call mattrp(nbf,xk)
           do a=1,nvir
             ia = (i-1)*nvir + a
             aa = (nclosed+nact) + a
             do b=1,nvir
               jb = (j-1)*nvir + b
               bb = (nclosed+nact) + b
               xx = 4.d0*xk(bb,aa) - xj(bb,aa) - xk(aa,bb)
               h(ia,jb) = h(ia,jb) + 2.d0*xx
             enddo
           enddo
         enddo
       enddo
                 
       return
       end



c
c  Contribution from eqn A5 (Siegbahn, Almlof, Heiberg and Roos)
c
c
       subroutine hmat_ao_ivav( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)
c
       integer nvir,i,t,v,a,b,ia,bb,aa,nn,vi,avoff,vv,tb
       double precision xx,dd

       nvir = nbf - nclosed - nact
       avoff = nclosed*nvir
       nn = nbf*nbf
       do i=1,nclosed
         do v=1,nact
           vv = nclosed + v
           vi = (vv*(vv-1))/2 + i
           call dcopy(nn,xeri(1,vi),1,xj,1)
           call dcopy(nn,xeri(1,vi+nnocc),1,xk,1)
           do a=1,nvir
             ia = (i-1)*nvir + a
             aa = (nclosed+nact) + a
             do t=1,nact
               do b=1,nvir
                 tb = avoff + (t-1)*nvir + b
                 bb = (nclosed+nact) + b
                 dd = dm1(t,v)
                 xx = 4.d0*xk(aa,bb) - xk(bb,aa) - xj(bb,aa)
                 h(ia,tb) = h(ia,tb) + dd*xx
                 h(tb,ia) = h(tb,ia) + dd*xx
               enddo
             enddo
           enddo
         enddo
       enddo
                 
       return
       end







       subroutine hmat_fock( nbf, nclosed, nact, orlen, dm1,
     $                       fi, fa, fg, h )
       implicit none
       integer nbf, nclosed, nact, orlen
       double precision fi(nbf,nbf), fa(nbf,nbf), fg(nbf,nbf)
       double precision dm1(nact,nact)
       double precision h(orlen,orlen)
c
       integer i,j,a,ia,aa,u,t,uu,tt,it,ju,iaoff,nvir,ivoff
       double precision xx
c
       nvir = nbf - nclosed - nact
       iaoff = (nclosed+nact)*nvir
c
c   IA,IA part                            (eqn A1)
c       
       do i=1,nclosed
         do t=1,nact
           it = iaoff + (i-1)*nact + t
           tt = nclosed + t
           do j=1,nclosed
             do u=1,nact
               ju = iaoff + (j-1)*nact + u
               uu = nclosed + u
               xx = dm1(t,u)*fi(i,j)
               if (i.eq.j) xx = xx + 2.d0*fi(tt,uu) +
     $                          2.d0*fa(tt,uu) - fg(tt,uu)
               if (t.eq.u) xx = xx - 2.d0*(fi(i,j) + fa(i,j))
               h(it,ju) = h(it,ju) + xx
             enddo
           enddo
         enddo
       enddo
c
c   IA,IV part                            (eqn A2)
c
       do i=1,nclosed
         do t=1,nact
           it = iaoff + (i-1)*nact + t
           tt = nclosed + t
           do a=1,nvir
             ia = ivoff + (i-1)*nvir + a
             aa = (nclosed+nact) + a
             xx = 2*(fi(aa,tt)+fa(aa,tt)) - 0.5d0*fg(tt,aa)
             h(it,ia) = xx
             h(ia,it) = xx
           enddo
         enddo
       enddo
       return
       end

