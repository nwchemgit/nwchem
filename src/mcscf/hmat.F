c     
c     =========================================================
c     Set of routines to generate CASSCF orbital Hessian
c     explicitly. For debugging.
c     =========================================================
c     
c     
      subroutine hmat( nbf, nclosed, nact, orlen, dm1, dm2,
     $                 g_coul, g_exch, h )
      implicit none
      integer nbf, nclosed, nact, orlen
      double precision dm1(nact,nact)
      double precision dm2(nact,nact,nact,nact)
      integer g_coul, g_exch
      double precision h(orlen,orlen)
      double precision xx(1000),yy(1000)
      double precision hlp(1000000)
      integer i,j,nocc,nnocc
      logical odofock, odoao, odomo
      data odofock, odoao, odomo/.false.,.true.,.false./
c     
      call dfill((orlen*orlen),0.d0,h,1)
      nocc = nclosed + nact
      nnocc = (nocc*(nocc+1))/2
c     
c     AO section: P and K supermatrix contracted against 1-pdm
c     
      if (odoao) then
        call getaoints( nbf, nclosed, nact, nnocc, hlp )
        call hmat_ao_iaia( nbf, nclosed, nact, nnocc, orlen,
     $                     dm1, hlp, xx, yy, h )
      endif
c     
c     MO section: MO integrals contracted against 2-pdm
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
c     Scale
c     
      call dscal((orlen*orlen),2.d0,h,1)
c     
c     Write to file
c     
      open(unit=55,file='hessX',form='unformatted',
     $     status='unknown')
      do i=1,orlen
        write(55) (h(j,i),j=1,orlen)
      enddo
      close(55)
      
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
       call mcscf_hessv_2e00( geom, basis, nbf, nclosed, nact,
     $                        oskel,
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
c$$$       call mcscf_hessv_2e_ao( geom, basis, nbf, nclosed, nact, oskel, 
c$$$     $                         dm1, g_movecs, g_tmp1, g_tmp2, g_tmp3 )
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









       subroutine getaoints( nbf, nclosed, nact, nnocc, hlp )
       implicit none
       integer nbf, nclosed, nact, nnocc
       double precision hlp(nbf*nbf,nnocc*2)
       integer i

       open(unit=88,file='JKints',status='old',form='unformatted')
       do i=1,2*nnocc
         call hv_readv(88,(nbf*nbf),hlp(1,i))
         call moints_matprint(nbf,nbf,hlp(1,i))
       enddo
       close(88)
       return
       end

       






       subroutine hmat_ao_iaia( nbf, nclosed, nact, nnocc, orlen,
     $                          dm1, xeri, xk, xj, h )
       implicit none
       integer nbf, nclosed, nact, nnocc, orlen
       double precision dm1(nact,nact), xeri(nbf*nbf,nnocc*2)
       double precision xj(nbf,nbf), xk(nbf,nbf)
       double precision h(orlen,orlen)

       integer nvir
       integer i,j,t,u,v,it,ju
       integer iaoff,nn,ij
       double precision dd, xx

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
             do u=1,nact
               ju = iaoff + (j-1)*nact + u
               do v=1,nact
                 dd = -dm1(t,v)
                 if (v.eq.t) dd = dd + 1.d0
                 xx = 4.d0*xk(v,u) - xk(u,v) - xj(u,v)
               enddo
             enddo
           enddo
         enddo
       enddo
                 
       return
       end
