c     
c     =========================================================
c     Routines to generate CASSCF orbital Hessian
c     explicitly in-core. For debugging only.
c     =========================================================
c     
c     
      subroutine hmat( nbf, nclosed, nact, orlen, dm1, dm2,
     $                 g_ifock, g_afock, g_gfock,
     $                 g_coul, g_exch, h )
*
* $Id$
*
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
      data odofock, odoao, odomo/.true.,.true.,.true./
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





c
c  One-electron (Fock) contributions 
c  to the Hessian
c
       subroutine hmat_fock( nbf, nclosed, nact, orlen, dm1,
     $                       fi, fa, fg, h )
       implicit none
       integer nbf, nclosed, nact, orlen
       double precision fi(nbf,nbf), fa(nbf,nbf), fg(nbf,nbf)
       double precision dm1(nact,nact)
       double precision h(orlen,orlen)
c
       integer i,j,a,ia,aa,u,t
       integer ta,uu,tt,it,ju,ua
       integer b,bb,ja,ib,ub
       integer nvir,iaoff,avoff,ivoff
       double precision xx
c
       nvir = nbf - nclosed - nact
       ivoff = 0
       avoff = nclosed*nvir
       iaoff = avoff + nact*nvir
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
             h(it,ia) = h(it,ia) + xx
             h(ia,it) = h(ia,it) + xx
           enddo
         enddo
       enddo
c
c   IA,AV part                            (eqn A3)
c
       do i=1,nclosed
         do t=1,nact
           it = iaoff + (i-1)*nact + t
           tt = nclosed + t
           do a=1,nvir
             ia = ivoff + (i-1)*nvir + a
             ta = avoff + (t-1)*nvir + a
             aa = (nclosed+nact) + a
             xx = fi(aa,i) + fa(aa,i)
             h(it,ta) = h(it,ta) + xx
             h(ta,it) = h(ta,it) + xx
             do u=1,nact
               ua = avoff + (u-1)*nvir + a
               uu = nclosed + u
               xx = -dm1(t,u)*fi(i,aa)
               h(it,ua) = h(it,ua) + xx
               h(ua,it) = h(ua,it) + xx
             enddo
           enddo
         enddo
       enddo
c
c   IV,IV part                            (eqn A4)
c
       do i=1,nclosed
         do a=1,nvir
           ia = ivoff + (i-1)*nvir + a
           aa = nclosed + nact + a
           do j=1,nclosed
             ja = ivoff + (j-1)*nvir + a
             xx = -2.d0*(fi(i,j) + fa(i,j))
             h(ia,ja) = h(ia,ja) + xx
           enddo
           do b=1,nvir
             ib = ivoff + (i-1)*nvir + b
             bb = nclosed + nact + b
             xx = 2.d0*(fi(aa,bb) + fa(aa,bb))
             h(ia,ib) = h(ia,ib) + xx
           enddo
         enddo
       enddo  
c
c   IV,AV                                 (eqn A5)
c
       do i=1,nclosed
         do a=1,nvir
           ia = ivoff + (i-1)*nvir + a
           aa = nclosed + nact + a
           do t=1,nact
             ta = avoff + (t-1)*nvir + a
             tt = nclosed + t
             xx = -(fi(tt,i) + fa(tt,i)) - 0.5d0*fg(tt,i)
             h(ia,ta) = h(ia,ta) + xx
             h(ta,ia) = h(ta,ia) + xx
           enddo
         enddo
       enddo
c
c   AV,AV                                 (eqn A6)
c
       do t=1,nact
         do a=1,nvir
           ta = avoff + (t-1)*nvir + a
           tt = nclosed + t
           aa = nclosed + nact + a
           do u=1,nact
             uu = nclosed + u
             ua = avoff + (u-1)*nvir + a
             xx = -fg(tt,uu)
             h(ta,ua) = h(ta,ua) + xx
             do b=1,nvir
               ub = avoff + (u-1)*nvir + b
               bb = nclosed + nact + b
               xx = dm1(t,u)*fi(aa,bb)
               h(ta,ub) = h(ta,ub) + xx
             enddo
           enddo
         enddo
       enddo
         
       return
       end

