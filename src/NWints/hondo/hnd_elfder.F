c  Acknowledgment of AI-Assisted Development
c
c  Portions of this source code were developed with the assistance of
c  Claude (Anthropic), an AI assistant, specifically Claude Opus 5
c  (accessed via API/Claude.ai, 2026).
      subroutine hnd_elfder(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,nder,nint,elfder,scr,lscr,
     2  xyzpt, npt)
c
c $Id$
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npt)
      dimension scr(lscr) 
      dimension elfder(*)
c
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      idim = max(nder*3,1)
c
      call hnd_elfder1(xyzi,expi,coefi,i_nprim,i_ngen, Li, xyzj,expj,
     1 coefj, j_nprim, j_ngen, Lj, nder, nint, elfder, xyzpt, npt, idim,
     2 scr(elpt(1)),scr(elpt(2)),scr(elpt(3)),scr(elpt(4)),scr(elpt(5)),
     3 scr(elpt(6)),scr(elpt(7)),scr(elpt(8)),scr(elpt(9)))
c
      return
      end
c
c=============================================================================
c
      subroutine hnd_elfder1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1  expj,coefj, j_nprim, j_ngen, Lj, nder, nint, elfder, xyzpt, npt,
     2  idim,xv,yv,zv,dxv,dyv,dzv,ddxv,ddyv,ddzv)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "stdio.fh"
#include "hnd_tol.fh"
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
c
c---------------------------------------------------------------------------
c NOTE ON SCRATCH LAYOUT CHANGE
c
c   xv/yv/zv/dxv/.../ddzv are now dimensioned (nroots, Li+1, Lj+1)
c   instead of (Li+1, Lj+1, nroots).
c
c   Total element count is unchanged, so the elpt() offsets computed in
c   hnd_pointers.fh remain valid.  These arrays are scratch, private to
c   this routine, so no caller is affected.  The purpose is to make the
c   innermost reduction loop (over iroot) stride-1; see notes below.
c---------------------------------------------------------------------------
      dimension xv(*), yv(*), zv(*)
      dimension dxv(*), dyv(*), dzv(*)
      dimension ddxv(*), ddyv(*), ddzv(*)
c
      dimension w2(maxrys),w4(maxrys),w0(maxrys)
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npt),elfder(nint,idim,npt)
c
c     Automatic arrays: sized at entry, so no fixed upper bound to overflow.
      integer ixa(( Li+1)*(Li+2)/2), iya((Li+1)*(Li+2)/2)
      integer iza((Li+1)*(Li+2)/2)
      integer jxa(( Lj+1)*(Lj+2)/2), jya((Lj+1)*(Lj+2)/2)
      integer jza((Lj+1)*(Lj+2)/2)
c
c     Routine calculates the following integrals for a given shell:
c
c     nder =-2 : Electronic wave function
c     nder =-1 : Electronic density
c     nder = 0 : Electrostatic potential
c     nder = 1 : Electric field
c     nder = 2 : Electric field gradient
c
c     Maximum scratch size needs to be Max(nder*3,1))*(Li+1)*(Lj+1)*((Li+Lj)/2+1)
c     Data format elfder will be (nints,ipts,Max(nder*3,1))
c
c
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*max(nder*3,1)*npt,0.0d0,0,elfder,1)
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      lit = Li + 1
      maxi = lit*(lit+1)/2
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt = Lj + 1
      maxj = ljt*(ljt+1)/2
c
      dxij = xi-xj
      dyij = yi-yj
      dzij = zi-zj
      rr   = dxij*dxij + dyij*dyij + dzij*dzij
c
      nroots = (lit+ljt+nder-2)/2 + 1
      if (nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_elfgrd: need higher Rys root',nroots,INT_ERR)
      endif
c
c     Strides for the (nroots, lit, ljt) scratch layout.
      nr   = nroots
      nrli = nr*lit
c
c     ----- CHANGE 1: hoist getNxyz completely out of all loops -----
c     These depend only on Li,Lj and the component index, so they are
c     loop invariant with respect to ig, jg, ipt and iroot.
      do 15 i = 1,maxi
         call getNxyz(Li,i,Nxyz)
         ixa(i) = Nxyz(1) + 1
         iya(i) = Nxyz(2) + 1
         iza(i) = Nxyz(3) + 1
   15 continue
      do 16 j = 1,maxj
         call getNxyz(Lj,j,Nxyz)
         jxa(j) = Nxyz(1) + 1
         jya(j) = Nxyz(2) + 1
         jza(j) = Nxyz(3) + 1
   16 continue
c
c     ----- i primitive -----
      do 7000 ig = 1,i_nprim
      ai   = expi(ig)
      arri = ai*rr
      axi  = ai*xi
      ayi  = ai*yi
      azi  = ai*zi
      csi  = coefi(ig,i_ngen)
c
c     ----- j primitive -----
      do 6000 jg = 1,j_nprim
      aj  = expj(jg)
      aa  = ai + aj
      aa1 = 1.0d0/aa
      dum = aj*arri*aa1
      if (dum.gt.tol) go to 6000
      fac = exp(-dum)
      csj = coefj(jg,j_ngen)
      ax  = (axi + aj*xj)*aa1
      ay  = (ayi + aj*yj)*aa1
      az  = (azi + aj*zj)*aa1
c
c     ----- density factor -----
c
      dij = csi*fac*csj
c
c     Electronic density integrals and electronic wave
c     function are done differently
c     then the electric field, electric field gradient, and
c     electrostatic potential integrals
c     Switch between the two classes here:
c
      if (nder.ge.0) goto 399
c
c     ----- electronic wave function / density integral -----
c
c     For wave function we simply set j to 0.0d0 and reset the
c     density factor appropriately
c
c=========================================================================
c     nder < 0 : electronic density (-1) / wave function (-2)
c=========================================================================
      ljtl  = ljt
      maxjl = maxj
      if (nder.eq.-2) then
         dij   = csi
         ax    = xi
         ay    = yi
         az    = zi
         aa    = ai
         ljtl  = 1
         maxjl = 1
      endif
c
      do 300 ipt = 1,npt
         x0  = xyzpt(1,ipt)
         y0  = xyzpt(2,ipt)
         z0  = xyzpt(3,ipt)
         px  = x0 - ax
         py  = y0 - ay
         pz  = z0 - az
         dum = aa*(px*px + py*py + pz*pz)
         if (dum.gt.tol) go to 300
c        ----- CHANGE 2: fold fac into dij once, not once per (i,j) -----
         fdij = exp(-dum)*dij
c
c        Single root here, so the scratch is indexed (1,i,j).
         do 270 j = 1,ljtl
            nj = j
            do 260 i = 1,lit
               ni = i
               call hnd_denxyz
               idx = 1 + (i-1)*nr + (j-1)*nrli
               xv(idx) = xint
               yv(idx) = yint
               zv(idx) = zint
  260       continue
  270    continue
c
         ij = 0
         do 290 j = 1,maxjl
            jx = jxa(j)
            jy = jya(j)
            jz = jza(j)
            joffx = (jx-1)*nrli
            joffy = (jy-1)*nrli
            joffz = (jz-1)*nrli
            do 280 i = 1,maxi
               ij = ij + 1
               elfder(ij,1,ipt) = elfder(ij,1,ipt) + fdij
     1            * xv(1 + (ixa(i)-1)*nr + joffx)
     2            * yv(1 + (iya(i)-1)*nr + joffy)
     3            * zv(1 + (iza(i)-1)*nr + joffz)
  280       continue
  290    continue
  300 continue
      goto 6000
c
c=========================================================================
c     nder >= 0 : potential (0) / field (1) / field gradient (2)
c=========================================================================
  399 dijs = dij*pi212*aa1
      aax  = aa*ax
      aay  = aa*ay
      aaz  = aa*az
c
      do 500 ipt = 1,npt
         znuc = 1.0d0
         cx = xyzpt(1,ipt)
         cy = xyzpt(2,ipt)
         cz = xyzpt(3,ipt)
         pcx = ax - cx
         pcy = ay - cy
         pcz = az - cz
         yy  = aa*(pcx*pcx + pcy*pcy + pcz*pcz)
         call hnd_droot
c
c        ----- CHANGE 3: nder branch hoisted out of the i/j fill loops -----
c        The original used  goto (402,401) nder+1  inside loop 410, which
c        relies on the F77 out-of-range computed-GOTO fall-through for
c        nder=2.  Three explicit nests reproduce exactly the same call
c        sequence and store order, with no branch in the loop body:
c           nder=0 : hnd_sxyz
c           nder=1 : hnd_dervxyz(1), hnd_sxyz
c           nder=2 : hnd_dervxyz(2), hnd_dervxyz(1), hnd_sxyz
c        This loop stays SERIAL: hnd_sxyz/hnd_dervxyz pass arguments
c        through common/hnd_xyzder/, so no iteration of it is independent.
c
         if (nder.eq.2) then
c
            do 420 iroot = 1,nroots
               uu = u9(iroot)*aa
               ww = w9(iroot)*znuc
               w0(iroot) = ww
               w2(iroot) = ww*uu*2.0d0
               w4(iroot) = ww*uu*uu*4.0d0
               tt = 1.0d0/(aa + uu)
               tx = sqrt(tt)
               x0 = (aax + uu*cx)*tt
               y0 = (aay + uu*cy)*tt
               z0 = (aaz + uu*cz)*tt
               do 418 i = 1,lit
                  ni = i
                  do 416 j = 1,ljt
                     nj  = j
                     idx = iroot + (i-1)*nr + (j-1)*nrli
                     call hnd_dervxyz(2)
                     ddxv(idx) = xint
                     ddyv(idx) = yint
                     ddzv(idx) = zint
                     call hnd_dervxyz(1)
                     dxv(idx) = xint
                     dyv(idx) = yint
                     dzv(idx) = zint
                     call hnd_sxyz
                     xv(idx) = xint
                     yv(idx) = yint
                     zv(idx) = zint
  416             continue
  418          continue
  420       continue
c
            ij = 0
            do 440 j = 1,maxj
               jox = (jxa(j)-1)*nrli
               joy = (jya(j)-1)*nrli
               joz = (jza(j)-1)*nrli
               do 435 i = 1,maxi
c                 ----- CHANGE 4: base offsets hoisted; iroot is now the
c                 fastest-varying index, so all nine loads below are
c                 stride-1 and the compiler can vectorise unaided.
                  mx = (ixa(i)-1)*nr + jox
                  my = (iya(i)-1)*nr + joy
                  mz = (iza(i)-1)*nr + joz
                  ij = ij + 1
                  dumxx = 0.0d0
                  dumyy = 0.0d0
                  dumzz = 0.0d0
                  dumxy = 0.0d0
                  dumxz = 0.0d0
                  dumyz = 0.0d0
                  do 430 k = 1,nroots
                     xvv = xv(k+mx)
                     yvv = yv(k+my)
                     zvv = zv(k+mz)
                     w4r = w4(k)
c                    ----- CHANGE 5: xvv*yvv*zvv*w2 is one shared term -----
                     t   = xvv*yvv*zvv*w2(k)
                     dumxx = dumxx - t + ddxv(k+mx)*yvv*zvv*w4r
                     dumyy = dumyy - t + ddyv(k+my)*xvv*zvv*w4r
                     dumzz = dumzz - t + ddzv(k+mz)*xvv*yvv*w4r
                     dumxy = dumxy + dxv(k+mx)*dyv(k+my)*zvv*w4r
                     dumxz = dumxz + dxv(k+mx)*dzv(k+mz)*yvv*w4r
                     dumyz = dumyz + dyv(k+my)*dzv(k+mz)*xvv*w4r
  430             continue
                  elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumxx*dijs
                  elfder(ij,2,ipt) = elfder(ij,2,ipt) + dumyy*dijs
                  elfder(ij,3,ipt) = elfder(ij,3,ipt) + dumzz*dijs
                  elfder(ij,4,ipt) = elfder(ij,4,ipt) + dumxy*dijs
                  elfder(ij,5,ipt) = elfder(ij,5,ipt) + dumxz*dijs
                  elfder(ij,6,ipt) = elfder(ij,6,ipt) + dumyz*dijs
  435          continue
  440       continue
c
         elseif (nder.eq.1) then
c
            do 460 iroot = 1,nroots
               uu = u9(iroot)*aa
               ww = w9(iroot)*znuc
               w0(iroot) = ww
               w2(iroot) = ww*uu*2.0d0
               tt = 1.0d0/(aa + uu)
               tx = sqrt(tt)
               x0 = (aax + uu*cx)*tt
               y0 = (aay + uu*cy)*tt
               z0 = (aaz + uu*cz)*tt
               do 458 i = 1,lit
                  ni = i
                  do 456 j = 1,ljt
                     nj  = j
                     idx = iroot + (i-1)*nr + (j-1)*nrli
                     call hnd_dervxyz(1)
                     dxv(idx) = xint
                     dyv(idx) = yint
                     dzv(idx) = zint
                     call hnd_sxyz
                     xv(idx) = xint
                     yv(idx) = yint
                     zv(idx) = zint
  456             continue
  458          continue
  460       continue
c
            ij = 0
            do 480 j = 1,maxj
               jox = (jxa(j)-1)*nrli
               joy = (jya(j)-1)*nrli
               joz = (jza(j)-1)*nrli
               do 475 i = 1,maxi
                  mx = (ixa(i)-1)*nr + jox
                  my = (iya(i)-1)*nr + joy
                  mz = (iza(i)-1)*nr + joz
                  ij = ij + 1
                  dumx = 0.0d0
                  dumy = 0.0d0
                  dumz = 0.0d0
                  do 470 k = 1,nroots
                     xvv = xv(k+mx)
                     yvv = yv(k+my)
                     zvv = zv(k+mz)
                     w2r = w2(k)
                     dumx = dumx + dxv(k+mx)*yvv*zvv*w2r
                     dumy = dumy + xvv*dyv(k+my)*zvv*w2r
                     dumz = dumz + xvv*yvv*dzv(k+mz)*w2r
  470             continue
                  elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumx*dijs
                  elfder(ij,2,ipt) = elfder(ij,2,ipt) + dumy*dijs
                  elfder(ij,3,ipt) = elfder(ij,3,ipt) + dumz*dijs
  475          continue
  480       continue
c
         else
c
            do 490 iroot = 1,nroots
               uu = u9(iroot)*aa
               w0(iroot) = w9(iroot)*znuc
               tt = 1.0d0/(aa + uu)
               tx = sqrt(tt)
               x0 = (aax + uu*cx)*tt
               y0 = (aay + uu*cy)*tt
               z0 = (aaz + uu*cz)*tt
               do 488 i = 1,lit
                  ni = i
                  do 486 j = 1,ljt
                     nj  = j
                     idx = iroot + (i-1)*nr + (j-1)*nrli
                     call hnd_sxyz
                     xv(idx) = xint
                     yv(idx) = yint
                     zv(idx) = zint
  486             continue
  488          continue
  490       continue
c
            ij = 0
            do 496 j = 1,maxj
               jox = (jxa(j)-1)*nrli
               joy = (jya(j)-1)*nrli
               joz = (jza(j)-1)*nrli
               do 495 i = 1,maxi
                  mx = (ixa(i)-1)*nr + jox
                  my = (iya(i)-1)*nr + joy
                  mz = (iza(i)-1)*nr + joz
                  ij = ij + 1
                  dumx = 0.0d0
                  do 494 k = 1,nroots
                     dumx = dumx
     1                    + xv(k+mx)*yv(k+my)*zv(k+mz)*w0(k)
  494             continue
                  elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumx*dijs
  495          continue
  496       continue
c
         endif
c
  500 continue
c
 6000 continue
 7000 continue
c
      return
 9997 format(' in -elfgrd- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots = ',3i3)
      end
