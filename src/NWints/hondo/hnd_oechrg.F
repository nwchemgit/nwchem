      subroutine hnd_oechrg(ixyz,iprim,icoef,i_prim,i_gen,li,
     &                      jxyz,jprim,jcoef,j_prim,j_gen,lj,
     &                      kxyz,kprim,kcoef,k_prim,k_gen,lk,
     &                      lxyz,lprim,lcoef,l_prim,l_gen,ll)
c
c $Id$
c
      implicit none
c
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
c
c     Input parameters
c
      integer li,i_prim,i_gen
      integer lj,j_prim,j_gen
      integer lk,k_prim,k_gen
      integer ll,l_prim,l_gen
      double precision ixyz(3),iprim(i_prim),icoef(i_prim,i_gen)
      double precision jxyz(3),jprim(j_prim),jcoef(j_prim,j_gen)
      double precision kxyz(3),kprim(k_prim),kcoef(k_prim,k_gen)
      double precision lxyz(3),lprim(l_prim),lcoef(l_prim,l_gen)
c
c     Local variables
c
      integer ia,jb,jbmax
      double precision ai,aj,arri,aa1,dum
      double precision daexpa,cci,axi,ayi,azi,dtol,rtol,rri,aa,rrk
c
      rtol = 2.30258d00*itol
      dtol = 1.0d+01**(-itol)
c
c     ----- -ij- charge distribution -----
c
      lit = li + 1
      xi = ixyz(1)
      yi = ixyz(2)
      zi = ixyz(3)
      ljt = lj + 1
      xj = jxyz(1)
      yj = jxyz(2)
      zj = jxyz(3)
      rri=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
      expndi = lit.ge.ljt
      if (lit.ge.ljt) then
         xc=xi
         yc=yi
         zc=zi
         dxij=xi-xj
         dyij=yi-yj
         dzij=zi-zj
      else
         xc=xj
         yc=yj
         zc=zj
         dxij=xj-xi
         dyij=yj-yi
         dzij=zj-zi
      endif
c
c     ----- -giao- factors -----
c
      tijx = xi-xj
      tijy = yi-yj
      tijz = zi-zj
      qijx = yi*zj-zi*yj
      qijy = zi*xj-xi*zj
      qijz = xi*yj-yi*xj
c
c     ----- - i- primitive           -----
c
      nij=0
      do 40 ia=1,i_prim
         ai=iprim(ia)
         arri=ai*rri
         axi=ai*xi
         ayi=ai*yi
         azi=ai*zi
         cci=icoef(ia,i_gen)
c
c     ----- - j- primitive           -----
c
         jbmax=j_prim
         if(iieqjj) jbmax=ia
         do 30 jb=1,jbmax
            aj=jprim(jb)
            aa=ai+aj
            aa1=1.0d0/aa
            dum=aj*arri*aa1
            if(dum.gt.rtol) go to 30
            daexpa=cci*jcoef(jb,j_gen)* exp(-dum)*aa1
            dum= abs(daexpa)
            if(dum.le.dtol) go to 30
c
            nij=nij+1
            acharg( 1,nij)= daexpa
            if(iieqjj.and.ia.ne.jb) acharg( 1,nij)=2.0d0*daexpa
            acharg( 2,nij)= aa
            acharg( 3,nij)=(axi+aj*xj)*aa1
            acharg( 4,nij)=(ayi+aj*yj)*aa1
            acharg( 5,nij)=(azi+aj*zj)*aa1
c
   30    continue
   40 continue
c
c     ----- -kl- charge distribution -----
c
      lkt = lk + 1
      xk = kxyz(1)
      yk = kxyz(2)
      zk = kxyz(3)
      lmt = ll + 1
      xl = lxyz(1)
      yl = lxyz(2)
      zl = lxyz(3)
      rri=((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
      expndk = lkt.ge.lmt
      if (lkt.ge.lmt) then
         xd=xk
         yd=yk
         zd=zk
         dxkl=xk-xl
         dykl=yk-yl
         dzkl=zk-zl
      else
         xd=xl
         yd=yl
         zd=zl
         dxkl=xl-xk
         dykl=yl-yk
         dzkl=zl-zk
      endif
c
c     ----- -giao- factors -----
c
      tklx = xk-xl
      tkly = yk-yl
      tklz = zk-zl
      qklx = yk*zl-zk*yl
      qkly = zk*xl-xk*zl
      qklz = xk*yl-yk*xl
c
c     ----- - k- primitive           -----
c
      nkl=0
      do 60 ia=1,k_prim
         ai=kprim(ia)
         arri=ai*rri
         axi=ai*xk
         ayi=ai*yk
         azi=ai*zk
         cci=kcoef(ia,k_gen)
c
c     ----- - l- primitive           -----
c
         jbmax=l_prim
         if(kkeqll) jbmax=ia
         do 50 jb=1,jbmax
            aj=lprim(jb)
            aa=ai+aj
            aa1=1.0d0/aa
            dum=aj*arri*aa1
            if(dum.gt.rtol) go to 50
            daexpa=cci*lcoef(jb,l_gen)* exp(-dum)*aa1
            dum= abs(daexpa)
            if(dum.le.dtol) go to 50
c
            nkl=nkl+1
            bcharg( 1,nkl)= daexpa
            if(kkeqll.and.ia.ne.jb) bcharg( 1,nkl)=2.0d0*daexpa
            bcharg( 2,nkl)= aa
            bcharg( 3,nkl)=(axi+aj*xl)*aa1
            bcharg( 4,nkl)=(ayi+aj*yl)*aa1
            bcharg( 5,nkl)=(azi+aj*zl)*aa1
c
   50    continue
   60 continue
      return
      end
