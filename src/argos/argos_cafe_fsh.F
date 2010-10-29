      subroutine argos_cafe_fsh(nangls,indexl,msh,msp,iang,ang,rang,rub,
     + natoms,ndim,igan,isgm,imol,idyn,ichg,chg,xs,fs,ess,lpbc,lpbcs,
     + lupden,lupdti,dera,lseq)
c
c $Id$
c
c     argos_cafe_fsh returns forces and energies for solute angles
c
c     =======================================================
c
c     description of arguments
c     ------------------------
c
c     in: integer nangls     = number of angles to consider
c         integer indexl     = index list
c
c         integer idsh(msh)  = global atom id i
c         integer jdsh(msh)  = global atom id j
c         integer kdsh(msh)  = global atom id k
c         real*8 cdsh(msh,6) = angle force constants
c         real*8 ddsh(msh,6) = angle reference value
c         real*8 rdsh(msh)   = angle value
c
c         integer natoms     = number of atoms in arrays
c         integer ndim       = leading dimension atom arrays
c         integer igan(ndim) = global atom numbers
c         integer imol(ndim) = atom molecule fraction
c         integer idyn(ndim) = atom dynamics type
c         real*8 qs(ndim)    = atomic charges
c         real*8 xs(ndim,3)  = atom coordinates
c
c         logical lupden     = if .true. energies are updated
c
c    out: real*8 fs(ndim,3)  = atom forces (ACCUMULATED)
c         real*8 ush(msh)    = angle energies
c
      implicit none
c
#include "argos_cafe_common.fh"
c
c     declaration of arguments
c     ------------------------
c
      integer msh,msp
      integer iang(msh,4)
      real*8 ang(msh,msp,mset),rang(msh,2),rub(msh,2)
      integer isgm(msa),lseq(mseq)
c
      integer nangls
      integer indexl(nangls)
c
c      integer idsh(msh),jdsh(msh),kdsh(msh)
c      real*8 cdsh(msh,6),ddsh(msh,6),rdsh(msh),ush(msh)
c
      integer natoms,ndim
      integer igan(ndim),imol(ndim),idyn(ndim),ichg(ndim)
      real*8 chg(mqt,mqp,mset)
      real*8 xs(ndim,3),fs(ndim,3),ess(msf,msf,mpe,2)
c
      logical lpbc,lpbcs,lupden,lupdti
      real*8 dera(6,nsatot)
c
c     declaration of local variables
c     ------------------------------
c
      integer i,j,ish,isa,jsa,ksa,isf,jsf,ksf,ibset
      integer ifacu
      real*8 angle,dangle,for,dfor,dfs,phi,cphi,sphi,factu,dercon
      real*8 xsijx,xskjx,xsijy,xskjy,xsijz,xskjz
      real*8 rsij2,rskj2,rsij2i,rskj2i,rsikji,dx(3),eterm
      real*8 qij,xs1,xs2,xs3,rss,rsi,ferfc,fderfc,dfs1,dfs2,dfs3
      real*8 qijp2,qijp3
      real*8 rb,bond,rss2,rssi,rss2i,dbond
c
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
c     solute angles
c     =============
c
      do 1 i=1,nangls
c
c     find index into list of angles
c     ------------------------------
c
      ish=indexl(i)
c
c     find local atom numbers involved in this angle
c     ----------------------------------------------
c
      isa=0
      jsa=0
      ksa=0
c
      do 2 j=1,natoms
      if(iang(ish,1).eq.igan(j)) isa=j
      if(iang(ish,2).eq.igan(j)) jsa=j
      if(iang(ish,3).eq.igan(j)) ksa=j
    2 continue
c
c     get solute molecule numbers involved in this angle
c     --------------------------------------------------
c
      isf=imol(isa)
      jsf=imol(jsa)
      ksf=imol(ksa)
c
c     determine the factor for the energies depending on
c     atoms being dynamic or fixed
c     --------------------------------------------------
c
      ifacu=0
      if(iand(idyn(isa),mdynam).eq.ldynam) ifacu=ifacu+1
      if(iand(idyn(jsa),mdynam).eq.ldynam) ifacu=ifacu+1
      if(iand(idyn(ksa),mdynam).eq.ldynam) ifacu=ifacu+1
c      factu=dble(ifacu)/three
      factu=one
      if(ifacu.eq.0) factu=zero
      if(includ.eq.1) factu=one
c
c     get reference angle and force constant
c     --------------------------------------
c
      if(nfhop.eq.0) then
      angle=ang(ish,1,iset)
      for=ang(ish,2,iset)
      else
      angle=ang(ish,1,lseq(isgm(jsa)))
      for=ang(ish,2,lseq(isgm(jsa)))
      endif
c
c     determine the angle
c     -------------------
c
      xsijx=xs(isa,1)-xs(jsa,1)
      xskjx=xs(ksa,1)-xs(jsa,1)
      xsijy=xs(isa,2)-xs(jsa,2)
      xskjy=xs(ksa,2)-xs(jsa,2)
      xsijz=xs(isa,3)-xs(jsa,3)
      xskjz=xs(ksa,3)-xs(jsa,3)
c
c     periodic boundary conditions
c
      if(lpbc.or.lpbcs) then
      dx(1)=xsijx
      dx(2)=xsijy
      dx(3)=xsijz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xsijx=dx(1)
      xsijy=dx(2)
      xsijz=dx(3)
      dx(1)=xskjx
      dx(2)=xskjy
      dx(3)=xskjz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xskjx=dx(1)
      xskjy=dx(2)
      xskjz=dx(3)
      endif
c
      rsij2=xsijx*xsijx+xsijy*xsijy+xsijz*xsijz
      rskj2=xskjx*xskjx+xskjy*xskjy+xskjz*xskjz
      cphi=xsijx*xskjx+xsijy*xskjy+xsijz*xskjz
      rsij2i=one/rsij2
      rskj2i=one/rskj2
      rsikji=one/sqrt(rsij2*rskj2)
      cphi=cphi*rsikji
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      phi=acos(cphi)
      rang(ish,1)=phi
      dangle=phi-angle
c
c     evaluate energies and forces
c     ----------------------------
c
      rang(ish,2)=half*for*dangle*dangle
      eterm=zero
      if(lupden) then
      ess(isf,isf,2,1)=ess(isf,isf,2,1)+third*factu*rang(ish,2)
      ess(jsf,jsf,2,1)=ess(jsf,jsf,2,1)+third*factu*rang(ish,2)
      ess(ksf,ksf,2,1)=ess(ksf,ksf,2,1)+third*factu*rang(ish,2)
      endif
      eterm=factu*rang(ish,2)
      sphi=sin(phi)
      if(abs(sphi).lt.small) sphi=small
      dfor=for*dangle/sphi
      dfs=dfor*(xskjx*rsikji-xsijx*rsij2i*cphi)
      fs(isa,1)=fs(isa,1)+dfs
      fs(jsa,1)=fs(jsa,1)-dfs
      dfs=dfor*(xsijx*rsikji-xskjx*rskj2i*cphi)
      fs(ksa,1)=fs(ksa,1)+dfs
      fs(jsa,1)=fs(jsa,1)-dfs
      dfs=dfor*(xskjy*rsikji-xsijy*rsij2i*cphi)
      fs(isa,2)=fs(isa,2)+dfs
      fs(jsa,2)=fs(jsa,2)-dfs
      dfs=dfor*(xsijy*rsikji-xskjy*rskj2i*cphi)
      fs(ksa,2)=fs(ksa,2)+dfs
      fs(jsa,2)=fs(jsa,2)-dfs
      dfs=dfor*(xskjz*rsikji-xsijz*rsij2i*cphi)
      fs(isa,3)=fs(isa,3)+dfs
      fs(jsa,3)=fs(jsa,3)-dfs
      dfs=dfor*(xsijz*rsikji-xskjz*rskj2i*cphi)
      fs(ksa,3)=fs(ksa,3)+dfs
      fs(jsa,3)=fs(jsa,3)-dfs
c
c     evaluate the hessian
c     -------------------
c
cx      if(ihess.gt.0) then
cx      endif
c
      if(lupdti) then
c
c     for thermodynamic perturbations evaluate the energies
c     using 'perturbed' parameters of set 2 and/or set 3
c     -----------------------------------------------------
c
      if(ip2(20))
     + ep2(1)=ep2(1)-eterm+factu*half*ang(ish,2,2)*(phi-ang(ish,1,2))**2
      if(ip3(20))
     + ep3(1)=ep3(1)-eterm+factu*half*ang(ish,2,3)*(phi-ang(ish,1,3))**2
c
c     for thermodynamic integrations evaluate the derivative
c     ------------------------------------------------------
c
      if(ith(20)) then
      dercon=dangle*(half*dangle*ang(ish,2,4)-for*ang(ish,1,4))
      deriv(20,1)=deriv(20,1)+dercon
      if(npgdec.gt.1) then
      dera(5,iang(ish,1))=dera(5,iang(ish,1))+third*dercon
      dera(5,iang(ish,2))=dera(5,iang(ish,2))+third*dercon
      dera(5,iang(ish,3))=dera(5,iang(ish,3))+third*dercon
      endif
      endif
      endif
c
c     particle-mesh Ewald correction energy and forces
c     ------------------------------------------------
c
      if(ipme.ne.0) then
      isf=imol(isa)
      ksf=imol(ksa)
      qij=chg(ichg(isa),1,1)*chg(ichg(ksa),1,1)
      xs1=xs(isa,1)-xs(ksa,1)
      xs2=xs(isa,2)-xs(ksa,2)
      xs3=xs(isa,3)-xs(ksa,3)
      rss=sqrt(xs1**2+xs2**2+xs3**2)
      rsi=one/rss
      ferfc=one-erfc(ealpha*rss)
      fderfc=-(ealpha*derfc(ealpha*rss))
      epmecs=epmecs-ferfc*qij*rsi
      if(lupden) then
      ess(isf,isf,9,1)=ess(isf,isf,9,1)-half*ferfc*qij*rsi
      ess(ksf,ksf,9,1)=ess(ksf,ksf,9,1)-half*ferfc*qij*rsi
      if(ipert2) then
      qijp2=chg(ichg(isa),1,2)*chg(ichg(ksa),1,2)
      ess(isf,isf,10,1)=ess(isf,isf,10,1)-half*ferfc*qijp2*rsi
      ess(ksf,ksf,10,1)=ess(ksf,ksf,10,1)-half*ferfc*qijp2*rsi
      endif
      if(ipert2) then
      qijp3=chg(ichg(isa),1,3)*chg(ichg(ksa),1,3)
      ess(isf,isf,11,1)=ess(isf,isf,11,1)-half*ferfc*qijp3*rsi
      ess(ksf,ksf,11,1)=ess(ksf,ksf,11,1)-half*ferfc*qijp3*rsi
      endif
      endif
      dfor=-(qij*rsi*rsi*(ferfc*rsi-fderfc))
      dfs1=dfor*xs1
      dfs2=dfor*xs2
      dfs3=dfor*xs3
      fs(isa,1)=fs(isa,1)-dfs1
      fs(ksa,1)=fs(ksa,1)+dfs1
      fs(isa,2)=fs(isa,2)-dfs2
      fs(ksa,2)=fs(ksa,2)+dfs2
      fs(isa,3)=fs(isa,3)-dfs3
      fs(ksa,3)=fs(ksa,3)+dfs3
      vpmeb(1)=vpmeb(1)+dfs1*xs1
      vpmeb(2)=vpmeb(2)+dfs2*xs1
      vpmeb(3)=vpmeb(3)+dfs3*xs1
      vpmeb(4)=vpmeb(4)+dfs2*xs2
      vpmeb(5)=vpmeb(5)+dfs3*xs2
      vpmeb(6)=vpmeb(6)+dfs3*xs3
      endif
c
    1 continue
c
      if(iffld.eq.2) then
c
c     Urey-Bradley solute angles
c     ==========================
c
      do 3 i=1,nangls
c
c     find index into list of angles
c     ------------------------------
c
      ish=indexl(i)
c
c     find local atom numbers involved in this angle
c     ----------------------------------------------
c
      isa=0
      ksa=0
c
      do 4 j=1,natoms
      if(iang(ish,1).eq.igan(j)) isa=j
      if(iang(ish,3).eq.igan(j)) ksa=j
    4 continue
      isf=imol(isa)
      ksf=imol(ksa)
c
      if(nfhop.eq.0) then
      ibset=iset
      else
      ibset=lseq(isgm(isa))
      endif
c
c     find reference bond length and force constant
c
      bond=ang(ish,3,ibset)
      for=ang(ish,4,ibset)
c
      if(for.gt.0.0d0) then
c
c     determine actual distance between the atoms
c
      xs1=xs(isa,1)-xs(ksa,1)
      xs2=xs(isa,2)-xs(ksa,2)
      xs3=xs(isa,3)-xs(ksa,3)
c
c     periodic boundary conditions
c
      if(lpbc.or.lpbcs) then
      dx(1)=xs1
      dx(2)=xs2
      dx(3)=xs3
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xs1=dx(1)
      xs2=dx(2)
      xs3=dx(3)
      endif
c
      rss2=xs1**2+xs2**2+xs3**2
      if(rss2.gt.tiny) then
      rss=sqrt(rss2)
      rssi=one/rss
      rss2i=rssi*rssi
      else
      rss=zero
      rssi=one
      rss2i=one
      endif
c
c     determine fraction of energy to be counted
c     this depends on the atoms being dynamic or fixed
c
      factu=zero
      if(iand(idyn(isa),mdynam).eq.ldynam.or.
     + iand(idyn(ksa),mdynam).eq.ldynam) factu=one
c      if((iand(idyn(isa),mdynam).eq.ldynam.and.
c     + iand(idyn(ksa),mdynam).ne.ldynam) .or.
c     + (iand(idyn(isa),mdynam).ne.ldynam.and.
c     + iand(idyn(ksa),mdynam).eq.ldynam)) factu=half
      if(includ.eq.1) factu=one
c
      dbond=rss-bond
c
c     evaluate energies and forces
c
      rb=half*for*dbond*dbond
      eterm=zero
      if(lupden) then
      ess(isf,isf,13,1)=ess(isf,isf,13,1)+half*factu*rb
      ess(ksf,ksf,13,1)=ess(ksf,ksf,13,1)+half*factu*rb
      endif
      eterm=factu*rb
      dfor=for*dbond*rssi
      dfs1=dfor*xs1
      dfs2=dfor*xs2
      dfs3=dfor*xs3
      fs(isa,1)=fs(isa,1)-dfs1
      fs(ksa,1)=fs(ksa,1)+dfs1
      fs(isa,2)=fs(isa,2)-dfs2
      fs(ksa,2)=fs(ksa,2)+dfs2
      fs(isa,3)=fs(isa,3)-dfs3
      fs(ksa,3)=fs(ksa,3)+dfs3
c
      rub(ish,1)=rss
      rub(ish,2)=eterm
c
      if(lupdti) then
c
c     for thermodynamic perturbations evaluate the energies using
c     the 'perturbed' parameters in set 2 and/or 3
c
      if(ip2(18))
     + ep2(1)=ep2(1)-eterm+factu*half*ang(ish,4,2)*(rss-ang(ish,3,2))**2
      if(ip3(18))
     + ep3(1)=ep3(1)-eterm+factu*half*ang(ish,4,3)*(rss-ang(ish,3,3))**2
c
c     for thermodynamic integrations evaluate the derivative
c
      if(ith(18)) then
      dercon=dbond*(half*dbond*ang(ish,4,4)-for*ang(ish,3,4))
      deriv(20,1)=deriv(20,1)+dercon
      if(npgdec.gt.1) then
      dera(5,iang(ish,1))=dera(5,iang(ish,1))+half*dercon
      dera(5,iang(ish,3))=dera(5,iang(ish,3))+half*dercon
      endif
      endif
c
      endif
c
      endif
c
    3 continue
c
      endif
c
      return
      end
