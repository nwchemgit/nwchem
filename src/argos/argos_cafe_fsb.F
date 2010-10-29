      subroutine argos_cafe_fsb(nbonds,indexl,msb,msp,ibnd,bnd,rbnd,
     + natoms,ndim,igan,isgm,imol,idyn,ichg,chg,xs,fs,ess,lpbc,lpbcs,
     + lupden,lupdti,dera,lseq)
c
c $Id$
c
c     argos_cafe_fsb returns forces and energies for solute bonds
c
c     =====================================================
c
c     description of arguments
c     ------------------------
c
c     in: integer nbonds     = number of bonds to consider
c         integer indexl     = index list
c
c         integer idsb(msb)  = global atom id i
c         integer jdsb(msb)  = global atom id j
c         integer isbs(msb)  = bond type
c         real*8 cdsb(msb,6) = bond force constants
c         real*8 ddsb(msb,6) = bond reference value
c         real*8 rdsb(msb)   = bond value
c
c         integer natoms     = number of atoms in arrays
c         integer ndim       = leading dimension atom arrays
c         integer igan(ndim) = global atom numbers
c         integer imol(ndim) = atom molecule fraction
c         integer idyn(ndim) = atom dynamics type
c         real*8 qs(ndim)    = atomic charges
c         real*8 xs(ndim,3)  = atom coordinates
c
c    out: real*8 fs(ndim,3)  = atom forces (ACCUMULATED)
c         real*8 usb(msb)    = bond energies
c
      implicit none
c
#include "argos_cafe_common.fh"
c
c     declaration of arguments
c     ------------------------
c
      integer msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,mset),rbnd(msb,2)
      integer isgm(msa),lseq(mseq)
c
      integer nbonds
      integer indexl(nbonds)
c
      logical lpbc,lpbcs,lupden,lupdti
      real*8 dera(6,nsatot)
c
      integer natoms,ndim
      integer igan(ndim),imol(ndim),idyn(ndim),ichg(ndim)
      real*8 chg(mqt,mqp,mset)
      real*8 xs(ndim,3),fs(ndim,3)
      real*8 dx(3)
c
c     declaration of local variables
c     ------------------------------
c
      integer i,j,isb,isa,jsa,isf,jsf,ibset
      real*8 factu,dercon,qij,ferfc,fderfc,qijp2,qijp3
      real*8 bond,dbond,for,dfor,dfs1,dfs2,dfs3,eterm
      real*8 xs1,xs2,xs3,rss,rss2,rssi,rss2i,ess(msf,msf,mpe,2)
c
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
c     solute bonds
c     ============
c
cx      write(*,'(4i7)') (i,(ibnd(i,j),j=1,3),i=1,msb)
c
cx      write(*,'(10i7)') (igan(j),j=1,natoms)
cx      write(*,'(a,i7)') 'bonds  ',nbonds
cx      write(*,'(10i7)') (indexl(j),j=1,nbonds)
c
      do 1 i=1,nbonds
c
c     find index into list of bonds
c
      isb=indexl(i)
c
c     find local atom numbers involved in this bond
c
      isa=0
      jsa=0
      do 2 j=1,natoms
      if(ibnd(isb,1).eq.igan(j)) isa=j
      if(ibnd(isb,2).eq.igan(j)) jsa=j
    2 continue
c
      if(nfhop.eq.0) then
      ibset=iset
      else
      ibset=lseq(isgm(isa))
      endif
c
c      write(*,'(a,5i5)') 'bond ',i,nbonds,isb,isa,jsa
c
c     find solute molecule numbers involved in this constrained
c
      isf=imol(isa)
      jsf=imol(jsa)
c      write(*,'(a,6i5)') 'bond ',i,nbonds,isa,jsa,isf,jsf
c
c     determine actual distance between the atoms
c
      xs1=xs(isa,1)-xs(jsa,1)
      xs2=xs(isa,2)-xs(jsa,2)
      xs3=xs(isa,3)-xs(jsa,3)
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
      rbnd(isb,1)=rss
c
c     if bond not constrained or pme
c
      if(iand(ibnd(isb,3),icnstr).eq.0.or.ipme.ne.0) then
c
c     if bond not constrained
c
      if(iand(ibnd(isb,3),icnstr).eq.0) then
c
c     determine fraction of energy to be counted
c     this depends on the atoms being dynamic or fixed
c
      factu=zero
      if(iand(idyn(isa),mdynam).eq.ldynam.or.
     + iand(idyn(jsa),mdynam).eq.ldynam) factu=one
c      if((iand(idyn(isa),mdynam).eq.ldynam.and.
c     + iand(idyn(jsa),mdynam).ne.ldynam) .or.
c     + (iand(idyn(isa),mdynam).ne.ldynam.and.
c     + iand(idyn(jsa),mdynam).eq.ldynam)) factu=half
      if(includ.eq.1) factu=one
c
c     find reference bond length and force constant
c
      bond=bnd(isb,1,ibset)
      for=bnd(isb,2,ibset)
c
      dbond=rss-bond
c
c     evaluate energies and forces
c
      rbnd(isb,2)=half*for*dbond*dbond
      eterm=zero
      if(lupden) then
      ess(isf,isf,1,1)=ess(isf,isf,1,1)+half*factu*rbnd(isb,2)
      ess(jsf,jsf,1,1)=ess(jsf,jsf,1,1)+half*factu*rbnd(isb,2)
      endif
      eterm=factu*rbnd(isb,2)
      dfor=for*dbond*rssi
      dfs1=dfor*xs1
      dfs2=dfor*xs2
      dfs3=dfor*xs3
      fs(isa,1)=fs(isa,1)-dfs1
      fs(jsa,1)=fs(jsa,1)+dfs1
      fs(isa,2)=fs(isa,2)-dfs2
      fs(jsa,2)=fs(jsa,2)+dfs2
      fs(isa,3)=fs(isa,3)-dfs3
      fs(jsa,3)=fs(jsa,3)+dfs3
c
c     evaluate hessian
 
cx      if(ihess.gt.0) then
cx      isag=igan(isa)
cx      jsag=igan(jsa)
c
cx      hess=for*(one-bond*rssi*(one+xs1*xs1*rss2i))
cx      hs(isa,1,1,isag)=hs(isa,1,1,isag)+hess
cx      hs(isa,1,1,jsag)=hs(isa,1,1,jsag)-hess
cx      hs(jsa,1,1,jsag)=hs(jsa,1,1,jsag)-hess
cx      hs(jsa,1,1,isag)=hs(jsa,1,1,isag)+hess
c
cx      hess=for*(one-bond*rssi*(one+xs2*xs2*rss2i))
cx      hs(isa,2,2,isag)=hs(isa,2,2,isag)+hess
cx      hs(isa,2,2,jsag)=hs(isa,2,2,jsag)-hess
cx      hs(jsa,2,2,jsag)=hs(jsa,2,2,jsag)-hess
cx      hs(jsa,2,2,isag)=hs(jsa,2,2,isag)+hess
c
cx      hess=for*(one-bond*rssi*(one+xs3*xs3*rss2i))
cx      hs(isa,3,3,isag)=hs(isa,3,3,isag)+hess
cx      hs(isa,3,3,jsag)=hs(isa,3,3,jsag)-hess
cx      hs(jsa,3,3,jsag)=hs(jsa,3,3,jsag)-hess
cx      hs(jsa,3,3,isag)=hs(jsa,3,3,isag)+hess
c
cx      hess=for*bond*xs1*xs2*rss2i*rssi
cx      hs(isa,1,2,isag)=hs(isa,1,2,isag)+hess
cx      hs(isa,2,1,isag)=hs(isa,2,1,isag)+hess
cx      hs(isa,1,2,jsag)=hs(isa,1,2,jsag)-hess
cx      hs(isa,2,1,jsag)=hs(isa,2,1,jsag)-hess
cx      hs(jsa,1,2,jsag)=hs(jsa,1,2,jsag)-hess
cx      hs(jsa,2,1,jsag)=hs(jsa,2,1,jsag)-hess
cx      hs(jsa,1,2,isag)=hs(jsa,1,2,isag)+hess
cx      hs(jsa,2,1,isag)=hs(jsa,2,1,isag)+hess
c
cx      hess=for*bond*xs1*xs3*rss2i*rssi
cx      hs(isa,1,3,isag)=hs(isa,1,3,isag)+hess
cx      hs(isa,3,1,isag)=hs(isa,3,1,isag)+hess
cx      hs(isa,1,3,jsag)=hs(isa,1,3,jsag)-hess
cx      hs(isa,3,1,jsag)=hs(isa,3,1,jsag)-hess
cx      hs(jsa,1,3,jsag)=hs(jsa,1,3,jsag)-hess
cx      hs(jsa,3,1,jsag)=hs(jsa,3,1,jsag)-hess
cx      hs(jsa,1,3,isag)=hs(jsa,1,3,isag)+hess
cx      hs(jsa,3,1,isag)=hs(jsa,3,1,isag)+hess
c
cx      hess=for*bond*xs2*xs3*rss2i*rssi
cx      hs(isa,2,3,isag)=hs(isa,2,3,isag)+hess
cx      hs(isa,3,2,isag)=hs(isa,3,2,isag)+hess
cx      hs(isa,2,3,jsag)=hs(isa,2,3,jsag)-hess
cx      hs(isa,3,2,jsag)=hs(isa,3,2,jsag)-hess
cx      hs(jsa,2,3,jsag)=hs(jsa,2,3,jsag)-hess
cx      hs(jsa,3,2,jsag)=hs(jsa,3,2,jsag)-hess
cx      hs(jsa,2,3,isag)=hs(jsa,2,3,isag)+hess
cx      hs(jsa,3,2,isag)=hs(jsa,3,2,isag)+hess
c
cx      endif
c
      if(lupdti) then
c
c     for thermodynamic perturbations evaluate the energies using
c     the 'perturbed' parameters in set 2 and/or 3
c
      if(ip2(18))
     + ep2(1)=ep2(1)-eterm+factu*half*bnd(isb,2,2)*(rss-bnd(isb,1,2))**2
      if(ip3(18))
     + ep3(1)=ep3(1)-eterm+factu*half*bnd(isb,2,3)*(rss-bnd(isb,1,3))**2
c
c     for thermodynamic integrations evaluate the derivative
c
      if(ith(18)) then
      dercon=dbond*(half*dbond*bnd(isb,2,4)-for*bnd(isb,1,4))
      deriv(18,1)=deriv(18,1)+dercon
      if(npgdec.gt.1) then
      dera(5,ibnd(isb,1))=dera(5,ibnd(isb,1))+half*dercon
      dera(5,ibnd(isb,2))=dera(5,ibnd(isb,2))+half*dercon
      endif
      endif
c
      endif
      endif
c
      if(ipme.ne.0) then
      qij=chg(ichg(isa),1,ibset)*chg(ichg(jsa),1,ibset)
      ferfc=one-erfc(ealpha*rss)
      fderfc=-(ealpha*derfc(ealpha*rss))
      epmecs=epmecs-ferfc*qij*rssi
      if(lupden) then
      ess(isf,isf,9,1)=ess(isf,isf,9,1)-half*ferfc*qij*rssi
      ess(jsf,jsf,9,1)=ess(jsf,jsf,9,1)-half*ferfc*qij*rssi
      if(ipert2) then
      qijp2=chg(ichg(isa),1,2)*chg(ichg(jsa),1,2)
      ess(isf,isf,10,1)=ess(isf,isf,10,1)-half*ferfc*qijp2*rssi
      ess(jsf,jsf,10,1)=ess(jsf,jsf,10,1)-half*ferfc*qijp2*rssi
      endif
      if(ipert3) then
      qijp3=chg(ichg(isa),1,3)*chg(ichg(jsa),1,3)
      ess(isf,isf,11,1)=ess(isf,isf,11,1)-half*ferfc*qijp3*rssi
      ess(jsf,jsf,11,1)=ess(jsf,jsf,11,1)-half*ferfc*qijp3*rssi
      endif
      endif
      dfor=-(qij*rssi*rssi*(ferfc*rssi-fderfc))
      dfs1=dfor*xs1
      dfs2=dfor*xs2
      dfs3=dfor*xs3
      fs(isa,1)=fs(isa,1)-dfs1
      fs(jsa,1)=fs(jsa,1)+dfs1
      fs(isa,2)=fs(isa,2)-dfs2
      fs(jsa,2)=fs(jsa,2)+dfs2
      fs(isa,3)=fs(isa,3)-dfs3
      fs(jsa,3)=fs(jsa,3)+dfs3
      vpmeb(1)=vpmeb(1)+dfs1*xs1
      vpmeb(2)=vpmeb(2)+dfs2*xs1
      vpmeb(3)=vpmeb(3)+dfs3*xs1
      vpmeb(4)=vpmeb(4)+dfs2*xs2
      vpmeb(5)=vpmeb(5)+dfs3*xs2
      vpmeb(6)=vpmeb(6)+dfs3*xs3
      endif
c
      endif
c
    1 continue
c
      return
      end
