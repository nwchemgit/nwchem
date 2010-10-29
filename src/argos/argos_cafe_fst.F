      subroutine argos_cafe_fst(nthrds,indexl,mst,idtn,vdw,chg,
     + natoms,ndim,
     + iatt,igan,isgm,imol,idyn,ichg,isss,xs,fs,ess,lpbc,lpbcs,dera,
     + esa,lseq)
c
c $Id$
c
c     argos_cafe_fst returns forces and energies for solute third neighbors
c
c     ================================================================
c
c     description of arguments
c     ------------------------
c
c     in: integer nthrds       = number of third neighbor pairs to consider
c         integer indexl       = index list
c
c         integer idst(mst)    = global atom id i
c         integer jdst(mst)    = global atom id j
c
c         real*8 cb6(matt,matt,3)
c         real*8 cb12(matt,matt,3)
c
c         real*8 cdso(mso,6)   = dihedral angle force constants
c         real*8 ddso(mso,6)   = dihedral angle reference value
c         real*8 rdso(mso)     = dihedral angle value
c
c         integer natoms       = number of atoms in arrays
c         integer ndim         = leading dimension atom arrays
c         integer igan(ndim)   = global atom numbers
c         integer imol(ndim)   = atom molecule fraction
c         integer idyn(ndim)   = atom dynamics type
c         real*8 qs(ndim,4)    = atomic charges
c         real*8 xs(ndim,3)    = atom coordinates
c
c    out: real*8 fs(ndim,3)    = atom forces (ACCUMULATED)
c         real*8 uso(mso)      = dihedral angle energies
c
      implicit none
c
#include "argos_cafe_common.fh"
c
c     declaration of arguments
c     ------------------------
c
      integer nthrds
      integer indexl(nthrds)
      integer isgm(msa),lseq(mseq)
c
      integer mst
      integer idtn(0:mst,2)
c      integer idst(mst),jdst(mst)
c      real*8 cb6(mat,mat,6),cb12(mat,mat,6)
c
      integer natoms,ndim
      real*8 vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer iatt(ndim),igan(ndim),imol(ndim),idyn(ndim),ichg(ndim)
      integer isss(msa)
      real*8 xs(ndim,3),fs(ndim,3)
      real*8 ess(msf,msf,mpe,2),esa(nsa)
c
      real*8 dera(6,nsatot)
c
      logical lpbc,lpbcs
c
c     declaration of local variables
c     ------------------------------
c
      integer i,j,ist,isa,jsa,isf,jsf,isad,istt,jstt,isrx
      real*8 factu,etermq,eterml
      real*8 rxx,rxy,rxz,r2,r2i,r1i,dfs,dfsu,r6i,c6,c12,cf6,cf12,q
      real*8 ferfc,fderfc,qij,dfsp,c64,c124,dercon,qi,qi4,qj,qj4,q14
      real*8 dx(3)
c
c      real*8 ang,dang,dangp,for,dfor,dfs,phi,cphi,cphii,sphi,sphii
c      real*8 rmul,factu,dercon
c      real*8 xsijx,xskjx,xsijy,xskjy,xsijz,xskjz
c      real*8 xsklx,xsjlx,xskly,xsjly,xsklz,xsjlz
c      real*8 xsikx,xsiky,xsikz,xmx,xmy,xmz,xnx,xny,xnz,xdx,xdy,xdz
c      real*8 xex,xey,xez,xox,xoy,xoz,xpx,xpy,xpz
c      real*8 dfsix,dfsiy,dfsiz,dfsjx,dfsjy,dfsjz
c      real*8 dfskx,dfsky,dfskz,dfslx,dfsly,dfslz
c      real*8 rm2i,rn2i,rmni,s,rpa
c      real*8 rsij2,rskj2,rsij2i,rskj2i,rsikji,rkj,rkjo,rkjp
c      real*8 qij,xs1,xs2,xs3,rss,rsi,ferfc,fderfc,dfs1,dfs2,dfs3
c
#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
      isad=0
      dfsp=zero
      etermq=zero
c
      do 1 i=1,nthrds
c
c     find index into third neighbor list
c
      ist=indexl(i)
c
    3 continue
      if(ist.gt.idtn(isad,1)) then
      isad=isad+1
      goto 3
      endif
c
c     find local atoms involved
c
      isa=0
      jsa=0
      do 2 j=1,natoms
      if(isad.eq.igan(j)) isa=j
      if(idtn(ist,2).eq.igan(j)) jsa=j
    2 continue
c
c     find solute molecule numbers involved
c
      isf=imol(isa)
      jsf=imol(jsa)
c
      if(lssscl) then
c
      istt=iand(isss(isa),48)
      jstt=iand(isss(jsa),48)
      if(isf.ne.jsf) then
      if(istt.eq.16.or.jstt.eq.16) isrx=-1
      if(istt.eq.32.or.jstt.eq.32) isrx=1
      endif
c
      istt=iand(isss(isa),384)
      jstt=iand(isss(jsa),384)
      if(istt.eq.128.or.jstt.eq.128) isrx=-2
      if(istt.eq.256.or.jstt.eq.256) isrx=2
c
      istt=iand(isss(isa),384)
      jstt=iand(isss(jsa),384)
      if(istt.eq.128.and.jstt.eq.256) isrx=999
      if(istt.eq.256.and.jstt.eq.128) isrx=999
c
c      write(*,'(5i5)')
c     + isga(isfr+isa),isga(lssptr),istt,jstt,isrx(nax+jnum)
c
      endif
c
c     determine energy factor depending on atoms being dynamic
c     or fixed
c
      factu=zero
      if(iand(idyn(isa),mdynam).eq.ldynam.or.
     + iand(idyn(jsa),mdynam).eq.ldynam) factu=one
c      if((iand(idyn(isa),mdynam).eq.ldynam.and.
c     + iand(idyn(jsa),mdynam).ne.ldynam) .or.
c     + (iand(idyn(isa),mdynam).ne.ldynam.and.
c     + iand(idyn(jsa),mdynam).eq.ldynam) ) factu=half
      if(includ.eq.1) factu=one
c
c     evaluate the distance
c
      rxx=xs(isa,1)-xs(jsa,1)
      rxy=xs(isa,2)-xs(jsa,2)
      rxz=xs(isa,3)-xs(jsa,3)
c
c     periodic boundary conditions
c
      if(lpbc.or.lpbcs) then
      dx(1)=rxx
      dx(2)=rxy
      dx(3)=rxz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      rxx=dx(1)
      rxy=dx(2)
      rxz=dx(3)
      endif
c
      r2=rxx*rxx+rxy*rxy+rxz*rxz
      r2i=one/r2
c
      if(lssscl) then
      if(isrx.eq.999) then
      r2i=zero
      isrx=0
      endif
      if(isrx.gt.0) then
      r2i=one/(one/r2i+shift0(1))
      elseif(isrx.lt.0) then
      r2i=one/(one/r2i+shift1(1))
      endif
      endif
c
      r1i=sqrt(r2i)
      dfs=zero
      dfsu=zero
      r6i=r2i*r2i*r2i
c
c     get force constants
c
c      c6=cb6(iatt(isa),iatt(jsa),iset)
c      c12=cb12(iatt(isa),iatt(jsa),iset)
      if(nfhop.eq.0) then
      c6=vdw(iatt(isa),iatt(jsa),2,iset)
      c12=vdw(iatt(isa),iatt(jsa),4,iset)
      else
      c6=vdw(iatt(isa),iatt(jsa),2,lseq(isgm(jsa)))
      c12=vdw(iatt(isa),iatt(jsa),4,lseq(isgm(jsa)))
      endif
      cf6=six*c6
      cf12=twelve*c12
c
c     calculate the Lennard-Jones force contribution
c
      dfs=(cf12*r6i-cf6)*r6i*r2i
      dfsu=(c12*r6i-c6)*r6i
c
c     calculate the Lennard-Jones energies
c
      eterml=factu*dfsu
      ess(isf,jsf,7,1)=ess(isf,jsf,7,1)+eterml
      if(ipme.eq.0.or.isolvo.ne.0) then
c
c     calculate the electrostatic energies
c
c      q=q14fac*qs(isa,iset)*qs(jsa,iset)
      if(nfhop.eq.0) then
      q=q14fac*chg(ichg(isa),1,iset)*chg(ichg(jsa),1,iset)
      else
      q=q14fac*chg(ichg(isa),1,lseq(isgm(isa)))*
     + chg(ichg(jsa),1,lseq(isgm(jsa)))
      endif
      etermq=factu*q*r1i
      ess(isf,jsf,8,1)=ess(isf,jsf,8,1)+etermq
      if(npener.gt.0) then
      esa(igan(isa))=esa(igan(isa))+half*(eterml+etermq)
      esa(igan(jsa))=esa(igan(jsa))+half*(eterml+etermq)
      endif
c
c     calculate the electrostatic force contribution
c
      dfs=dfs+q*r1i*r2i
      else
      ferfc=erfc(ealpha/r1i)
      fderfc=ealpha*derfc(ealpha/r1i)
c
c     calculate the electrostatic energies
c
      if(nfhop.eq.0) then
      q=q14fac*chg(ichg(isa),1,iset)*chg(ichg(jsa),1,iset)
      else
      q=q14fac*chg(ichg(isa),1,lseq(isgm(isa)))*
     + chg(ichg(jsa),1,lseq(isgm(jsa)))
      endif
      qij=(one-q14fac)*q
      epmecs=epmecs-(one-ferfc)*qij*r1i*factu
      ess(isf,isf,8,1)=ess(isf,isf,8,1)-half*(one-ferfc)*q*r1i*factu
      ess(isf,jsf,8,1)=ess(jsf,jsf,8,1)-half*(one-ferfc)*q*r1i*factu
      etermq=ferfc*factu*q*r1i
      if(npener.gt.0) then
      esa(igan(isa))=esa(igan(isa))+half*(one-ferfc)*q*r1i*factu
      esa(igan(jsa))=esa(igan(jsa))+half*(one-ferfc)*q*r1i*factu
      endif
c
c     calculate the electrostatic force contribution
c
      dfsp=-(qij*r2i*((one-ferfc)*r1i+fderfc))
      dfs=dfs+q*r2i*(ferfc*r1i-fderfc)-dfsp
      endif
c
c     reaction field contributions
c
      if(ireact.ne.0) then
      ess(isf,jsf,8,1)=ess(isf,jsf,8,1)+factu*q*rffss/r2i
      if(npener.gt.0) then
      esa(igan(isa))=esa(igan(isa))+half*factu*q*rffss/r2i
      esa(igan(jsa))=esa(igan(jsa))+half*factu*q*rffss/r2i
      endif
      dfs=dfs-two*q*rffss
      endif
c
c     accumulate the forces
c
      fs(isa,1)=fs(isa,1)+dfs*rxx
      fs(jsa,1)=fs(jsa,1)-dfs*rxx
      fs(isa,2)=fs(isa,2)+dfs*rxy
      fs(jsa,2)=fs(jsa,2)-dfs*rxy
      fs(isa,3)=fs(isa,3)+dfs*rxz
      fs(jsa,3)=fs(jsa,3)-dfs*rxz
c
      if(ipme.ne.0) then
      vpmeb(1)=vpmeb(1)+dfsp*rxx*rxx
      vpmeb(2)=vpmeb(2)+dfsp*rxy*rxx
      vpmeb(3)=vpmeb(3)+dfsp*rxz*rxx
      vpmeb(4)=vpmeb(4)+dfsp*rxy*rxy
      vpmeb(5)=vpmeb(5)+dfsp*rxz*rxy
      vpmeb(6)=vpmeb(6)+dfsp*rxz*rxz
      endif
c
c     for thermodynamic integration evaluate the derivatives
c
      if(ithint) then
      if(ith(14)) then
c      c64=cb6(iatt(isa),iatt(jsa),4)
c      c124=cb12(iatt(isa),iatt(jsa),4)
      c64=vdw(iatt(isa),iatt(jsa),2,4)
      c124=vdw(iatt(isa),iatt(jsa),4,4)
      dercon=(c124*r6i-c64)*r6i
      if(isrx.gt.0) then
      c64=three*vdw(iatt(isa),iatt(jsa),1,iset)
      c124=six*vdw(iatt(isa),iatt(jsa),3,iset)
      dercon=dercon+shift0(4)*r2i*r6i*(c64-c124*r6i)
      elseif(isrx.lt.0) then
      c64=three*vdw(iatt(isa),iatt(jsa),1,iset)
      c124=six*vdw(iatt(isa),iatt(jsa),3,iset)
      dercon=dercon+shift1(4)*r2i*r6i*(c64-c124*r6i)
      endif
      deriv(15,1)=deriv(15,1)+dercon
      if(npgdec.gt.1) then
      dera(3,igan(isa))=dera(3,igan(isa))+half*dercon
      dera(3,igan(jsa))=dera(3,igan(jsa))+half*dercon
      endif
      endif
      if(ith(16)) then
c      qi=qs(isa,iset)
c      qi4=qs(isa,4)
c      qj=qs(jsa,iset)
c      qj4=qs(jsa,4)
      qi=chg(ichg(isa),1,iset)
      qi4=chg(ichg(isa),1,4)
      qj=chg(ichg(jsa),1,iset)
      qj4=chg(ichg(jsa),1,4)
      if(ipme.eq.0) then
      dercon=q14fac*(qi*qj4+qj*qi4)*r1i
      else
      dercon=q14fac*(qi*qj4+qj*qi4)*r1i
      endif
      if(isrx.gt.1) then
      dercon=dercon-half*shift0(4)*
     + chg(ichg(isa),1,iset)*chg(ichg(jsa),1,iset)*r2i
      elseif(isrx.lt.-1) then
      dercon=dercon-half*shift1(4)*
     + chg(ichg(isa),1,iset)*chg(ichg(jsa),1,iset)*r2i
      endif
      deriv(17,1)=deriv(17,1)+dercon
      if(ireact.ne.0) then
      deriv(17,1)=deriv(17,1)+q14fac*(qi*qj4+qj*qi4)*rffss/r2i
      endif
      if(npgdec.gt.1) then
      dera(4,igan(isa))=dera(4,igan(isa))+half*dercon
      dera(4,igan(jsa))=dera(4,igan(jsa))+half*dercon
      endif
      endif
      endif
c
c     thermodynamic perturbation 1
c
      if(ipert2) then
      if(ip2(14)) then
      ep2(1)=ep2(1)-eterml+factu*
     + (vdw(iatt(isa),iatt(jsa),4,2)*r6i-
     + vdw(iatt(isa),iatt(jsa),2,2))*r6i
      endif
      if(ip2(16).or.ip2(17)) then
      rxx=xs(isa,1)-xs(jsa,1)
      rxy=xs(isa,2)-xs(jsa,2)
      rxz=xs(isa,3)-xs(jsa,3)
      r2=rxx*rxx+rxy*rxy+rxz*rxz
      r2i=one/r2
      r1i=sqrt(r2i)
      if(ipme.eq.0) then
      q14=chg(ichg(isa),1,2)*chg(ichg(jsa),1,2)*q14fac
c      q14=qs(isa,2)*qs(jsa,2)*q14fac
      else
      q14=chg(ichg(isa),1,2)*chg(ichg(jsa),1,2)*q14fac*erfc(ealpha/r1i)
c      q14=qs(isa,2)*qs(jsa,2)*q14fac*erfc(ealpha/r1i)
      endif
      ep2(1)=ep2(1)-etermq+factu*q14*r1i
      if(ireact.ne.0) then
      ep2(1)=ep2(1)-etermq+factu*q14*rffss/r2i
      endif
      endif
      endif
c
c     thermodynamic perturbation 2
c
      if(ipert3) then
      if(ip3(14)) then
      ep3(1)=ep3(1)-eterml+factu*
     + (vdw(iatt(isa),iatt(jsa),4,3)*r6i-
     + vdw(iatt(isa),iatt(jsa),2,3))*r6i
      endif
      if(ip2(16).or.ip2(17)) then
      rxx=xs(isa,1)-xs(jsa,1)
      rxy=xs(isa,2)-xs(jsa,2)
      rxz=xs(isa,3)-xs(jsa,3)
      r2=rxx*rxx+rxy*rxy+rxz*rxz
      r2i=one/r2
      r1i=sqrt(r2i)
      if(ipme.eq.0) then
      q14=chg(ichg(isa),1,3)*chg(ichg(jsa),1,3)*q14fac
      else
      q14=chg(ichg(isa),1,3)*chg(ichg(jsa),1,3)*q14fac*erfc(ealpha/r1i)
      endif
      ep3(1)=ep3(1)-etermq+factu*q14*r1i
      if(ireact.ne.0) then
      ep3(1)=ep3(1)-etermq+factu*q14*rffss/r2i
      endif
      endif
      endif
    1 continue
c
      return
      end
