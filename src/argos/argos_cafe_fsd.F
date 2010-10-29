      subroutine argos_cafe_fsd(ndihes,indexl,msd,msp,idih,dih,rdih,
     + natoms,ndim,igan,isgm,imol,idyn,xs,fs,ess,lpbc,lpbcs,lupden,
     + lupdti,dera,lseq)
c
c $Id$
c
c     argos_cafe_fsd returns forces and energies for solute angles
c
c     =======================================================
c
c     description of arguments
c     ------------------------
c
c     in: integer ndihes     = number of angles to consider
c         integer indexl     = index list
c
c         integer idsd(msd)  = global atom id i
c         integer jdsd(msd)  = global atom id j
c         integer kdsd(msd)  = global atom id k
c         integer ldsd(msd)  = global atom id l
c         real*8 cdsd(msd,6) = dihedral angle force constants
c         real*8 ddsd(msd,6) = dihedral angle reference value
c         real*8 edsd(msd,6) = dihedral angle multiplicity
c         real*8 rdsd(msd)   = dihedral angle value
c
c         integer natoms     = number of atoms in arrays
c         integer ndim       = leading dimension atom arrays
c         integer igan(ndim) = global atom numbers
c         integer imol(ndim) = atom molecule fraction
c         integer idyn(ndim) = atom dynamics type
c         real*8 xs(ndim,3)  = atom coordinates
c
c         logical lupden     = if .true. energies are updated
c
c    out: real*8 fs(ndim,3)  = atom forces (ACCUMULATED)
c         real*8 usd(msd)    = dihedral angle energies
c
      implicit none
c
#include "argos_cafe_common.fh"
c
c     declaration of arguments
c     ------------------------
c
      integer msd,msp
      integer idih(msd,5)
      real*8 dih(msd,msp,mset),rdih(msd,2)
      integer isgm(msa),lseq(mseq)
c
      integer ndihes
      integer indexl(ndihes)
c
      integer natoms,ndim
      integer igan(ndim),imol(ndim),idyn(ndim)
      real*8 xs(ndim,3),fs(ndim,3),ess(msf,msf,mpe,2)
c
      logical lpbc,lpbcs,lupden,lupdti
      real*8 dera(6,nsatot)
c
c     declaration of local variables
c     ------------------------------
c
      integer i,j,isd,isa,jsa,ksa,lsa,isf,jsf,ksf,lsf
      integer ifacu
      real*8 ang,for,dfor,phi,cphi,cphii,sphi,sphii
      real*8 rmul,factu,dercon
      real*8 xsijx,xskjx,xsijy,xskjy,xsijz,xskjz
      real*8 xsklx,xsjlx,xskly,xsjly,xsklz,xsjlz
      real*8 xsikx,xsiky,xsikz,xmx,xmy,xmz,xnx,xny,xnz,xdx,xdy,xdz
      real*8 xex,xey,xez,xox,xoy,xoz,xpx,xpy,xpz
      real*8 dfsix,dfsiy,dfsiz,dfsjx,dfsjy,dfsjz
      real*8 dfskx,dfsky,dfskz,dfslx,dfsly,dfslz
      real*8 rm2i,rn2i,rmni,s,rpa,dx(3),eterm
      real*8 rkj,rkjo,rkjp
c
c#include "argos_cafe_funcs_dec.fh"
#include "bitops_decls.fh"
c#include "argos_cafe_funcs_sfn.fh"
#include "bitops_funcs.fh"
c
      do 1 i=1,ndihes
c
c     find index into list of dihedrals
c     ---------------------------------
c
      isd=indexl(i)
c
c     find local atom numbers involved in dihedral
c     --------------------------------------------
c
      isa=0
      jsa=0
      ksa=0
      lsa=0
c
      do 2 j=1,natoms
      if(idih(isd,1).eq.igan(j)) isa=j
      if(idih(isd,2).eq.igan(j)) jsa=j
      if(idih(isd,3).eq.igan(j)) ksa=j
      if(idih(isd,4).eq.igan(j)) lsa=j
    2 continue
c
c     find solute molecule numbers involved in dihedral
c     -------------------------------------------------
c
      isf=imol(isa)
      jsf=imol(jsa)
      ksf=imol(ksa)
      lsf=imol(lsa)
c
c     find energy factor that depends on the atoms involved
c     being dynamic or fixed
c     -----------------------------------------------------
c
      ifacu=0
      if(iand(idyn(isa),mdynam).eq.ldynam) ifacu=ifacu+1
      if(iand(idyn(jsa),mdynam).eq.ldynam) ifacu=ifacu+1
      if(iand(idyn(ksa),mdynam).eq.ldynam) ifacu=ifacu+1
      if(iand(idyn(lsa),mdynam).eq.ldynam) ifacu=ifacu+1
c      factu=dble(ifacu)/four
      factu=one
      if(ifacu.eq.0) factu=zero
      if(includ.eq.1) factu=one
c
c     find reference angle and force constants
c     ----------------------------------------
c
      if(nfhop.eq.0) then
      ang=dih(isd,2,iset)
      for=dih(isd,3,iset)
      rmul=dih(isd,1,iset)
      else
      ang=dih(isd,2,lseq(isgm(jsa)))
      for=dih(isd,3,lseq(isgm(jsa)))
      rmul=dih(isd,1,lseq(isgm(jsa)))
      endif
c
c     determine the dihedral angle
c     ----------------------------
c
      xsijx=xs(isa,1)-xs(jsa,1)
      xskjx=xs(ksa,1)-xs(jsa,1)
      xsklx=xs(ksa,1)-xs(lsa,1)
      xsikx=xsijx-xskjx
      xsjlx=xsklx-xskjx
      xsijy=xs(isa,2)-xs(jsa,2)
      xskjy=xs(ksa,2)-xs(jsa,2)
      xskly=xs(ksa,2)-xs(lsa,2)
      xsiky=xsijy-xskjy
      xsjly=xskly-xskjy
      xsijz=xs(isa,3)-xs(jsa,3)
      xskjz=xs(ksa,3)-xs(jsa,3)
      xsklz=xs(ksa,3)-xs(lsa,3)
      xsikz=xsijz-xskjz
      xsjlz=xsklz-xskjz
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
      dx(1)=xsikx
      dx(2)=xsiky
      dx(3)=xsikz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xsikx=dx(1)
      xsiky=dx(2)
      xsikz=dx(3)
      dx(1)=xskjx
      dx(2)=xskjy
      dx(3)=xskjz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xskjx=dx(1)
      xskjy=dx(2)
      xskjz=dx(3)
      dx(1)=xsklx
      dx(2)=xskly
      dx(3)=xsklz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xsklx=dx(1)
      xskly=dx(2)
      xsklz=dx(3)
      dx(1)=xsjlx
      dx(2)=xsjly
      dx(3)=xsjlz
      call argos_cafe_pbc(1,dx,1,dx,1,0,1,1)
      xsjlx=dx(1)
      xsjly=dx(2)
      xsjlz=dx(3)
      endif
c
      xmx=xsijy*xskjz-xskjy*xsijz
      xmy=xsijz*xskjx-xskjz*xsijx
      xmz=xsijx*xskjy-xskjx*xsijy
      xnx=xskjy*xsklz-xskly*xskjz
      xny=xskjz*xsklx-xsklz*xskjx
      xnz=xskjx*xskly-xsklx*xskjy
      rm2i=one/(xmx*xmx+xmy*xmy+xmz*xmz)
      rn2i=one/(xnx*xnx+xny*xny+xnz*xnz)
      rmni=sqrt(rm2i*rn2i)
      cphi=(xmx*xnx+xmy*xny+xmz*xnz)*rmni
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      phi=acos(cphi)
      s=xskjx*(xmy*xnz-xmz*xny) +xskjy*(xmz*xnx-xmx*xnz)
     + +xskjz*(xmx*xny-xmy*xnx)
      if(s.lt.zero) phi=-phi
      rdih(isd,1)=phi
      sphi=sin(phi)
      rpa=rmul*phi-ang
c
c     evaluate energies
c     -----------------
c
      rdih(isd,2)=for*(one+cos(rpa))
c
      eterm=zero
      if(lupden) then
      ess(isf,isf,3,1)=ess(isf,isf,3,1)+fourth*factu*rdih(isd,2)
      ess(jsf,jsf,3,1)=ess(jsf,jsf,3,1)+fourth*factu*rdih(isd,2)
      ess(ksf,ksf,3,1)=ess(ksf,ksf,3,1)+fourth*factu*rdih(isd,2)
      ess(lsf,lsf,3,1)=ess(lsf,lsf,3,1)+fourth*factu*rdih(isd,2)
      endif
      eterm=factu*rdih(isd,2)
c
      dfor=(-for)*rmul*sin(rpa)
c
c     for thermodynamic perturbations evaluate the energies using
c     the 'perturbed' parameters in set 2 and/or set 3
c     -----------------------------------------------------------
c
      if(ip2(21)) ep2(1)=ep2(1)-eterm+
     + factu*dih(isd,3,2)*(one+cos(dih(isd,1,2)*phi-dih(isd,2,2)))
      if(ip3(21)) ep3(1)=ep3(1)-eterm+
     + factu*dih(isd,3,3)*(one+cos(dih(isd,1,3)*phi-dih(isd,2,3)))
c
c     for normal angles use cosine equation to get the forces
c     -------------------------------------------------------
c
      if(abs(sphi).gt.small) then
      sphii=one/sphi
      xdx=(-dfor)*sphii*(rmni*xnx-cphi*rm2i*xmx)
      xex=(-dfor)*sphii*(rmni*xmx-cphi*rn2i*xnx)
      xdy=(-dfor)*sphii*(rmni*xny-cphi*rm2i*xmy)
      xey=(-dfor)*sphii*(rmni*xmy-cphi*rn2i*xny)
      xdz=(-dfor)*sphii*(rmni*xnz-cphi*rm2i*xmz)
      xez=(-dfor)*sphii*(rmni*xmz-cphi*rn2i*xnz)
      dfsix=xskjy*xdz-xskjz*xdy
      dfsiy=xskjz*xdx-xskjx*xdz
      dfsiz=xskjx*xdy-xskjy*xdx
      dfsjx=xsiky*xdz-xsikz*xdy-xskly*xez+xsklz*xey
      dfsjy=xsikz*xdx-xsikx*xdz-xsklz*xex+xsklx*xez
      dfsjz=xsikx*xdy-xsiky*xdx-xsklx*xey+xskly*xex
      dfskx=xsjly*xez-xsjlz*xey-xsijy*xdz+xsijz*xdy
      dfsky=xsjlz*xex-xsjlx*xez-xsijz*xdx+xsijx*xdz
      dfskz=xsjlx*xey-xsjly*xex-xsijx*xdy+xsijy*xdx
      dfslx=xskjy*xez-xskjz*xey
      dfsly=xskjz*xex-xskjx*xez
      dfslz=xskjx*xey-xskjy*xex
      else
c
c     for small angles use the sine equations to get forces
c     -----------------------------------------------------
c
      cphii=one/cphi
      rkj=sqrt(xskjx*xskjx+xskjy*xskjy+xskjz*xskjz)
      xox=xskly*xsikz-xsklz*xsiky
      xoy=xsklz*xsikx-xsklx*xsikz
      xoz=xsklx*xsiky-xskly*xsikx
      rkjo=(xskjx*xox+xskjy*xoy+xskjz*xoz)/rkj
      xpx=xsijy*xsjlz-xsijz*xsjly
      xpy=xsijz*xsjlx-xsijx*xsjlz
      xpz=xsijx*xsjly-xsijy*xsjlx
      rkjp=(xskjx*xpx+xskjy*xpy+xskjz*xpz)/rkj
      dfsix=dfor*cphii* (rkj*rmni*xnx-sphi*rm2i*(xskjy*xmz-xskjz*xmy))
      dfsiy=dfor*cphii* (rkj*rmni*xny-sphi*rm2i*(xskjz*xmx-xskjx*xmz))
      dfsiz=dfor*cphii* (rkj*rmni*xnz-sphi*rm2i*(xskjx*xmy-xskjy*xmx))
      dfsjx=dfor*cphii*((-rmni)*(rkjo*xskjx+rkj*xox)-
     + sphi*rm2i*(xsiky*xmz-xsikz*xmy)+
     + sphi*rn2i*(xskly*xnz-xsklz*xny))
      dfsjy=dfor*cphii*((-rmni)*(rkjo*xskjy+rkj*xoy)-
     + sphi*rm2i*(xsikz*xmx-xsikx*xmz)+
     + sphi*rn2i*(xsklz*xnx-xsklx*xnz))
      dfsjz=dfor*cphii*((-rmni)*(rkjo*xskjz+rkj*xoz)-
     + sphi*rm2i*(xsikx*xmy-xsiky*xmx)+
     + sphi*rn2i*(xsklx*xny-xskly*xnx))
      dfskx=dfor*cphii*((-rmni)*(rkjp*xskjx+rkj*xpx)-
     + sphi*rn2i*(xsjly*xnz-xsjlz*xny)+
     + sphi*rm2i*(xsijy*xmz-xsijz*xmy))
      dfsky=dfor*cphii*((-rmni)*(rkjp*xskjy+rkj*xpy)-
     + sphi*rn2i*(xsjlz*xnx-xsjlx*xnz)+
     + sphi*rm2i*(xsijz*xmx-xsijx*xmz))
      dfskz=dfor*cphii*((-rmni)*(rkjp*xskjz+rkj*xpz)-
     + sphi*rn2i*(xsjlx*xny-xsjly*xnx)+
     + sphi*rm2i*(xsijx*xmy-xsijy*xmx))
      dfslx=dfor*cphii*
     + ((-rkj)*rmni*xmx-sphi*rn2i*(xskjy*xnz-xskjz*xny))
      dfsly=dfor*cphii*
     + ((-rkj)*rmni*xmy-sphi*rn2i*(xskjz*xnx-xskjx*xnz))
      dfslz=dfor*cphii*
     + ((-rkj)*rmni*xmz-sphi*rn2i*(xskjx*xny-xskjy*xnx))
      endif
c
c     accumulate the forces
c     ---------------------
c
      fs(isa,1)=fs(isa,1)-dfsix
      fs(jsa,1)=fs(jsa,1)-dfsjx
      fs(ksa,1)=fs(ksa,1)-dfskx
      fs(lsa,1)=fs(lsa,1)-dfslx
      fs(isa,2)=fs(isa,2)-dfsiy
      fs(jsa,2)=fs(jsa,2)-dfsjy
      fs(ksa,2)=fs(ksa,2)-dfsky
      fs(lsa,2)=fs(lsa,2)-dfsly
      fs(isa,3)=fs(isa,3)-dfsiz
      fs(jsa,3)=fs(jsa,3)-dfsjz
      fs(ksa,3)=fs(ksa,3)-dfskz
      fs(lsa,3)=fs(lsa,3)-dfslz
c
      if(lupdti) then
c
c     for thermodynamic integration evaluate the derivative
c     -----------------------------------------------------
c
      if(ith(21)) then
      dercon=(one+cos(rpa))*dih(isd,3,4)
     + -for*sin(rpa)*(phi*dih(isd,1,4)-dih(isd,2,4))
      deriv(21,1)=deriv(21,1)+dercon
      if(npgdec.gt.1) then
      dera(5,idih(isd,1))=dera(5,idih(isd,1))+fourth*dercon
      dera(5,idih(isd,2))=dera(5,idih(isd,2))+fourth*dercon
      dera(5,idih(isd,3))=dera(5,idih(isd,3))+fourth*dercon
      dera(5,idih(isd,4))=dera(5,idih(isd,4))+fourth*dercon
      endif
      endif
      endif
c
    1 continue
c
      return
      end
