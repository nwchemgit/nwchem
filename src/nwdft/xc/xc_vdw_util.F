c
c     $Id$
c

c     Interpolation of C6 to account for coordination number dependence
c
C> \brief Interpolate \f$C_6\f$ to account for coordination numbers
C>
C> Interpolate the \f$C_6\f$ coefficients according to Eq.(16) of [1]:
C> \f{eqnarray*}{
C>   C_6^{AB}\left(\mathrm{CN}^A(R),\mathrm{CN}^B(R)\right)
C>   &=& \frac{\sum_i^{N_A}\sum_j^{N_B}C_{6,\mathrm{ref}}^{AB}\left(\mathrm{CN}_i^A,\mathrm{CN}_j^B\right)L_{ij}^{AB}(R)}{\sum_i^{N_A}\sum_j^{N_B}L_{ij}^{AB}(R)} \\\\
C>   L_{ij}^{AB}(R) &=& e^{-k_3[(\mathrm{CN}^A(R)-\mathrm{CN}^A_i)^2
C>                             +(\mathrm{CN}^B(R)-\mathrm{CN}^B_j)^2]}
C> \f}
C> where \f$\mathrm{CN}^A(R)\f$ and \f$\mathrm{CN}^B(R)\f$ are geometry
C> dependent coordination numbers for atoms \f$A\f$ and \f$B\f$.
C> The quantities \f$\mathrm{CN}^A_i\f$ and \f$\mathrm{CN}^B_j\f$
C> are coordination numbers of atoms of the same type as \f$A\f$ and
C> \f$B\f$ in reference molecules \f$i\f$ and \f$j\f$. These latter 
C> quantities as well as \f$C_{6,\mathrm{ref}}^{AB}\f$ are precomputed
C> constants, independent of the current geometry.
C> 
C> \return The coordination number dependent \f$C_6\f$ coefficient
C>
C> ### References ###
C>
C> [1] S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
C>     "A consistent and accurate ab initio parametrization of density
C>      functional dispersion correction (DFT-D) for the 94 elements
C>      H-Pu",
C>     J. Chem. Phys. (2010) <b>132</b>, 154104, DOI:
C>     <a href="https://doi.org/10.1063/1.3382344"> 
C>     10.1063/1.3382344</a>.
C>
      double precision function c6cn(iat,jat,cni,cnj)
      implicit none
#include "xc_vdw.fh"
      integer iat !< [Input] The atomic number of \f$A\f$ (i.e. \f$Z_A\f$)
      integer jat !< [Input] The atomic number of \f$B\f$ (i.e. \f$Z_B\f$)
      double precision cni !< [Input] \f$\mathrm{CN}^A(R)\f$
      double precision cnj !< [Input] \f$\mathrm{CN}^B(R)\f$
      double precision top,bottom,dist,c6_ref
      double precision cna,cnb
      integer i,j
      double precision eps
      parameter (eps=1d-30)
c
      c6cn=0.0d0
      top=0.0d0
      bottom=0.0d0
      do i=1,maxcn(iat)
         do j=1,maxcn(jat)
            c6_ref=c6ab(iat,jat,i,j,1)
            if (c6_ref.gt.0.0d0) then
               cna=c6ab(iat,jat,i,j,2)
               cnb=c6ab(iat,jat,i,j,3)
               dist=(cna-cni)**2+(cnb-cnj)**2
               top=top+dexp(k3*dist)*c6_ref
               bottom=bottom+dexp(k3*dist)
            endif
         enddo
      enddo
c
      if (bottom.gt.eps) then
         c6cn=top/bottom
      else
         c6cn=0.0d0
      endif
      return
      end
c
c     Analytical gradient of coordination number dependent C6
c
c  Comput derivatives of the C6 coefficients with respect to
c  the coordination numbers.
c
c  The nuclear gradient will be assembled on a later stage
c
C> 
C> ### References ###
C>
C> [1] S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
C>     "A consistent and accurate ab initio parametrization of density
C>      functional dispersion correction (DFT-D) for the 94 elements
C>      H-Pu",
C>     J. Chem. Phys. (2010) <b>132</b>, 154104, DOI:
C>     <a href="https://doi.org/10.1063/1.3382344"> 
C>     10.1063/1.3382344</a>.
C>
      subroutine c6_grad(iat,jat,x,z,n,dc6i,dc6j,cni,cnj)
      implicit none
#include "xc_vdw.fh"
      integer n   !< [Input] The number of atoms
      integer iat !< [Input] Atom \f$A\f$
      integer jat !< [Input] Atom \f$B\f$
      integer z(n) !< [Input] The atomic numbers of the atoms
      double precision x(3,n)
      double precision dc6i, dc6j
      double precision cni,cnj
      double precision crd_nr
      double precision c6ref
      external crd_nr
      double precision t1,t2,nomi,nomj,denomi,denomj
      double precision tmp1,tmp2,tmp3,tmp4,fac1,fac2
      double precision eps
      parameter (eps=1d-30)
      integer i,j
c
      t1=0.0d0
      t2=0.0d0
      nomi = 0d0
      nomj = 0d0
      denomi = 0d0
      denomj = 0d0

      do i=1,maxcn(z(iat))
        do j=1,maxcn(z(jat))
          c6ref = c6AB(z(iat),z(jat),i,j,1)
          if (c6ref.gt.0d0) then
            tmp1=c6AB(z(iat),z(jat),i,j,3)-cnj
            tmp2=c6AB(z(iat),z(jat),i,j,2)-cni
            tmp3=dexp(k3*(tmp1*tmp1+tmp2*tmp2))
            t1=t1+c6AB(z(iat),z(jat),i,j,1)*tmp3
            t2=t2+tmp3
            fac1=tmp3*k3*2.0d0
            fac2=fac1*c6ref
            nomi = nomi - fac2*tmp2
            denomi = denomi - fac1*tmp2
            nomj = nomj - fac2*tmp1
            denomj = denomj -  fac1*tmp1
          endif
        enddo
      enddo

      if (t2.gt.eps) then
        dc6i = ((nomi*t2 - denomi*t1))/(t2*t2)
        dc6j = ((nomj*t2 - denomj*t1))/(t2*t2)
      else
        dc6i = 0d0
        dc6j = 0d0
      endif

      end
c
      double precision function c6ij_sk(i,j,z)
      implicit none
#include "xc_vdw.fh"
      integer i,j
      integer z(*)
c     
      integer zi,zj
      double precision ni,nj
      external xc_neff
      double precision xc_neff,ci,cj
c
      c6ij_sk=0.0d0
      zi=z(i)
      zj=z(j)
      if (zi.le.0.or.zj.le.0) then 
         c6ij_sk=0.0d0
         return
      end if
c
      ci=c6(zi)
      cj=c6(zj)
c
c     write(6,*) zi,zj,ci,cj
c
      if(ci.le.0.d0.or.cj.le.0.d0) then
         c6ij_sk=0.0d0
         return
      endif
c
      if(ivdw.eq.1)then
c version 1 
c     S. Grimme J Comp Chem 25, 1463 (2004)
          c6ij_sk =2.d0*ci*cj/(ci+cj)
      else if(ivdw.eq.2) then
c version 2
c     S. Grimme J Comp Chem 27, 1787 (2006)
          c6ij_sk =sqrt(ci*cj)
      else if(ivdw.eq.3) then
c version 3
c     S. Grimme J Chem Phys 132, 154104 (2010)
      else if(ivdw.eq.0) then
          ni=xc_neff(zi)
          nj=xc_neff(zj)
          c6ij_sk=2d0*(ci*ci*cj*cj*ni*nj)**(1d0/3d0)/
     /       ((nj*nj*ci)**(1d0/3d0)+(ni*ni*cj)**(1d0/3d0))
      endif
      return
      end
c
c     Derivatives of coordination number function
c
C> \brief Fractional coordination number gradients
C> 
C> This routine evaluates the gradients of the fractional coordination
C> numbers of the atoms. The fractional coordination numbers are
C> given by Eq.(15) of [1]:
C> \f{eqnarray*}{
C>    \mathrm{CN}^A &=& \sum_{B\ne A}^{N_{atom}} 
C>    \frac{1}{1+e^{-k_1\left(k_2\frac{R_{A,\mathrm{cov}}+R_{B,\mathrm{cov}}}{R_{AB}}-1\right)}}
C> \f}
C> In the actual implementation the factor \f$k_2\f$ has been absorbed
C> in \f$R_{A,\mathrm{cov}}\f$ and those scaled covalent radii are
C> stored in `cov_table`.
C>
C> The derivatives wrt. atoms \f$B\f$ are calculated separately and
C> stored for every atom \f$A\f$. The derivatives wrt. atom \f$A\f$ are
C> calculated exploiting the translational invariance of the gradient.
C>
C> ### References ###
C>
C> [1] S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
C>     "A consistent and accurate ab initio parametrization of density
C>      functional dispersion correction (DFT-D) for the 94 elements
C>      H-Pu",
C>     J. Chem. Phys. (2010) <b>132</b>, 154104, DOI:
C>     <a href="https://doi.org/10.1063/1.3382344">
C>     10.1063/1.3382344</a>.
C>
      subroutine crd_nr_der(n,xyz,iz,cnij,cnijk)
      implicit none
#include "xc_vdw.fh"
      integer n      !< [Input] The number of atoms \f$N_{atom}\f$
      integer iz(n)  !< [Input] The atomic numbers of the atoms
      double precision xyz(3,n) !< [Input] The atomic coordinates
      double precision cnij(3,n) !< [Output] \f$\mathrm{cnij}(i,A)=\frac{\partial\mathrm{CN}^A}{\partial \left(\vec{R}_A\right)_i}\f$
      double precision cnijk(3,n,n) !< [Output] \f$\mathrm{cnijk}(i,B,A)=\frac{\partial\mathrm{CN}^A}{\partial \left(\vec{R}_B\right)_i}\f$
      double precision dx,dy,dz,r,cov_rad,expf,fac2,fac3
      integer iat !< Atom B
      integer i   !< Atom A
c
      do i=1,n
         cnij(1,i)=0.0d0
         cnij(2,i)=0.0d0
         cnij(3,i)=0.0d0
         if (iz(i).ne.0) then
            do iat=1,n
               if (i.ne.iat.and.iz(iat).ne.0) then
                  dx=xyz(1,iat)-xyz(1,i)
                  dy=xyz(2,iat)-xyz(2,i)
                  dz=xyz(3,iat)-xyz(3,i)
                  r=dsqrt(dx*dx+dy*dy+dz*dz)
                  cov_rad=cov_table(iz(iat))+cov_table(iz(i))
                  expf=dexp(-k1*((cov_rad/r)-1.0d0))
                  fac2=1.0d0/(expf+1.0d0)
                  fac3=k1*cov_rad*expf*fac2*fac2/(r*r*r) 
                  cnij(1,i)=cnij(1,i)-fac3*dx
                  cnij(2,i)=cnij(2,i)-fac3*dy
                  cnij(3,i)=cnij(3,i)-fac3*dz
                  cnijk(1,iat,i)=fac3*dx
                  cnijk(2,iat,i)=fac3*dy
                  cnijk(3,iat,i)=fac3*dz
               else
                  cnijk(1,iat,i)=0.0d0
                  cnijk(2,iat,i)=0.0d0
                  cnijk(3,iat,i)=0.0d0
               endif
            enddo
         else
            do iat=1,n
               cnijk(1,iat,i)=0.0d0
               cnijk(2,iat,i)=0.0d0
               cnijk(3,iat,i)=0.0d0
            enddo
         endif
      enddo
      return
      end

      double precision function dcndrij(iat,jat,r,iz)
      implicit none
#include "xc_vdw.fh"
      integer iat,jat
      double precision r
      double precision rcov,expterm
      integer iz(*)

      rcov = cov_table(iz(iat))+cov_table(iz(jat))
      expterm = exp(-k1*(rcov/r-1d0))
      dcndrij = -k1*rcov*expterm/(r*(expterm+1d0))**2

      end
c
c     Coordination numbers based on inverse damping function 
c
C> \brief Fractional coordination numbers
C> 
C> This routine evaluates the fractional coordination of an atom
C> according to Eq.(15) of [1]:
C> \f{eqnarray*}{
C>    \mathrm{CN}^A &=& \sum_{B\ne A}^{N_{atom}} 
C>    \frac{1}{1+e^{-k_1\left(k_2\frac{R_{A,\mathrm{cov}}+R_{B,\mathrm{cov}}}{R_{AB}}-1\right)}}
C> \f}
C> In the actual implementation the factor \f$k_2\f$ has been absorbed
C> in \f$R_{A,\mathrm{cov}}\f$ and those scaled covalent radii are
C> stored in `cov_table`.
C>
C> \return The value of \f$\mathrm{CN}^A\f$.
C>
C> ### References ###
C>
C> [1] S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
C>     "A consistent and accurate ab initio parametrization of density
C>      functional dispersion correction (DFT-D) for the 94 elements
C>      H-Pu",
C>     J. Chem. Phys. (2010) <b>132</b>, 154104, DOI:
C>     <a href="https://doi.org/10.1063/1.3382344">
C>     10.1063/1.3382344</a>.
C>
      double precision function crd_nr(iat,nat,xyz,iz)
      implicit none
#include "xc_vdw.fh"
      integer iat    !< [Input] The atom \f$A\f$
      integer nat    !< [Input] The number of atoms \f$N_{atom}\f$
      integer iz(nat) !< [Input] The atomic numbers of the atoms
      double precision xyz(3,nat) !< The atomic coordinates
      double precision dx,dy,dz,r,rcov
      integer i
c       
      crd_nr=0.0d0
      if (iz(iat).eq.0) return
      do i=1,nat
         if(iat.ne.i.and.iz(i).ne.0)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=dsqrt(dx*dx+dy*dy+dz*dz)
            rcov=cov_table(iz(i))+cov_table(iz(iat))
            crd_nr=crd_nr+1.d0/(1.d0+dexp(-k1*(rcov/r-1.0d0)))
         endif
      enddo
      return
      end
c
c      f1dmp(x1,x2,ff)=(alpha/x2)*ff*ff*exp(-alpha*((x1/x2)-1.0d0))
c
      double precision function f1dmp(x1,x2,ff)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2,ff
      f1dmp=(alpha/x2)*ff*ff*dexp(-alpha*((x1/x2)-1.0d0))
      return
      end
c
c      fdmp3(x1,x2)=1.0d0/(1.0d0+6.0d0*(x1/x2)**-a)
c
      double precision function fdmp3(x1,x2,a)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2,a
      fdmp3=1.0d0/(1.0d0+6.0d0*((x1/x2)**(-a)))
      return
      end
C>
C> \brief Evaluate the Becke-Johnson damping factor
C>
C> Evaluate the Becke-Johnson damping factor [1,2] for the dispersion
C> correction:
C> \f{eqnarray}{
C>   \mathrm{xc\_fdmpbj}(R,R^0,a_1,a_2,n) &=&
C>   \frac{1}{R^n+f(R^0,a_1,a_2)^n} \\\\
C>   f(R^0,a_1,a_2) &=& a_1 R^0 + a_2
C> \f}
C>
C> ### References ###
C>
C> [1] A.D. Becke, E.R. Johnson,
C>     "A unified density-functional treatment of dynamical,
C>      nondynamical and dispersion correlations",
C>     J. Chem. Phys. (2007) <b>127</b> 124108, DOI:
C>     <a href="https://doi.org/10.1063/1.2768530">
C>     10.1063/1.2768530</a> (See appendix C).
C>
C> [2] S. Grimme, S. Ehrlich, L. Goerigk,
C>     "Effect of the damping function in dispersion corrected
C>      density functional theory", J. Comput. Chem. (2011)
C>     <b>32</b>, pp. 1456-1465, DOI:
C>     <a href="https://doi.org/10.1002/jcc.21759">
C>     10.1002/jcc.21759</a> (See Eqs.(5-6)).
C>
      double precision function xc_fdmpbj(r,r0,a1,a2,n)
c
      implicit none
c
      integer          n  !< [Input] The order of the term
      double precision r  !< [Input] Interactomic distance \f$R\f$
      double precision r0 !< [Input] Parameter \f$R^0=\sqrt{C_8/C_6}\f$
      double precision a1 !< [Input] Parameter \f$a_1\f$
      double precision a2 !< [Input] Parameter \f$a_2\f$
c
      double precision f
      f=a1*r0+a2
      xc_fdmpbj=1.0d0/(r**n+f**n)
      return
      end
C>
C> \brief Evaluate the gradient of the Becke-Johnson damping factor
C>
C> The Becke-Johnson damping factor [1,2] for the dispersion
C> correction is:
C> \f{eqnarray}{
C>   \mathrm{xc\_fdmpbj}(R,R^0,a_1,a_2,n) &=&
C>   \frac{1}{R^n+f(R^0,a_1,a_2)^n} \\\\
C>   f(R^0,a_1,a_2) &=& a_1 R^0 + a_2
C> \f}
C> This function evaluates the gradient of the Becke-Johnson damping
C> wrt. \f$R\f$:
C> \f{eqnarray}{
C>   \frac{\partial\mathrm{xc\_fdmpbj}(R,R^0,a_1,a_2,n)}{\partial R} &=&
C>   \frac{-n R^{n-1}}{\mathrm{xc\_fdmpbj}(R,R^0,a_1,a_2,n)^2}
C> \f}
C>
C> ### References ###
C>
C> [1] A.D. Becke, E.R. Johnson,
C>     "A unified density-functional treatment of dynamical,
C>      nondynamical and dispersion correlations",
C>     J. Chem. Phys. (2007) <b>127</b> 124108, DOI:
C>     <a href="https://doi.org/10.1063/1.2768530">
C>     10.1063/1.2768530</a> (See appendix C).
C>
C> [2] S. Grimme, S. Ehrlich, L. Goerigk,
C>     "Effect of the damping function in dispersion corrected
C>      density functional theory", J. Comput. Chem. (2011)
C>     <b>32</b>, pp. 1456-1465, DOI:
C>     <a href="https://doi.org/10.1002/jcc.21759">
C>     10.1002/jcc.21759</a> (See Eqs.(5-6)).
C>
      double precision function xc_fdmpbj_d1(r,r0,a1,a2,n)
c
      implicit none
c
      integer          n  !< [Input] The order of the term
      double precision r  !< [Input] Interactomic distance \f$R\f$
      double precision r0 !< [Input] Parameter \f$R^0=\sqrt{C_8/C_6}\f$
      double precision a1 !< [Input] Parameter \f$a_1\f$
      double precision a2 !< [Input] Parameter \f$a_2\f$
c
      double precision xc_fdmpbj
      external         xc_fdmpbj
c
      xc_fdmpbj_d1=-n*r**(n-1)*(xc_fdmpbj(r,r0,a1,a2,n)**2)
c
      return
      end
c
c      fdmp(x1,x2)=1.0d0/(1.0d0+exp(-alpha*((x1/x2)-1.0d0)))
c
      double precision function fdmp(x1,x2)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2
      fdmp=1.0d0/(1.0d0+dexp(-alpha*((x1/x2)-1.0d0)))
      return
      end
C>
C> \brief Fetch the functional dependent scaling factors
C>
C> Different scaling factors are used for the different dispersion
C> model. These scaling factors are also functional dependent. 
C> This routine looks up the model and functional dependent
C> factors. The models support are DFT-D1, DFT-D2, DFT-D3 and
C> DFT-D3BJ [1] (see Table 1).
C>
C> ### References ###
C>
C> [1] A. Hoffmann, R. Grunzke, S. Herres-Pawlis,
C>     "Insights into the influence of dispersion correction in the
C>      theoretical treatment of Guanidine-Quinoline Copper(I)
C>      complexes", J. Comput. Chem. (2014) <b>35</b>, pp. 1943-1950,
C>      DOI:
C>      <a href="https://doi.org/10.1002/jcc.23706">
C>      10.1002/jcc.23706</a>
C>
      subroutine get_scaling_fac(scales6,scales8,scalesr6,scalesr8,
     +                           a1,a2)
      implicit none
#include "cdft.fh"
#include "xc_vdw.fh"
#include "stdio.fh"
#include "global.fh"
      double precision scales6,scales8,scalesr6,scalesr8,a1,a2
      logical read_ok
      read_ok=.false.
c
c DFT-D1 version
c
      if( ivdw.eq.1) then
        scales8  = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
          scales6 = 1.40d0
c PBE / PBE96
        else if(xccomb(58).or.(xfac(7).eq.1d0.and.cfac(12).eq.1d0)) then
         scales6 = 0.70d0
c BP86
         else if(xfac(3).eq.1.0d0.and.cfac(4).eq.1.d0) then
         scales6 = 1.30d0
c default
        else 
         scales6 = 1.00d0
         if (ga_nodeid().eq.0) then
           write(luout,*)'WARNING:'
           write(luout,*)'you are using HF or an untested DF together'
           write(luout,*)'with the dispersion correction!'
         endif
        endif
       else if ( ivdw.eq.2) then
c
c DFT-D2 version
c
        scales8  = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scales6 = 1.20d0
c PBE / PBE96
        else if(xccomb(58).or.(xfac(7).eq.1d0.and.cfac(12).eq.1d0)) then
           scales6 = 0.75d0
c TPSS
        else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scales6 = 1.00d0
c B3LYP
        else if (xccomb(2)) then
           scales6 = 1.05d0
c B97-D
        else if (xccomb(37).or.
     X          (xfac(39).eq.1.0d0.and.cfac(41).eq.1.d0)) then
           scales6 = 1.25d0
c SSB-D
        else if (xccomb(46)) then
           scales6 = 0.847455d0
c BP86
        else if(xfac(3).eq.1.0d0.and.cfac(4).eq.1.d0) then
          scales6 = 1.05d0
c B2-PLYP
        else if(xccomb(2).and.xfac(1).eq.0.53d0.and.
     &          xfac(3).eq.0.47d0.and.cfac(2).eq.0.73d0) then
          scales6 = 0.55d0
c revPBE
c         scales6 = 1.25d0
c PBE0
        else if (xccomb(10)) then
          scales6 = 0.60d0
c PW6B95
c         scales6 = 0.50d0
c B2GP-PLYP
c         scales6 = 0.40d0
         else
c default
           scales6 =1.20d0
           write(6,*)'WARNING:'
           write(6,*)'you are using HF or an untested DF together'
           write(6,*)'with the dispersion correction!'
        endif
      else if (ivdw.eq.3) then
c
c DFT-D3 version
c
        scales6  = 1.0d0
        scalesr8 = 1.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scalesr6 = 1.094d0
           scales8  = 1.682d0
c BP86 
        else if(xfac(3).eq.1.0d0.and.cfac(4).eq.1.d0) then
           scalesr6 = 1.139d0
           scales8  = 1.683d0
c B97-D
        else if (xccomb(37).or.
     X          (xfac(39).eq.1.0d0.and.cfac(41).eq.1.d0)) then
           scalesr6 = 0.892d0
           scales8  = 0.909d0
c revPBE
        else if (xfac(31).eq.1.0d0.and.cfac(12).eq.1.d0) then
          scalesr6 = 0.923d0
          scales8  = 1.010d0
c PBE/PBE96
        else if(xccomb(58).or.(xfac(7).eq.1d0.and.cfac(12).eq.1d0)) then
           scalesr6 = 1.217d0
           scales8  = 0.722d0
c PBESOL
c          scalesr6 = 1.320d0
c          scales8  = 0.581d0
c RPW86-PBE
c          scalesr6 = 1.224d0
c          scales8  = 0.901d0
c rPBE
c          scalesr6 = 0.872d0
c          scales8  = 0.514d0
c TPSS
        else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scalesr6 = 1.166d0
           scales8  = 1.105d0
c B2-PLYP
        else if(xccomb(2).and.xfac(1).eq.0.53d0.and.
     &          xfac(3).eq.0.47d0.and.cfac(2).eq.0.73d0) then
           scales6  = 0.640d0
           scalesr6 = 1.094d0
           scales8  = 1.682d0
c B3LYP
        else if (xccomb(2)) then
           scalesr6 = 1.261d0
           scales8  = 1.703d0
        else if (xccomb(10)) then
c PBE0
           scalesr6 = 1.287d0
           scales8  = 0.928d0
c revPBE38
c          scalesr6 = 1.021d0
c          scales8  = 0.862d0
c PW6B95
        else if (xccomb(25)) then
           scalesr6 = 1.532d0
           scales8  = 1.682d0
c TPSS0
c          scalesr6 = 1.252d0
c          scales8  = 1.242d0
c PWPB95
c          scales6  = 0.820d0
c          scalesr6 = 1.557d0
c          scales8  = 0.705d0
c B2GP-PLYP
c          scales6  = 0.560d0
c          scalesr6 = 1.586d0
c          scales8  = 0.760d0
c PTPSS
c          scales6  = 0.750d0
c          scalesr6 = 1.541d0
c          scales8  = 0.879d0
c HF   
c          scalesr6 = 1.158d0
c          scales8  = 1.746d0
c mpwlyp
c          scalesr6 = 1.239d0
c          scales8  = 1.098d0
c bpbe
       else if (xfac(3).eq.1d0.and.cfac(12).eq.1d0) then
          scalesr6 = 1.087d0
          scales8  = 2.033d0
c bhlyp
        else if (xccomb(55)) then
           scalesr6 = 1.370d0
           scales8  = 1.442d0
c tpssh
        else if (xccomb(18)) then
           scalesr6 = 1.223d0
           scales8  = 1.219d0
c pwb6k
        else if (xccomb(26)) then
           scalesr6 = 1.660d0
           scales8  = 0.550d0
c b1b95
        else if (xccomb(23)) then
           scalesr6 = 1.613d0
           scales8  = 1.868d0
c bop
       else if (xccomb(27)) then
          scalesr6 = 0.929d0
          scales8  = 1.975d0
c o-lyp
       else if (xfac(16).eq.1d0.and.cfac(2).eq.1d0) then
          scalesr6 = 0.806d0
          scales8  = 1.764d0
c o-pbe
       else if (xfac(16).eq.1d0.and.cfac(12).eq.1d0) then
          scalesr6 = 0.837d0
          scales8  = 1.868d0
c ssb
        else if (xccomb(46)) then
           scalesr6 = 1.215d0
           scales8  = 0.663d0
c revssb
c          scalesr6 = 1.221d0
c          scales8  = 0.560d0
c otpss
c          scalesr6 = 1.128d0
c          scales8  = 1.494d0
c b3pw91
        else if (xccomb(57)) then
           scalesr6 = 1.176d0
           scales8  = 1.775d0
c s12g
        else if (xccomb(60)) then
           scalesr6 = 1.17755954d0
           scales8  = 0.84432515d0
c s12h
        else if (xccomb(61)) then
           scalesr6 = 1.07735222d0
           scales8  = 0.37705816d0
c cam-s12g
        else if (xccomb(62)) then
           scalesr6 = 1.20250451d0
           scales8  = 0.86124355d0
c cam-s12h
        else if (xccomb(63)) then
           scalesr6 = 1.08034183d0
           scales8  = 0.37999939d0
c revpbe0
c          scalesr6 = 0.949d0
c          scales8  = 0.792d0
c pbe38
c          scalesr6 = 1.333d0
c          scales8  = 0.998d0
c mpw1b95
        else if (xccomb(21)) then
           scalesr6 = 1.605d0
           scales8  = 1.118d0
c mpwb1k
        else if (xccomb(22)) then
           scalesr6 = 1.671d0
           scales8  = 1.061d0
c bmk
c          scalesr6 = 1.931d0
c          scales8  = 2.168d0
c CAM-B3LYP
         else if (abs(xfac(40)-1.0d0).lt.1d-5.and.cfac(2).ne.0d0.
     A          and.xfac(1).ne.0d0) then
          scalesr6 = 1.378d0
          scales8  = 1.217d0
c LC-wPBE
         else if (xfac(43).eq.1.0d0.and.cfac(12).eq.1d0.
     A          and.xfac(1).eq.1d0) then
          scalesr6 = 1.355d0
          scales8  = 1.279d0
c m05
        else if (xccomb(29)) then
           scalesr6 = 1.373d0
           scales8  = 0.595d0
c m052x
        else if (xccomb(30)) then
           scalesr6 = 1.417d0
           scales8  = 0.000d0
c m06l
        else if (xccomb(33)) then
           scalesr6 = 1.581d0
           scales8  = 0.000d0
c m06
        else if (xccomb(35)) then
           scalesr6 = 1.325d0
           scales8  = 0.000d0
c m062x
        else if (xccomb(36)) then
           scalesr6 = 1.619d0
           scales8  = 0.000d0
c m06hf
        else if (xccomb(34)) then
           scalesr6 = 1.446d0
           scales8  = 0.000d0
c wb97-d3
        else if (xccomb(78)) then
           scalesr6 = 1.281d0
           scalesr8 = 1.094d0
           scales6  = 1d0
           scales8  = 1d0
c hf
        else if (xfac(1).eq.1d0) then
           scalesr6 = 1.158d0
           scales8  = 1.746d0
c dftb
c          scalesr6 = 1.699d0
c          scales8  = 1.504d0
        else
c
c     backdoor to  input parameters, e.g. for B2PLYP
c    set dft:vdw 0.64
c    set dft:vdw_sr6 1.427
c    set dft:vdw_s8 1.022
           read_ok=abs(scales6).gt.0d0.and.
     A      abs(scales8).gt.0d0.and.
     A      abs(scalesr6).gt.0d0
           if(.not.read_ok)
     c      call errquit(' vdw 3: scaling factors not available',0,0)
        endif
      else if (ivdw.eq.4) then
c
c DFT-D3BJ version
c
        scales6  = 1.0d0
cold        scalesr6 = 0.0d0
cold        scalesr8 = 0.0d0
cold        a1=0d0
cold        a2=0d0
c B2-PLYP
        if(xccomb(2).and.xfac(1).eq.0.53d0.and.
     &       xfac(3).eq.0.47d0.and.cfac(2).eq.0.73d0) then
          scales6 = 0.640d0
          a1      = 0.3065d0
          scales8 = 0.9147d0
          a2      = 5.0570d0
c B3LYP
        else if (xccomb(2)) then
           scales8  = 1.9889d0
           a1       = 0.3981d0
           a2       = 4.4211d0
c B97-D
       else if (xccomb(37)) then
          a1      = 0.5545d0
          scales8 = 2.2609d0
          a2      = 3.2297d0
c BLYP
       elseif(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scales8  = 2.6996d0
           a1       = 0.4298d0
           a2       = 4.2359d0
c BP86 
       else if(xfac(3).eq.1.0d0.and.cfac(4).eq.1.d0) then
          a1      = 0.3946d0
          scales8 = 3.2822d0
          a2      = 4.8516d0
       else if (xccomb(10)) then
c PBE0
          a1      = 0.4145d0
          scales8 = 1.2177d0
          a2      = 4.8593d0
c PBE/PBE96
       else if(xccomb(58).or.
     O         (cfac(12).eq.1d0.and.xfac(7).eq.1d0)) then
          a1      = 0.4289d0
          scales8 = 0.7875d0
          a2      = 4.4407d0
c PW6B95
       else if (xccomb(25)) then
          a1      = 0.2076d0
          scales8 = 0.7257d0
          a2      = 6.3750d0
c TPSS
        else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scales8  = 1.9435d0
           a1       = 0.4535d0
           a2       = 4.4752d0
c revPBE
        else if(xfac(31).eq.1.0d0.and.cfac(12).eq.1.d0) then
          a1      = 0.5238d0
          scales8 = 2.3550d0
          a2      = 3.5016d0
c bhlyp
        else if (xccomb(55)) then
           scales8  = 1.0354d0
           a1       = 0.2793d0
           a2       = 4.9615d0
c tpssh
        else if (xccomb(18)) then
           scales8  = 2.2382d0
           a1       = 0.4529d0
           a2       = 4.6550d0
c b3pw91
       else if (xccomb(57)) then
           scales8  = 0.4312d0
           a1       = 2.8524d0
           a2       = 4.4693d0
c pwb6k
       else if (xccomb(26)) then
           a1       = 0.1805d0
           scales8  = 0.9383d0
           a2       = 7.7627d0
c b1b95
       else if (xccomb(23)) then
           a1       = 0.2092d0
           scales8  = 1.4507d0
           a2       = 5.5545d0
c CAM-B3LYP
         else if (abs(xfac(40)-1.0d0).lt.1d-5.and.cfac(2).ne.0d0.
     A          and.xfac(1).ne.0d0) then
          a1      = 0.3708d0
          scales8 = 2.0674d0
          a2      = 5.4743d0
c LC-wPBE
       else if (abs(xfac(43)-1d0).lt.1d-5.and.
     c          abs(cfac(12)-1d0).lt.1d-5.and.
     A          abs(xfac(1) -1d0).lt.1d-5) then
          a1      = 0.3919d0
          scales8 = 1.8541d0
          a2      = 5.0897d0
c hcth120
       else if (xccomb(7)) then
           a1       = 0.3563d0
           scales8  = 1.0821d0
           a2       = 4.3359d0
c mpw1b95
       else if (xccomb(21)) then
           a1       = 0.1955d0
           scales8  = 1.0508d0
           a2       = 6.4177d0
c bop
       else if (xccomb(27)) then
           a1       = 0.4870d0
           scales8  = 3.2950d0
           a2       = 3.5043d0
c o-lyp
       else if (xfac(16).eq.1d0.and.cfac(2).eq.1d0) then
           a1       = 0.5299d0
           scales8  = 2.6205d0
           a2       = 2.8065d0
c bpbe
       else if (xfac(3).eq.1d0.and.cfac(12).eq.1d0) then
           a1       = 0.4567d0
           scales8  = 4.0728d0
           a2       = 4.3908d0
c o-pbe
       else if (xfac(16).eq.1d0.and.cfac(12).eq.1d0) then
           a1       = 0.5512d0
           scales8  = 3.3816d0
           a2       = 2.9444d0
c ssb
        else if (xccomb(46)) then
           a1       = -0.0952d0
           scales8  = -0.1744d0
           a2       =  5.2170d0
c mpwb1k
        else if (xccomb(22)) then
           scales8  =  0.9499d0
           a1       =  0.1474d0
           a2       =  6.6223d0
c hf
        else if (xfac(1).eq.1d0) then
           scales8  =  0.9171d0
           a1       =  0.3385d0
           a2       =  2.8830d0
c scan
        else if (xccomb(66)) then
           a1       =  0.5380d0
           scales8  =  0.0000d0
           a2       =  5.4200d0
c r-scan
chttps://github.com/awvwgk/r2scan-d4-paper
        else if (xccomb(71)) then
           a1      = 0.47023427d0
           scales8 = 1.08859014d0
           a2      = 5.73408312d0
c r2scan
chttps://github.com/awvwgk/r2scan-d4-paper
        else if (xccomb(73)) then
           a1      = 0.49484001d0
           scales8 = 0.78981345d0
           a2      = 5.73083694
        else
c     try to see if set rtdb was used
           read_ok=abs(scales8).gt.0d0.and.
     A      abs(a1).gt.0d0.and.
     A      abs(a2).gt.0d0
           if(.not.read_ok)
     c    call errquit(' vdw 4: scaling factors not available',0,0)
        endif
      else
        scales6 = 0.0d0
        scales8 = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
        a1 = 0.0d0
        a2 = 0.0d0
      endif
      return
      end
c
c     Small function to extract coordination number representation number
c     and atomic number from column two and three of the table
c
      subroutine xc_vdwlimit(iat,jat,iadr,jadr)
      implicit none
      integer iat,jat,iadr,jadr
c
      iadr=1
      jadr=1
810   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 810
      endif
820   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 820
      endif
      return
      end
c
      double precision function xc_neff(z)
      implicit none
      integer z
c
      integer nv
c
      if(z.eq.1) then
         xc_neff=0.8d0
         return
      endif
      nv=0
      if(z.le.2) nv=z
      if(z.le.10) nv=z-2
      if(z.gt.10) stop 'neff stop '
      xc_neff=1.17d0+0.33d0*nv
      return
      end
