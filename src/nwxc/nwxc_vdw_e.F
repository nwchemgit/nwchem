c
c     empirical dispersion: energy
c
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the dispersion energy
C>
C> This function evaluates the dispersion energy based on an empirical
C> expression. This routine supports multiple expressions commonly 
C> used in density functional theory. 
C>
C> The DFT-D3 correction with BJ damping [4,5] is given by
C> \f{eqnarray*}{
C>   E_{\mathrm{disp}}^{D3(BJ)} &=& -\frac{1}{2}\sum_{A\ne B}
C>     s_6\frac{C_6^{AB}}{R_{AB}^6+[f(R^0_{AB})]^6} +
C>     s_8\frac{C_8^{AB}}{R_{AB}^8+[f(R^0_{AB})]^8} \\\\
C>   f(R^0_{AB}) &=& a_1R^0_{AB}+a_2 \\\\
C>   R^0_{AB} &=& \sqrt{\frac{C_8^{AB}}{C_6^{AB}}}
C> \f}
C>
C> \return The dispersion energy \f$E_{\mathrm{disp}}\f$.
C>
C> ### References ###
C>
C> [1] S. Grimme,
C>     "Accurate description of van der Waals complexes by density
C>      functional theory including empirical corrections",
C>     J. Comp. Chem. (2004) <b>25</b>, 1463-1473, DOI:
C>     <a href="https://doi.org/10.1002/jcc.20078">
C>     10.1002/jcc.20078</a>.
C>
C> [2] U. Zimmerli, M. Parrinello, P. Koumoutsakos,
C>     "Dispersion corrections to density functionals for water
C>      aromatic interactions",
C>     J. Chem. Phys. (2004) <b>120</b>, 2693, DOI:
C>     <a href="https://doi.org/10.1063/1.1637034">
C>     10.1063/1.1637034</a>.
C>
C> [3] Q. Wu, W. Yang,
C>     "Empirical correction to density functional theory for van der
C>      Waals interactions",
C>     J. Chem. Phys. (2002) <b>116</b>, 515, DOI:
C>     <a href="https://doi.org/10.1063/1.1424928">
C>     10.1063/1.1424928</a>.
C>
C> [4] A.D. Becke, E.R. Johnson,
C>     "A unified density-functional treatment of dynamical,
C>      nondynamical and dispersion correlations",
C>     J. Chem. Phys. (2007) <b>127</b> 124108, DOI:
C>     <a href="https://doi.org/10.1063/1.2768530">
C>     10.1063/1.2768530</a> (See appendix C).
C>
C> [5] S. Grimme, S. Ehrlich, L. Goerigk,
C>     "Effect of the damping function in dispersion corrected
C>      density functional theory", J. Comput. Chem. (2011)
C>     <b>32</b>, pp. 1456-1465, DOI:
C>     <a href="https://doi.org/10.1002/jcc.21759">
C>     10.1002/jcc.21759</a> (See Eqs.(5-6)).
C>
      double precision function nwxc_vdw_e(s6,s8,sr6,sr8,a1,a2,n,x,z)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
c
#include "nwxc_vdw.fh"
c
      double precision s6     !< [Input] The \f$s_6\f$ coefficient
      double precision s8     !< [Input] The \f$s_8\f$ coefficient
      double precision sr6    !< [Input] The \f$s_{r,6}\f$ coefficient
      double precision sr8    !< [Input] The \f$s_{r,8}\f$ coefficient
      double precision a1     !< [Input] The \f$a_1\f$ coefficient
      double precision a2     !< [Input] The \f$a_2\f$ coefficient
      integer n               !< [Input] The number of atoms
      double precision x(3,n) !< [Input] The atomic coordinates
      integer z(n)            !< [Input] The atomic numbers of the atoms
c
      integer i,j
      double precision nwxc_fdmp, nwxc_fdmp3, nwxc_fdmpbj, cni, cnj,c6d3
      double precision nwxc_c6ij_sk,rij,nwxc_c6cn,nwxc_crd_nr,e6,e8
      double precision r0bj
      external nwxc_c6ij_sk,nwxc_c6cn,nxtask,nwxc_crd_nr,
     &         nwxc_fdmp,nwxc_fdmp3,nwxc_fdmpbj
c
      nwxc_vdw_e=0d0
      e6=0.0d0
      e8=0.0d0
c
c DFT-D1 / DFT-D2
c
      if (ivdw.le.2) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
               nwxc_vdw_e=nwxc_vdw_e-nwxc_c6ij_sk(i,j,z)*
     *            nwxc_fdmp(rij,r0(z(i))+r0(z(j)))*
     *            (rij)**(-6.0d0)
            enddo
          endif
        enddo
        nwxc_vdw_e=nwxc_vdw_e*s6
c
c DFT-D3
c
c As off August, 2011 Grimme states: "Adding three-body corrections is
c currently not recommended, as very little is known about the three-
c body behaviour of common DFs in overlapping density regions."
c http://toc.uni-muenster.de/DFTD3/data/man.pdf, section 1.3.
c Hence the three-body terms have not been implemented.
c
c The reference to three-center derivatives in the gradient code
c refers to contributions that come from differentiating the
c coordination dependent dispersion coefficients.
c
      else if (ivdw.eq.3) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
               cni=nwxc_crd_nr(i,n,x,z)
               cnj=nwxc_crd_nr(j,n,x,z)
               c6d3=nwxc_c6cn(z(i),z(j),cni,cnj)
               c8=3.0d0*c6d3*Qatom(z(i))*Qatom(z(j))
               e6=e6-c6d3*nwxc_fdmp3(rij,r0AB(z(i),z(j))*sr6,alpha)*
     *            (rij)**(-6.0d0)
               e8=e8-c8*nwxc_fdmp3(rij,r0AB(z(i),z(j))*sr8,alpha+2.0d0)*
     *            (rij)**(-8.0d0)
            enddo
          endif
        enddo
        nwxc_vdw_e=e6*s6+e8*s8
      else if (ivdw.eq.4) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
               cni=nwxc_crd_nr(i,n,x,z)
               cnj=nwxc_crd_nr(j,n,x,z)
               c6d3=nwxc_c6cn(z(i),z(j),cni,cnj)
               c8=3.0d0*c6d3*Qatom(z(i))*Qatom(z(j))
               r0bj=dsqrt(c8/c6d3)
               e6=e6-c6d3*nwxc_fdmpbj(rij,r0bj,a1,a2,6)
               e8=e8-c8*nwxc_fdmpbj(rij,r0bj,a1,a2,8)
            enddo
          endif
        enddo
        nwxc_vdw_e=e6*s6+e8*s8
      endif
c
      return
      end
C>
C> @}
c $Id$
