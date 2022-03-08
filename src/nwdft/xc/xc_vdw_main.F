c
c     $Id$
c
      subroutine xc_vdw_der(s6,s8,sr6,sr8,a1,a2,n,x,z,force)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "xc_vdw.fh"
#include "global.fh"
      double precision s6,s8,sr6,sr8,a1,a2
      integer n
      double precision x(3,n),force(3,n)
      integer z(n)
c
      integer i,j,k,A,l_cnij,k_cnij,l_cnijk,k_cnijk
      double precision c6ij_sk
      external c6ij_sk
      double precision drajdxa
      double precision ff1,rr,ff
      double precision fdmp,f1dmp,cnA,cnj
      external c6cn,crd_nr
      double precision c6cn,crd_nr
      double precision fac6,fac8,fdmp6,fdmp6a,fdmp8,fdmp8a,Qfac
      double precision rAj,rAk,rjk,r0aj,r0ak,r0jk,c6Aj,grad_c6(3)
      double precision dxAj,dyAj,dzAj,dxAk,dyAk,dzAk,dxjk,dyjk,dzjk
      double precision tmp6,tmp6a,tmp8,tmp8a
      double precision xc_fdmpbj, xc_fdmpbj_d1
      external         xc_fdmpbj, xc_fdmpbj_d1
      double precision dEdrij
      double precision dc6(n),cns(n),dc6i,dc6j
      double precision dEdc6,rAj2,dcn,rij(3)
      double precision dcndrij
      external dcndrij

c
c     Derivatives of Grimme dispersion term
c
c  DFT-D1 / DFT-D2
c
      if (ivdw.le.2) then
         do A=1,n
            force(1,A)=0d0
            force(2,A)=0d0
            force(3,A)=0d0
            if (Z(A).ne.0) then
              do j=1,n
                 if(A.ne.j) then
                    rAj=sqrt(
     +                 (x(1,A)-x(1,j))**2 +
     +                 (x(2,A)-x(2,j))**2 +
     +                 (x(3,A)-x(3,j))**2)
c     protect from NaNs caused by bqs
               if(raj.lt.1d30.or.abs(raj).lt.1d-8) then
                    r0aj=r0(z(A))+r0(z(j))
                    ff= fdmp(rAj,r0aj)
                    ff1= f1dmp(rAj,r0aj,ff)
                    rr=c6ij_sk(A,j,z)/(rAj**6)*
     *               ((-6d0*ff/rAj)+ff1)
                    do i=1,3
                       drAjdxa=(x(i,A)-x(i,j))/rAj
                       force(i,A)=force(i,A)-rr*drAjdxa
                    enddo
                    endif
                 endif
              enddo
            endif
         enddo
         if(abs(s6-1d0).gt.1d-9) 
     F        call dscal(3*n,s6,force,1)
c
c DFT-D3
c
      else if (ivdw.eq.3) then
c
c Precompute coordination numbers
c
         do A=1,n
           if (Z(A).eq.0) cycle
           cns(A) = crd_nr(A,n,x,z)
         enddo
c
c Get derivatives with respect to CNs
c
         dc6(:) = 0d0
         do A=2,n
           if (Z(A).eq.0) cycle
           do j=1,A-1
             if (Z(j).eq.0) cycle
             rAj = sum((x(:,j) - x(:,A))**2)
             if (rAj.gt.9000.0d0) cycle
             rAj = dsqrt(rAj)
c
c            Factors
c
             r0aj=r0AB(z(A),z(j))
             Qfac=Qatom(z(A))*Qatom(z(j))
             fac6=(rAj/(sr6*r0aj))**(-alpha)
             fac8=(rAj/(sr8*r0aj))**(-(alpha+2.0d0))
             fdmp6=1.0d0/(1.0d0+6.0d0*fac6)
             fdmp8=1.0d0/(1.0d0+6.0d0*fac8)

             ! Derivative of the C6 coeff with respect to CNs
             call c6_grad(A,j,x,z,n,dc6i,dc6j,cns(A),cns(j))

             ! Derivative of Edisp with respect to C6
             dEdc6 = fdmp6*s6/rAj**6 + 3d0*s8*Qfac*fdmp8/rAj**8

             ! Derivative of Edisp with respect to CNs
             dc6(A) = dc6(A) + dEdc6*dc6i
             dc6(j) = dc6(j) + dEdc6*dc6j
           enddo
         enddo
c
c Get the gradients
c
         force(:,:) = 0d0
         do A=2,n
           if (Z(A).eq.0) cycle
           do j=1,A-1
             if (Z(j).eq.0) cycle
             rij(:) = x(:,j) - x(:,A)

             rAj2 = sum(rij**2)
             rAj = dsqrt(rAj2)

             ! Derivative of the CNs with respect to rij
             if (rAj2.gt.1600.d0) then
               dcn = 0d0
             else
               dcn = dcndrij(A,j,rAj,z)
             endif

             ! Intermediate factors
             cnA = cns(A)
             cnj = cns(j)
             r0aj=r0AB(z(A),z(j))
             Qfac=Qatom(z(A))*Qatom(z(j))
             c6Aj=c6cn(z(A),z(j),cnA,cnj)
             c8=3.0d0*c6Aj*Qfac
             fac6=(rAj/(sr6*r0aj))**(-alpha)
             fac8=(rAj/(sr8*r0aj))**(-(alpha+2.0d0))
             fdmp6=1.0d0/(1.0d0+6.0d0*fac6)
             fdmp8=1.0d0/(1.0d0+6.0d0*fac8)
             tmp6=6.0d0*alpha*fac6*fdmp6**2/(rAj)
             tmp8=6.0d0*(alpha+2d0)*fac8*fdmp8**2/(rAj)

             ! Derivative of Edisp with respect to rij
             dEdrij = -6d0*s6*c6Aj*fdmp6/rAj**7 + s6*c6Aj*tmp6/rAj**6 -
     $                 24d0*s8*c6Aj*Qfac*fdmp8/rAj**9 + 
     $                 3d0*s8*c6Aj*Qfac*tmp8/rAj**8 +
     $                 dcn*(dc6(A) + dc6(j))

             ! Nuclear gradients
             force(:,A) = force(:,A) + dEdrij*rij/rAj
             force(:,j) = force(:,j) - dEdrij*rij/rAj
           enddo
         enddo
c
c DFT-D3BJ
c
      else if (ivdw.eq.4) then
c
c Precompute coordination numbers
c
         do A=1,n
           if (Z(A).eq.0) cycle
           cns(A) = crd_nr(A,n,x,z)
         enddo
c
c Get derivatives with respect to CNs
c
         dc6(:) = 0d0
         do A=2,n
           if (Z(A).eq.0) cycle
           do j=1,A-1
             if (Z(j).eq.0) cycle
             rAj = sum((x(:,j) - x(:,A))**2)
             if (rAj.gt.9000.0d0) cycle
             rAj = dsqrt(rAj)
             Qfac=Qatom(z(A))*Qatom(z(j))
             r0Aj=dsqrt(3.0d0*Qfac)

             ! Derivative of the C6 coeff with respect to CNs
             call c6_grad(A,j,x,z,n,dc6i,dc6j,cns(A),cns(j))

             ! Derivative of Edisp with respect to C6
             dEdc6 = s6*xc_fdmpbj(rAj,r0Aj,a1,a2,6) +
     $               s8*3d0*Qfac*xc_fdmpbj(rAj,r0Aj,a1,a2,8)

             ! Derivative of Edisp with respect to CNs
             dc6(A) = dc6(A) + dEdc6*dc6i
             dc6(j) = dc6(j) + dEdc6*dc6j
           enddo
         enddo
c
c Get the gradients
c
         force(:,:) = 0d0
         do A=2,n
           if (Z(A).eq.0) cycle
           do j=1,A-1
             if (Z(j).eq.0) cycle
             rij(:) = x(:,j) - x(:,A)

             rAj2 = sum(rij**2)
             rAj = dsqrt(rAj2)

             ! Derivative of the CNs with respect to rij
             if (rAj2.gt.1600.d0) then
               dcn = 0d0
             else
               dcn = dcndrij(A,j,rAj,z)
             endif

             ! Intermediate factors
             cnA = cns(A)
             cnj = cns(j)
             Qfac=Qatom(z(A))*Qatom(z(j))
             c6Aj=c6cn(z(A),z(j),cnA,cnj)
             c8=3.0d0*c6Aj*Qfac
             r0Aj=dsqrt(3.0d0*Qfac)

             ! Derivative of Edisp with respect to rij
             dEdrij =
     $       s6*c6Aj*xc_fdmpbj_d1(rAj,r0Aj,a1,a2,6) +
     $       s8*c8*xc_fdmpbj_d1(rAj,r0Aj,a1,a2,8) +
     $       dcn*(dc6(A) + dc6(j))

             ! Nuclear gradients
             force(:,A) = force(:,A) + dEdrij*rij/rAj
             force(:,j) = force(:,j) - dEdrij*rij/rAj
           enddo
         enddo
c
      endif
c
#ifdef DEBUG
      write(6,*) ' gradient vdw called'
#endif
      return
      end
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
      double precision function xc_vdw_e(s6,s8,sr6,sr8,a1,a2,n,x,z)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
#include "errquit.fh"
#include "xc_vdw.fh"
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
      integer i6, i8
      parameter(i6 = 6)
      parameter(i8 = 8)
      double precision fdmp, fdmp3, cni, cnj,c6d3,xc_fdmpbj
      double precision c6ij_sk,rij,c6cn,crd_nr,e6,e8,r0bj
      external c6ij_sk,c6cn,nxtask,crd_nr,fdmp3,xc_fdmpbj
c
      xc_vdw_e=0d0
      e6=0.0d0
      e8=0.0d0
c
c DFT-D1 / DFT-D2
c
      if (ivdw.le.2) then
        do i=1,n-1
          if (Z(i).ne.0) then
            if (r0(z(i)).le.0.0d0) then
              call errquit("xc_vdw_e: no Grimme parameters for element",
     +                     int(z(i)),UERR)
            endif
            do j=i+1,n
               if (Z(j).ne.0) then
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
c     protect from NaNs caused by bqs
               if(rij.gt.1d30.or.abs(rij).lt.1d-4) then
                   xc_vdw_e=0d0
               else
               xc_vdw_e=xc_vdw_e-c6ij_sk(i,j,z)*
     *            fdmp(rij,r0(z(i))+r0(z(j)))*
     *            (rij)**(-6.0d0)
               endif
            endif
            enddo
          endif
        enddo
        xc_vdw_e=xc_vdw_e*s6
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
               if (Z(j).ne.0) then
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
c     protect from NaNs caused by bqs
               if(rij.lt.1d30.or.abs(rij).lt.1d-8) then
               cni=crd_nr(i,n,x,z)
               cnj=crd_nr(j,n,x,z)
               c6d3=c6cn(z(i),z(j),cni,cnj)
               c8=3.0d0*c6d3*Qatom(z(i))*Qatom(z(j))
               e6=e6-c6d3*fdmp3(rij,r0AB(z(i),z(j))*sr6,alpha)*
     *            (rij)**(-6.0d0)
               e8=e8-c8*fdmp3(rij,r0AB(z(i),z(j))*sr8,alpha+2.0d0)*
     *            (rij)**(-8.0d0)
               endif
               endif
            enddo
          endif
        enddo
        xc_vdw_e=e6*s6+e8*s8
      else if (ivdw.eq.4) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               if (Z(j).ne.0) then
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
c     protect from NaNs caused by bqs
               if(rij.lt.1d30.or.abs(rij).lt.1d-8) then
               cni=crd_nr(i,n,x,z)
               cnj=crd_nr(j,n,x,z)
               c6d3=c6cn(z(i),z(j),cni,cnj)
               c8=3.0d0*c6d3*Qatom(z(i))*Qatom(z(j))
               r0bj=dsqrt(c8/c6d3)
               e6=e6-c6d3*xc_fdmpbj(rij,r0bj,a1,a2,i6)
               e8=e8-c8*xc_fdmpbj(rij,r0bj,a1,a2,i8)
               endif
               endif
            enddo
          endif
        enddo
        xc_vdw_e=e6*s6+e8*s8
      endif
      return
      end
c
      subroutine xc_vdw(rtdb,geom,exc,force,what)
      implicit none
      character *(*) what
      integer geom,rtdb
      double precision exc,force(*),s6,s8,sr6,sr8,a1,a2
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "xc_vdw.fh"
c
      integer n
      integer itags,ltags,i_xyz,l_xyz,icharge,lcharge,
     I     l_fvdw,i_fvdw, i_xyz2,l_xyz2,i_iz2,l_iz2
      external xc_vdw_e
      double precision xc_vdw_e,evdw,scalea
      integer iz,lz,i
      logical xc_vdw_init
      external xc_vdw_init
      logical oprint,oprinth
      logical stat
      logical use_nwxc_disp,out1
c
      double precision delta,delta_default
      double precision s6_in
c
      oprint = util_print('vdw', print_medium)
      oprinth = util_print('vdw high', print_high)
c
c     Allocate memory blocks
c
      if (.not. geom_ncent(geom, n)) 
     &   call errquit('xcvdw: geom_ncent failed',geom, GEOM_ERR)
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_xyz,i_xyz))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,n,'atns',lz,iz))
     &   call errquit('xcvdw: cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,n*16,'tags',ltags,itags))
     &   call errquit('xcvdw: cannot allocate tags',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,n,'charge',lcharge,icharge))
     &   call errquit('xcvdw: cannot allocate charge',0, MA_ERR)
      if (.not. geom_cart_get2(geom, n, Byte_MB(itags),
     &        Dbl_MB(i_xyz), Dbl_MB(icharge), int_mb(iz)))
     &   call errquit('xcvdw: geom_cart_get failed',74, GEOM_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('xcvdw: cannot pop stack',3, MA_ERR)
c
c     Which Grimme dispersion version
c
      if (.not.rtdb_get(rtdb,'dft:ivdw',MT_INT,1,ivdw))
     &      ivdw = 2
c
c     get rid of bqs
c      
      if (.not.MA_push_get(mt_dbl,n*3,'xyz',l_xyz2,i_xyz2))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_push_get(mt_int,n,'iz',l_iz2,i_iz2))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      call xc_vdw_nobqs(geom, Byte_MB(itags),
     N     n, dbl_mb(i_xyz), int_mb(iz),
     T     dbl_mb(i_xyz2),int_mb(i_iz2))
      if (.not.ma_chop_stack(l_xyz2)) 
     C   call errquit('xcvdw: cannot pop stack',14, MA_ERR)
c
      use_nwxc_disp = .false.
      if(util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_has_disp",out1)
         use_nwxc_disp = out1
         if (use_nwxc_disp) then
c       get ivdw here as this data is needed for xc_vdw_init
            call nwxc_get_disp(ivdw,s6,s8,sr6,sr8,a1,a2,alpha)
         endif
      endif
c
c     conversion factor angs 2 au
c
       if(.not.geom_get_ang2au(geom, scalea)) call
     S     errquit('xcvdw: gang2au failed',0,0)
c
c     Initialize are variables
c
      if(.not.xc_vdw_init(scalea))
     &   call errquit('xcvdw: vwdinit failed',0, 0)
c
c     Read in some user defined parameters
c
      if (.not.rtdb_get(rtdb,'dft:vdwalpha',MT_DBL,1,alpha))
     &      alpha = 20.0d0
      if (ivdw.eq.3) alpha = 14.0d0
      if (ivdw.eq.4) alpha = 14.0d0
c
c     Get proper scaling factors depending on Grimme dispersion version
c
      if (.not.use_nwxc_disp) then
        if(.not.rtdb_get(rtdb,'dft:vdw_s6',mt_dbl,1,s6)) s6=0d0
        if(.not.rtdb_get(rtdb,'dft:vdw_s8',mt_dbl,1,s8)) s8=0d0
        if(.not.rtdb_get(rtdb,'dft:vdw_sr6',mt_dbl,1,sr6)) sr6=0d0
        if(.not.rtdb_get(rtdb,'dft:vdw_sr8',mt_dbl,1,sr8)) sr8=0d0
        if(.not.rtdb_get(rtdb,'dft:vdw_a1',mt_dbl,1,a1)) a1=0d0
        if(.not.rtdb_get(rtdb,'dft:vdw_a2',mt_dbl,1,a2)) a2=0d0
        call get_scaling_fac(s6,s8,sr6,sr8,a1,a2)
        if(rtdb_get(rtdb, 'dft:vdw', mt_dbl, 1, s6_in)) then
           s6=s6_in
           if(ga_nodeid().eq.0)write(6,*) ' WARNING: vdw s6 set = ',s6
        endif
      else
c       get the scaling factors here (just being paranoid)
        call nwxc_get_disp(ivdw,s6,s8,sr6,sr8,a1,a2,alpha)
      endif
c
      if(what.eq.'energy') then
c
c     Compute energy contribution
c
        if(oprinth.and.ga_nodeid().eq.0) then
           write(luout,*) ' s6 =',s6
           write(luout,*) ' s8 =',s8
           write(luout,*) ' sr6 =',sr6
           write(luout,*) ' sr8 =',sr8
           write(luout,*) ' alpha  =',alpha
           write(luout,*) ' ivdw  =',ivdw
         endif
c
         evdw=xc_vdw_e(s6,s8,sr6,sr8,a1,a2,n,dbl_mb(i_xyz),int_mb(iz))
c
         if(oprint.and.ga_nodeid().eq.0) then
            write(luout,*) 
     D           '     Dispersion Parameters'
            write(luout,*)
     D           '     ---------------------'
           if (ivdw.eq.1.or.ivdw.eq.2) then
              write(luout,222) s6, evdw
 222  format(/
     &     '          s6 scale factor :', f22.12/
     &     '              vdW contrib :', f22.12/)
           endif
           if (ivdw.eq.3) then
              write(luout,223) s6, s8, sr6, sr8, evdw
 223  format(/
     &     '             DFT-D3 Model   ', /
     &     '          s6 scale factor  :', f22.12/         
     &     '          s8 scale factor  :', f22.12/
     &     '          sr6 scale factor :', f22.12/
     &     '          sr8 scale factor :', f22.12/
     &     '              vdW contrib  :', f22.12/)
           endif
           if (ivdw.eq.4) then
              write(luout,224) s6, s8, a1, a2, evdw
 224  format(/
     &     '           DFT-D3BJ Model   ', /
     &     '          s6 scale factor  :', f22.12/
     &     '          s8 scale factor  :', f22.12/
     &     '          a1 parameter     :', f22.12/
     &     '          a2 parameter     :', f22.12/
     &     '              vdW contrib  :', f22.12/)
           endif
         endif
c
c        Add contribution to Exc
c
         Exc=Exc+evdw
c
      elseif(what.eq.'forces') then
c
c     Gradient calculation
c
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_fvdw,i_fvdw))
     &   call errquit('xcvdw: cannot allocate forcev',0, MA_ERR)
c
         call xc_vdw_der(s6,s8,sr6,sr8,a1,a2,n,dbl_mb(i_xyz),int_mb(iz),
     D     dbl_mb(i_fvdw))
c
         if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' vdW contrib for S6=',s6
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               dbl_mb(i_fvdw+3*(i-1)),
     Y               dbl_mb(i_fvdw+3*(i-1)+1),
     Z        dbl_mb(i_fvdw+3*(i-1)+2)
             enddo
             write(luout,*) ' before vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo

          endif
c
c         Add to force matrix
c
          call daxpy(3*n,1d0,dbl_mb(i_fvdw),1,force,1)
c
          if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' after vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo
          endif
c
      elseif(what.eq.'hessian') then
c
c     Hessian calculation, numerical from analytical gradients
c
c     Get delta as used in a numerical hessian DFT calculation
c
      delta_default  =  0.01d0
      if (.not.rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,delta))
     &      delta = delta_default 
c
        call xc_vdw_hess(delta,s6,s8,sr6,sr8,a1,a2,n,
     &       dbl_mb(i_xyz),int_mb(iz))
c
        if(oprint.and.ga_nodeid().eq.0)  then
            write(luout,*) ' s6 = ',s6
            write(luout,*) ' vdw to hessian contribution is done'
        endif
      endif
c
c     Clean up
c
      if (.not.ma_chop_stack(l_xyz)) 
     C   call errquit('xcvdw: cannot pop stack',4, MA_ERR)
c
      return
      end
c
      subroutine xc_vdw_hess(delta,s6,s8,sr6,sr8,a1,a2,n,x,z)
      implicit none
c This function makes vdw empirical correction to the hessian
c must be called before thermochemical data and vibrational
c analysis are generated.
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      double precision s6,s8,sr6,sr8,a1,a2
      integer n                
      double precision x(3,n)
      integer z(n)
      double precision l_force(3,n)
      double precision r_force(3,n)
      double precision hessvdw(3,n,3,n)
      double precision dval0
      double precision dval1
      double precision delta
      integer i, j, A, B
      integer n3xyz 
      integer nat2
      integer nhesst
      integer lenhess 
      integer l_exybs,k_exybs 
      integer l_exybt,k_exybt
c
      character*(nw_max_path_len) filehess
c
c  dispersion contribution to hessian 
c
c in principle this task is very fast,
c so only master node works, and read-write to disk.
c 
      if (ga_nodeid().eq.0) then
c
      call util_file_name('hess', .false., .false.,filehess)

      lenhess=inp_strlen(filehess)
      n3xyz=3*n
      nhesst=n3xyz*(n3xyz+1)/2
      nat2=n3xyz*n3xyz
c
      if (.not.ma_push_get(mt_dbl,nat2,'xcvdwhess exybs ',
     *  l_exybs,k_exybs))
     &   call errquit('xcvdwhess: cannot allocate exybs',0, MA_ERR)
      call dfill(nat2,0.0d00,dbl_mb(k_exybs),1)
c
      if (.not.ma_push_get(mt_dbl,nhesst,'xcvdwhess exybt ',
     *  l_exybt,k_exybt))
     &   call errquit('xcvdwhess: cannot allocate exybt',0, MA_ERR)
      call dfill(nhesst,0.0d00,dbl_mb(k_exybt),1)
c
      write(LuOut,* ) 'Read old hessian file : ', filehess 
c
      call dfill(nhesst,0.d0,dbl_mb(k_exybt),1)
      !write(6,* ) 'leee '
      open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')
      do j = 0,(nhesst-1)
        read(69,*,err=99901,end=99902) dval0
        dbl_mb(k_exybt+j) = dval0
      !write(6,* ) 'dval ', dval0
      enddo
      close(unit=69,status='keep')
c
c     complete the square matrix from triangle values
c
      call trin2squa(dbl_mb(k_exybs),dbl_mb(k_exybt),n3xyz)
c
      write(LuOut,* ) 'vdW contribution to hessian '
c
      call output(dbl_mb(k_exybs), 1, n3xyz, 1, n3xyz, n3xyz, n3xyz, 1)
c
       do A = 1, n
         do i = 1, 3
           do B=1, n
            do j=1, 3
               r_force(j,B)=0.0d0
               l_force(j,B)=0.0d0
            enddo
           enddo
c right displacement
            x(i,A) = x(i,A) + delta
            call xc_vdw_der(s6,s8,sr6,sr8,a1,a2,n,x,z,r_force)
c left displacement
            x(i,A) = x(i,A) - 2.0d0*delta
            call xc_vdw_der(s6,s8,sr6,sr8,a1,a2,n,x,z,l_force)
c back to original position
            x(i,A) = x(i,A) + delta
c
            do B=1,n
              do  j=1, 3
                dval1 = (r_force(j,B)-l_force(j,B)) / (2.0d00*delta)
                hessvdw(i,A,j,B)=dval1
              enddo
            enddo
         enddo
      enddo
c
      call output(hessvdw,1,n3xyz,1,n3xyz,
     &            n3xyz,n3xyz,1)
c
      call daxpy(nat2,1d0,hessvdw,1,dbl_mb(k_exybs),1)
c
c:write the final hessian
c
      call output(dbl_mb(k_exybs),1,n3xyz,1,n3xyz,
     &            n3xyz,n3xyz,1)
      call  stpr_wrt_fd_from_sq(dbl_mb(k_exybs),n3xyz,filehess)
      write(LuOut,* ) 'New hessian file vdw corrected has been 
     .                      written:', filehess 
      if (.not.ma_chop_stack(l_exybs)) 
     C   call errquit('xcvdwhess: cannot pop stack exybs',4, MA_ERR)

      endif
c
      call ga_sync
c
      return
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99900', 911, DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99901', 911, DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99902', 911, DISK_ERR)
      end
c
      subroutine xc_vdw_to_hessian(rtdb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb
      integer  geom
      double precision dum
      double precision dum1(1)
      character*255 name
c
      if (.not. rtdb_cget(rtdb,'geometry', 1, name))
     $     name = 'geometry'
c
      if (.not. geom_create(geom, name))
     $ call errquit('xc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb, geom, name))
     $ call errquit('xc_vdw_to_hessian: no geometry load form rtdb', 0,
     $        GEOM_ERR)
c
      call xc_vdw(rtdb, geom, dum, dum1, 'hessian')
c
      if(.not. geom_destroy(geom))
     $ call errquit('xc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      return
      end
      subroutine xc_vdw_nobqs(geom, tags,
     N     n, xyz, z,
     T     tempxyz, tempz)
      implicit none
#include "geom.fh"
#include "inp.fh"
      double precision xyz(3,*),tempxyz(3,*)
      character*16 tags(*)
      integer n
      integer z(*),tempz(*)
      integer geom
c
      integer i,n_left
c
      do i=1,n
         tempz(i)=z(i)
         tempxyz(1,i)=xyz(1,i)
         tempxyz(2,i)=xyz(2,i)
         tempxyz(3,i)=xyz(3,i)
      enddo
      n_left=0
      do i=1,n
               n_left=n_left+1
               z(n_left)=tempz(i)
            if (inp_compare(.false.,tags(i)(1:2),'bq')) then 
c move bq at a crazy large distance to kill its contribution
               xyz(1,n_left)=1d30
               xyz(2,n_left)=1d30
               xyz(3,n_left)=1d30
            else
               xyz(1,n_left)=tempxyz(1,i)
               xyz(2,n_left)=tempxyz(2,i)
               xyz(3,n_left)=tempxyz(3,i)
            endif
      enddo
      n=n_left
      return
      end
