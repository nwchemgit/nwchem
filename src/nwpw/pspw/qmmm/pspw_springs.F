c
c $Id$
c


*     *************************************
*     *                                   *
*     *    pspw_qmmm_spring_cbond_frag    *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_spring_cbond_frag(ncbond,indx,Kr0,rion)
      implicit none
      integer ncbond,indx(4,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  r1,x1,y1,z1
      real*8  r2,x2,y2,z2
      real*8  E

      E = 0.0d0
      do b=1,ncbond
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)
 
        x1 = rion(1,ii)-rion(1,jj)
        y1 = rion(2,ii)-rion(2,jj)
        z1 = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x1,y1,z1)
        r1 = dsqrt( x1**2 + y1**2 + z1**2)

        x2 = rion(1,kk)-rion(1,ll)
        y2 = rion(2,kk)-rion(2,ll)
        z2 = rion(3,kk)-rion(3,ll)
        call lattice_min_difference(x2,y2,z2)
        r2 = dsqrt( x2**2 + y2**2 + z2**2)

        E = E + Kr0(3,b)*(r1-Kr0(1,b))*(r2-Kr0(2,b))
      end do

      pspw_qmmm_spring_cbond_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *    pspw_qmmm_spring_cbond_frag_fion    *
*     *                                        *
*     ******************************************
      subroutine pspw_qmmm_spring_cbond_frag_fion(ncbond,indx,Kr0,rion,
     >                                           fion)
      implicit none
      integer ncbond,indx(4,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  x1,y1,z1,r1
      real*8  x2,y2,z2,r2
      real*8  dE1,dE2

      do b=1,ncbond
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)

        x1 = rion(1,ii)-rion(1,jj)
        y1 = rion(2,ii)-rion(2,jj)
        z1 = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x1,y1,z1)
        r1 = dsqrt( x1**2 + y1**2 + z1**2)

        x2 = rion(1,kk)-rion(1,ll)
        y2 = rion(2,kk)-rion(2,ll)
        z2 = rion(3,kk)-rion(3,ll)
        call lattice_min_difference(x2,y2,z2)
        r2 = dsqrt( x2**2 + y2**2 + z2**2)

        dE1 =  Kr0(3,b)*(r2-Kr0(2,b))/r1
        dE2 =  Kr0(3,b)*(r1-Kr0(1,b))/r2

        fion(1,ii) = fion(1,ii) - x1*dE1
        fion(2,ii) = fion(2,ii) - y1*dE1
        fion(3,ii) = fion(3,ii) - z1*dE1

        fion(1,jj) = fion(1,jj) + x1*dE1
        fion(2,jj) = fion(2,jj) + y1*dE1
        fion(3,jj) = fion(3,jj) + z1*dE1

        fion(1,kk) = fion(1,kk) - x2*dE2
        fion(2,kk) = fion(2,kk) - y2*dE2
        fion(3,kk) = fion(3,kk) - z2*dE2

        fion(1,ll) = fion(1,ll) + x2*dE2
        fion(2,ll) = fion(2,ll) + y2*dE2
        fion(3,ll) = fion(3,ll) + z2*dE2

      end do

      return
      end



*     *************************************
*     *                                   *
*     *     pspw_qmmm_spring_bond_frag    *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_spring_bond_frag(nbond,indx,Kr0,rion)
      implicit none
      integer nbond,indx(2,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,b
      real*8  r,E,x,y,z

      E = 0.0d0
      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)
 
        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)
        E = E + Kr0(1,b)*(r-Kr0(2,b))**2
      end do

      pspw_qmmm_spring_bond_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *     pspw_qmmm_spring_bond_frag_fion    *
*     *                                        *
*     ******************************************
      subroutine pspw_qmmm_spring_bond_frag_fion(nbond,indx,Kr0,rion,
     >                                           fion)
      implicit none
      integer nbond,indx(2,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,b
      real*8  r,dE,x,y,z

      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)

        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)
        dE =  2.0d0*Kr0(1,b)*(r-Kr0(2,b))/r
c        dE =  2.0d0*Kr0(1,b)*(1.0d0-Kr0(2,b)/r)

      

        fion(1,ii) = fion(1,ii) - x*dE
        fion(2,ii) = fion(2,ii) - y*dE
        fion(3,ii) = fion(3,ii) - z*dE

        fion(1,jj) = fion(1,jj) + x*dE
        fion(2,jj) = fion(2,jj) + y*dE
        fion(3,jj) = fion(3,jj) + z*dE
      end do

      return
      end



*     *************************************
*     *                                   *
*     *     pspw_qmmm_spring_angle_frag   *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_spring_angle_frag(nangle,indx,Kr0,rion)
      implicit none
      integer nangle,indx(3,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,kk,a
      real*8  x1,y1,z1,r1,x2,y2,z2,r2,ctheta,denom,q,E

      E = 0.0d0
      do a=1,nangle
        ii = indx(1,a)
        jj = indx(2,a)
        kk = indx(3,a)
        x1 = rion(1,ii) - rion(1,jj)
        y1 = rion(2,ii) - rion(2,jj)
        z1 = rion(3,ii) - rion(3,jj)
        call lattice_min_difference(x1,y1,z1)
        r1 = dsqrt(x1**2 + y1**2 + z1**2)

        x2 = rion(1,kk) - rion(1,jj)
        y2 = rion(2,kk) - rion(2,jj)
        z2 = rion(3,kk) - rion(3,jj)
        call lattice_min_difference(x2,y2,z2)
        r2 = dsqrt(x2**2 + y2**2 + z2**2)
        denom = r1*r2
        if (denom.gt.1.0d-11) then
           ctheta = (x1*x2+y1*y2+z1*z2)/(denom)
           if (ctheta .gt.  1.0d0) ctheta = 1.0d0
           if (ctheta .lt. -1.0d0) ctheta = -1.0d0
           q      = dacos(ctheta)-Kr0(2,a)
           E      = E + Kr0(1,a)*q*q
        end if
      end do

      pspw_qmmm_spring_angle_frag = E
      return
      end

      

*     ******************************************
*     *                                        *
*     *     pspw_qmmm_spring_angle_frag_fion   *
*     *                                        *
*     ******************************************

      subroutine pspw_qmmm_spring_angle_frag_fion(nangle,indx,Kr0,
     >                                            rion,fion)
      implicit none
      integer nangle,indx(3,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,kk,a
      real*8  x1,y1,z1,r1,x2,y2,z2,r2,ctheta,denom,q,E
      real*8  r1sq,r2sq,stheta,tk,a11,a12,a22,aa
      real*8  vx1,vy1,vz1
      real*8  vx2,vy2,vz2

      E = 0.0d0
      do a=1,nangle
        ii = indx(1,a)
        jj = indx(2,a)
        kk = indx(3,a)
        x1 = rion(1,ii) - rion(1,jj)
        y1 = rion(2,ii) - rion(2,jj)
        z1 = rion(3,ii) - rion(3,jj)
        call lattice_min_difference(x1,y1,z1)
        r1sq = (x1**2 + y1**2 + z1**2)
        r1 = dsqrt(r1sq)

        x2 = rion(1,kk) - rion(1,jj)
        y2 = rion(2,kk) - rion(2,jj)
        z2 = rion(3,kk) - rion(3,jj)
        call lattice_min_difference(x2,y2,z2)
        r2sq = (x2**2 + y2**2 + z2**2)
        r2 = dsqrt(r2sq)

        denom = r1*r2
        if (denom.gt.1.0d-11) then
           ctheta = (x1*x2+y1*y2+z1*z2)/(denom)
           if (ctheta .gt.  1.0d0) ctheta =  1.0d0
           if (ctheta .lt. -1.0d0) ctheta = -1.0d0
           stheta = dsqrt(1.0d0 - ctheta*ctheta)
           if (stheta .lt. 0.001d0) stheta = 0.001d0
           stheta = 1.0d0/stheta

           q      = dacos(ctheta) - Kr0(2,a)
           tk     = Kr0(1,a)*q
           aa     = 2.0d0*tk*stheta
           a11    =  aa*ctheta/r1sq
           a12    = -aa/(denom)
           a22    =  aa*ctheta/r2sq

           vx1 = a11*x1 + a12*x2
           vx2 = a22*x2 + a12*x1

           vy1 = a11*y1 + a12*y2
           vy2 = a22*y2 + a12*y1

           vz1 = a11*z1 + a12*z2
           vz2 = a22*z2 + a12*z1

           fion(1,ii) = fion(1,ii) - vx1
           fion(2,ii) = fion(2,ii) - vy1
           fion(3,ii) = fion(3,ii) - vz1

           fion(1,jj) = fion(1,jj) + vx1 + vx2
           fion(2,jj) = fion(2,jj) + vy1 + vy2
           fion(3,jj) = fion(3,jj) + vz1 + vz2

           fion(1,kk) = fion(1,kk) - vx2
           fion(2,kk) = fion(2,kk) - vy2
           fion(3,kk) = fion(3,kk) - vz2
        end if
      end do

      return
      end

*     *************************************
*     *                                   *
*     *    pspw_qmmm_morse_bond_frag    *
*     *                                   *
*     *************************************

* Pair Ion-Ion parameters - 1 - 4*Aij*((Bij/rij)**12-(Bij/rij)**6),
*                           2 - Aij*exp(-rij/Bij) - Cij/rij**6,
*                           3 - Aij*exp(-rij/Bij) - Cij/rij**6 - Dij/rij**8
*                           4 - Aij*(rij-Bij)^Cij,
*                           5 - Aij*exp(-rij/Bij)
*                           6 - Aij*(1-exp(-Bij*(rij-Cij)))**2
*                           7 - ???
*                           8 - Aij*exp(-(rij/Bij)**2)
*                           9 - Aij*exp(-(rij/Bij)**Cij)

      real*8 function pspw_qmmm_morse_bond_frag(nbond,indx,p,rion)
      implicit none
      integer nbond,indx(3,*)
      real*8  p(6,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,b,ptype
      real*8  r,x,y,z,u,u6,u12
      real*8  E

*     **** external functions ****
      real*8   pspw_VPair_fdamp
      external pspw_VPair_fdamp

      E = 0.0d0
      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)
        ptype = indx(3,b)
 
        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)

        if (ptype.eq.1) then
            u   = (p(2,b)/r)
            u6  = u**6
            u12 = u6**2
            E = E + 4.0d0*p(1,b)*(u12-u6)
        else if (ptype.eq.2) then
           E = E + p(1,b)*dexp(-r/p(2,b)) - p(3,b)/r**6
        else if (ptype.eq.3) then
           E = E + p(1,b)*dexp(-r/p(2,b)) - p(3,b)/r**6 - p(4,b)/r**8
        else if (ptype.eq.4) then
           E = E + p(1,b)*( (r-p(2,b))**p(3,b) )
        else if (ptype.eq.5) then
           E = E + p(1,b)*dexp(-r/p(2,b))
        else if (ptype.eq.6) then
           !E = E + Kr0(1,b)*(1.0d0-dexp(-Kr0(2,b)*(r-Kr0(3,b))))**2
           E = E + p(1,b)*(1.0d0-dexp(-p(2,b)*(r-p(3,b))))**2
        else if (ptype.eq.7) then
            E = E + p(1,b)*dexp(-r/p(2,b))
     >            - p(3,b)/r**6 * pspw_VPair_fdamp(6,p(5,b),1.0d0,r)
     >            - p(4,b)/r**8 * pspw_VPair_fdamp(8,p(6,b),1.0d0,r)
        else if (ptype.eq.8) then
           E = E + p(1,b)*dexp(-(r/p(2,b))**2)
        else if (ptype.eq.9) then
           E = E + p(1,b)*dexp(-(r/p(2,b))**p(3,b))
        end if

      end do

      pspw_qmmm_morse_bond_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *    pspw_qmmm_morse_bond_frag_fion      *
*     *                                        *
*     ******************************************

* Pair Ion-Ion parameters - 1 - 4*Aij*((Bij/rij)**12-(Bij/rij)**6),
*                           2 - Aij*exp(-rij/Bij) - Cij/rij**6,
*                           3 - Aij*exp(-rij/Bij) - Cij/rij**6 - Dij/rij**8
*                           4 - Aij*(rij-Bij)^Cij,
*                           5 - Aij*exp(-rij/Bij)
*                           6 - Aij*(1-exp(-Bij*(rij-Cij)))**2
*                           7 - ???
*                           8 - Aij*exp(-(rij/Bij)**2)
*                           8 - Aij*exp(-(rij/Bij)**Cij)

      subroutine pspw_qmmm_morse_bond_frag_fion(nbond,indx,p,rion,fion)
      implicit none
      integer nbond,indx(3,*)
      real*8  p(6,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,b,ptype
      real*8  x,y,z,r,u,u6,u12
      real*8  dE

*     **** external functions ****
      real*8   pspw_VPair_fdamp,pspw_Vpair_dfdamp
      external pspw_VPair_fdamp,pspw_Vpair_dfdamp


      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)
        ptype = indx(3,b)

        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)


        if (ptype.eq.1) then
           u = (p(2,b)/r)
           u6  = u**6
           u12 = u6**2
           dE = -(4.0d0*p(1,b)/r)*(12.0d0*u12-6.0d0*u6)
        else if (ptype.eq.2) then
           dE= -p(1,b)/p(2,b)*dexp(-r/p(2,b))+6.0d0*p(3,b)/r**7
        else if (ptype.eq.3) then
           dE= -p(1,b)/p(2,b)*dexp(-r/p(2,b))+6.0d0*p(3,b)/r**7
     >                                       +8.0d0*p(4,b)/r**9
        else if (ptype.eq.4) then
           dE =  p(1,b)*p(3,b)*( (r-p(2,b))**(p(3,b)-1.0d0) )
        else if (ptype.eq.5) then
           dE = -p(1,b)/p(2,b)*dexp(-r/p(2,b))
        else if (ptype.eq.6) then
           dE = 2.0d0*p(1,b)*(1.0d0-dexp(-p(2,b)*(r-p(3,b))))
     >               *p(2,b)*dexp(-p(2,b)*(r-p(3,b)))
        else if (ptype.eq.7) then
           dE = -p(1,b)/p(2,b)*dexp(-r/p(2,b))
     >     + 6.0d0*p(3,b)/r**7 * pspw_VPair_fdamp(6,p(5,b),1.0d0,r)
     >     + 8.0d0*p(4,b)/r**9 * pspw_VPair_fdamp(8,p(6,b),1.0d0,r)
     >     - p(3,b)/r**6 * pspw_VPair_dfdamp(6,p(5,b),1.0d0,r)
     >     - p(4,b)/r**8 * pspw_VPair_dfdamp(8,p(6,b),1.0d0,r)
        else if (ptype.eq.8) then
           dE = -(2.0d0*r*p(1,b)/p(2,b)**2)*dexp(-(r/p(2,b))**2)
        else if (ptype.eq.9) then
           dE = -(p(1,b)*p(3,b)*(r/p(2,b))**p(3,b))
     >           *dexp(-(r/p(2,b))**p(3,b))/r
        else
           dE = 0.0d0
        end if

        fion(1,ii) = fion(1,ii) - (x/r)*dE 
        fion(2,ii) = fion(2,ii) - (y/r)*dE 
        fion(3,ii) = fion(3,ii) - (z/r)*dE 

        fion(1,jj) = fion(1,jj) + (x/r)*dE 
        fion(2,jj) = fion(2,jj) + (y/r)*dE 
        fion(3,jj) = fion(3,jj) + (z/r)*dE 

      end do

      return
      end


*     *************************************
*     *                                   *
*     *  pspw_qmmm_dihedral_value         *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_dihedral_value(indx,rion)
      implicit none
      integer indx(4)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  dx_ij,dy_ij,dz_ij
      real*8  dx_ik,dy_ik,dz_ik
      real*8  dx_jl,dy_jl,dz_jl
      real*8  dx_kj,dy_kj,dz_kj
      real*8  dx_kl,dy_kl,dz_kl
      real*8  dx_ijk,dy_ijk,dz_ijk
      real*8  dx_kjl,dy_kjl,dz_kjl
      real*8  dx_ijkl,dy_ijkl,dz_ijkl
      real*8  rr_ijk,rr_kjl
      real*8  E,cosphi,phi,sigphi

      ii = indx(1)
      jj = indx(2)
      kk = indx(3)
      ll = indx(4)

      !**** b1 ***
      dx_ij = rion(1,ii) - rion(1,jj)
      dy_ij = rion(2,ii) - rion(2,jj)
      dz_ij = rion(3,ii) - rion(3,jj)

      dx_ik = rion(1,ii) - rion(1,kk)
      dy_ik = rion(2,ii) - rion(2,kk)
      dz_ik = rion(3,ii) - rion(3,kk)

      dx_jl = rion(1,jj) - rion(1,ll)
      dy_jl = rion(2,jj) - rion(2,ll)
      dz_jl = rion(3,jj) - rion(3,ll)

      !**** b2 ***
      dx_kj = rion(1,kk) - rion(1,jj)
      dy_kj = rion(2,kk) - rion(2,jj)
      dz_kj = rion(3,kk) - rion(3,jj)

      !**** b3 ***
      dx_kl = rion(1,kk) - rion(1,ll)
      dy_kl = rion(2,kk) - rion(2,ll)
      dz_kl = rion(3,kk) - rion(3,ll)

      call lattice_min_difference(dx_ij,dy_ij,dz_ij)
      call lattice_min_difference(dx_ik,dy_ik,dz_ik)
      call lattice_min_difference(dx_jl,dy_jl,dz_jl)
      call lattice_min_difference(dx_kj,dy_kj,dz_kj)
      call lattice_min_difference(dx_kl,dy_kl,dz_kl)

      !*** n1=b1xb2
      dx_ijk = dy_ij * dz_kj - dy_kj * dz_ij
      dy_ijk = dz_ij * dx_kj - dz_kj * dx_ij
      dz_ijk = dx_ij * dy_kj - dx_kj * dy_ij

      !*** n2=b2xb3
      dx_kjl = dy_kj * dz_kl - dy_kl * dz_kj
      dy_kjl = dz_kj * dx_kl - dz_kl * dx_kj
      dz_kjl = dx_kj * dy_kl - dx_kl * dy_kj

      !**** n1xn2 
      dx_ijkl = dy_ijk * dz_kjl - dy_kjl * dz_ijk
      dy_ijkl = dz_ijk * dx_kjl - dz_kjl * dx_ijk
      dz_ijkl = dx_ijk * dy_kjl - dx_kjl * dy_ijk

      rr_ijk = 1.0d0 / dsqrt( dx_ijk * dx_ijk  +
     >                        dy_ijk * dy_ijk  +
     >                        dz_ijk * dz_ijk )

      rr_kjl = 1.0d0 / dsqrt( dx_kjl * dx_kjl  +
     >                        dy_kjl * dy_kjl  +
     >                        dz_kjl * dz_kjl )

      !**** n1*n2/(|n1||n2|)
      cosphi =  rr_ijk * rr_kjl * ( dx_ijk * dx_kjl  +
     >                              dy_ijk * dy_kjl  +
     >                              dz_ijk * dz_kjl )

      sigphi = dsign ( 1.0d0 , (dx_kj * dx_ijkl  +
     >                          dy_kj * dy_ijkl  +
     >                          dz_kj * dz_ijkl) )


      if (cosphi.gt.1.0d0)  cosphi =  1.0d0
      if (cosphi.lt.-1.0d0) cosphi = -1.0d0

      !phi    = -sigphi*dacos(cosphi)
      phi    =         dacos(cosphi)

      pspw_qmmm_dihedral_value = phi
      return
      end




*     *************************************
*     *                                   *
*     *  pspw_qmmm_dihedral_spring_frag   *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_dihedral_spring_frag(nangle,indx,Kr0,
     >                                               rion)
      implicit none
      integer nangle,indx(4,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  dx_ij,dy_ij,dz_ij
      real*8  dx_ik,dy_ik,dz_ik
      real*8  dx_jl,dy_jl,dz_jl
      real*8  dx_kj,dy_kj,dz_kj
      real*8  dx_kl,dy_kl,dz_kl
      real*8  dx_ijk,dy_ijk,dz_ijk
      real*8  dx_kjl,dy_kjl,dz_kjl
      real*8  dx_ijkl,dy_ijkl,dz_ijkl
      real*8  rr_ijk,rr_kjl
      real*8  E,cosphi,phi,sgn,sigphi

      E = 0.0d0
      do b=1,nangle
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)

        dx_ij = rion(1,ii) - rion(1,jj)
        dy_ij = rion(2,ii) - rion(2,jj)
        dz_ij = rion(3,ii) - rion(3,jj)

        dx_ik = rion(1,ii) - rion(1,kk)
        dy_ik = rion(2,ii) - rion(2,kk)
        dz_ik = rion(3,ii) - rion(3,kk)

        dx_jl = rion(1,jj) - rion(1,ll)
        dy_jl = rion(2,jj) - rion(2,ll)
        dz_jl = rion(3,jj) - rion(3,ll)

        dx_kj = rion(1,kk) - rion(1,jj)
        dy_kj = rion(2,kk) - rion(2,jj)
        dz_kj = rion(3,kk) - rion(3,jj)

        dx_kl = rion(1,kk) - rion(1,ll)
        dy_kl = rion(2,kk) - rion(2,ll)
        dz_kl = rion(3,kk) - rion(3,ll)

        call lattice_min_difference(dx_ij,dy_ij,dz_ij)
        call lattice_min_difference(dx_ik,dy_ik,dz_ik)
        call lattice_min_difference(dx_jl,dy_jl,dz_jl)
        call lattice_min_difference(dx_kj,dy_kj,dz_kj)
        call lattice_min_difference(dx_kl,dy_kl,dz_kl)

        dx_ijk = dy_ij * dz_kj - dy_kj * dz_ij
        dy_ijk = dz_ij * dx_kj - dz_kj * dx_ij
        dz_ijk = dx_ij * dy_kj - dx_kj * dy_ij

        dx_kjl = dy_kj * dz_kl - dy_kl * dz_kj
        dy_kjl = dz_kj * dx_kl - dz_kl * dx_kj
        dz_kjl = dx_kj * dy_kl - dx_kl * dy_kj

        dx_ijkl = dy_ijk * dz_kjl - dy_kjl * dz_ijk
        dy_ijkl = dz_ijk * dx_kjl - dz_kjl * dx_ijk
        dz_ijkl = dx_ijk * dy_kjl - dx_kjl * dy_ijk


        rr_ijk = 1.0d0 / dsqrt( dx_ijk * dx_ijk  +
     >                          dy_ijk * dy_ijk  +
     >                          dz_ijk * dz_ijk )

        rr_kjl = 1.0d0 / dsqrt( dx_kjl * dx_kjl  +
     >                          dy_kjl * dy_kjl  +
     >                          dz_kjl * dz_kjl )


        cosphi =  rr_ijk * rr_kjl * ( dx_ijk * dx_kjl  +
     >                                dy_ijk * dy_kjl  +
     >                                dz_ijk * dz_kjl )

        sigphi = dsign ( 1.0d0 , (dx_kj * dx_ijkl  +
     >                            dy_kj * dy_ijkl  +
     >                            dz_kj * dz_ijkl) )

        if (cosphi.gt.1.0d0)  cosphi =  1.0d0
        if (cosphi.lt.-1.0d0) cosphi = -1.0d0

        !phi = -sigphi*dacos(cosphi)
        phi = dacos(cosphi)
        E = E + Kr0(1,b)*(phi-Kr0(2,b))**2
      end do

      pspw_qmmm_dihedral_spring_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *  pspw_qmmm_dihedral_spring_frag_fion   *
*     *                                        *
*     ******************************************
      subroutine pspw_qmmm_dihedral_spring_frag_fion(nangle,indx,Kr0,
     >                                               rion,fion)
      implicit none
      integer nangle,indx(4,*)
      real*8  Kr0(2,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  dx_ij,dy_ij,dz_ij
      real*8  dx_ik,dy_ik,dz_ik
      real*8  dx_jl,dy_jl,dz_jl
      real*8  dx_kj,dy_kj,dz_kj
      real*8  dx_kl,dy_kl,dz_kl
      real*8  dx_ijk,dy_ijk,dz_ijk
      real*8  dx_kjl,dy_kjl,dz_kjl
      real*8  dx_ijkl,dy_ijkl,dz_ijkl
      real*8  rr_ijk,rr_kjl
      real*8  dE,cosphi,sinphi,sigphi,phi
      real*8  fxi,fyi,fzi
      real*8  fxj,fyj,fzj
      real*8  fxk,fyk,fzk
      real*8  fxl,fyl,fzl
      real*8  dx_a,dy_a,dz_a
      real*8  dx_b,dy_b,dz_b
c      real*8   dmax
c      external dmax


      do b=1,nangle
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)

        dx_ij = rion(1,ii) - rion(1,jj)
        dy_ij = rion(2,ii) - rion(2,jj)
        dz_ij = rion(3,ii) - rion(3,jj)

        dx_ik = rion(1,ii) - rion(1,kk)
        dy_ik = rion(2,ii) - rion(2,kk)
        dz_ik = rion(3,ii) - rion(3,kk)

        dx_jl = rion(1,jj) - rion(1,ll) 
        dy_jl = rion(2,jj) - rion(2,ll)
        dz_jl = rion(3,jj) - rion(3,ll)

        dx_kj = rion(1,kk) - rion(1,jj)
        dy_kj = rion(2,kk) - rion(2,jj)
        dz_kj = rion(3,kk) - rion(3,jj)

        dx_kl = rion(1,kk) - rion(1,ll)
        dy_kl = rion(2,kk) - rion(2,ll)
        dz_kl = rion(3,kk) - rion(3,ll)

        call lattice_min_difference(dx_ij,dy_ij,dz_ij)
        call lattice_min_difference(dx_ik,dy_ik,dz_ik)
        call lattice_min_difference(dx_jl,dy_jl,dz_jl)
        call lattice_min_difference(dx_kj,dy_kj,dz_kj)
        call lattice_min_difference(dx_kl,dy_kl,dz_kl)

        dx_ijk = dy_ij * dz_kj - dy_kj * dz_ij
        dy_ijk = dz_ij * dx_kj - dz_kj * dx_ij
        dz_ijk = dx_ij * dy_kj - dx_kj * dy_ij

        dx_kjl = dy_kj * dz_kl - dy_kl * dz_kj
        dy_kjl = dz_kj * dx_kl - dz_kl * dx_kj
        dz_kjl = dx_kj * dy_kl - dx_kl * dy_kj

        dx_ijkl = dy_ijk * dz_kjl - dy_kjl * dz_ijk
        dy_ijkl = dz_ijk * dx_kjl - dz_kjl * dx_ijk
        dz_ijkl = dx_ijk * dy_kjl - dx_kjl * dy_ijk

        rr_ijk = 1.0d0 / dsqrt( dx_ijk * dx_ijk  +
     >                          dy_ijk * dy_ijk  +
     >                          dz_ijk * dz_ijk )

        rr_kjl = 1.0d0 / dsqrt( dx_kjl * dx_kjl  +
     >                          dy_kjl * dy_kjl  +
     >                          dz_kjl * dz_kjl )


        cosphi =  rr_ijk * rr_kjl * ( dx_ijk * dx_kjl  +
     >                                dy_ijk * dy_kjl  +
     >                                dz_ijk * dz_kjl )


        sigphi = dsign ( 1.0d0 , (dx_kj * dx_ijkl  +
     >                            dy_kj * dy_ijkl  +
     >                            dz_kj * dz_ijkl) )


        if (cosphi.gt.1.0d0)  cosphi =  1.0d0
        if (cosphi.lt.-1.0d0) cosphi = -1.0d0

        !phi    = -sigphi*dacos(cosphi)
        phi    = dacos(cosphi)
        sinphi = dsin(phi)
        dE     = 2*Kr0(1,b) * (phi-Kr0(2,b))/max(1.0d-8,sinphi)

        dx_ijk = dx_ijk * rr_ijk
        dy_ijk = dy_ijk * rr_ijk
        dz_ijk = dz_ijk * rr_ijk

        dx_kjl = dx_kjl * rr_kjl
        dy_kjl = dy_kjl * rr_kjl
        dz_kjl = dz_kjl * rr_kjl

        dx_a = dx_kjl - cosphi * dx_ijk
        dy_a = dy_kjl - cosphi * dy_ijk
        dz_a = dz_kjl - cosphi * dz_ijk

        dx_b = dx_ijk - cosphi * dx_kjl
        dy_b = dy_ijk - cosphi * dy_kjl
        dz_b = dz_ijk - cosphi * dz_kjl

        fxi = rr_ijk * ( dy_kj * dz_a - dy_a * dz_kj )
        fyi = rr_ijk * ( dz_kj * dx_a - dz_a * dx_kj )
        fzi = rr_ijk * ( dx_kj * dy_a - dx_a * dy_kj )

        fxl = rr_kjl * ( dy_kj * dz_b - dy_b * dz_kj )
        fyl = rr_kjl * ( dz_kj * dx_b - dz_b * dx_kj )
        fzl = rr_kjl * ( dx_kj * dy_b - dx_b * dy_kj )

        fxj = (rr_ijk * ( dy_ik * dz_a - dy_a * dz_ik ) -
     >         rr_kjl * ( dy_kl * dz_b - dy_b * dz_kl ) )
        fyj = (rr_ijk * ( dz_ik * dx_a - dz_a * dx_ik ) -
     >         rr_kjl * ( dz_kl * dx_b - dz_b * dx_kl ) )
        fzj = (rr_ijk * ( dx_ik * dy_a - dx_a * dy_ik ) -
     >         rr_kjl * ( dx_kl * dy_b - dx_b * dy_kl ) )

        fxk = (rr_kjl * ( dy_jl * dz_b - dy_b * dz_jl ) -
     >         rr_ijk * ( dy_ij * dz_a - dy_a * dz_ij ) )
        fyk = (rr_kjl * ( dz_jl * dx_b - dz_b * dx_jl ) -
     >         rr_ijk * ( dz_ij * dx_a - dz_a * dx_ij ) )
        fzk = (rr_kjl * ( dx_jl * dy_b - dx_b * dy_jl ) -
     >         rr_ijk * ( dx_ij * dy_a - dx_a * dy_ij ) )

        fion(1,ii) = fion(1,ii) + fxi*dE
        fion(2,ii) = fion(2,ii) + fyi*dE
        fion(3,ii) = fion(3,ii) + fzi*dE 

        fion(1,jj) = fion(1,jj) + fxj*dE
        fion(2,jj) = fion(2,jj) + fyj*dE
        fion(3,jj) = fion(3,jj) + fzj*dE

        fion(1,kk) = fion(1,kk) + fxk*dE
        fion(2,kk) = fion(2,kk) + fyk*dE
        fion(3,kk) = fion(3,kk) + fzk*dE

        fion(1,ll) = fion(1,ll) + fxl*dE
        fion(2,ll) = fion(2,ll) + fyl*dE
        fion(3,ll) = fion(3,ll) + fzl*dE

      end do

      return
      end



*     *************************************
*     *                                   *
*     *     pspw_qmmm_dihedral_frag       *
*     *                                   *
*     *************************************

      real*8 function pspw_qmmm_dihedral_frag(nangle,indx,Kr0,rion)
      implicit none
      integer nangle,indx(4,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  dx_ij,dy_ij,dz_ij
      real*8  dx_ik,dy_ik,dz_ik
      real*8  dx_jl,dy_jl,dz_jl
      real*8  dx_kj,dy_kj,dz_kj
      real*8  dx_kl,dy_kl,dz_kl
      real*8  dx_ijk,dy_ijk,dz_ijk
      real*8  dx_kjl,dy_kjl,dz_kjl
      real*8  dx_ijkl,dy_ijkl,dz_ijkl
      real*8  rr_ijk,rr_kjl
      real*8  E,cosphi,phi 


      E = 0.0d0
      do b=1,nangle
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)
        dx_ij = rion(1,ii) - rion(1,jj)
        dy_ij = rion(2,ii) - rion(2,jj)
        dz_ij = rion(3,ii) - rion(3,jj)

        dx_ik = rion(1,ii) - rion(1,kk)
        dy_ik = rion(2,ii) - rion(2,kk)
        dz_ik = rion(3,ii) - rion(3,kk)

        dx_jl = rion(1,jj) - rion(1,ll)
        dy_jl = rion(2,jj) - rion(2,ll)
        dz_jl = rion(3,jj) - rion(3,ll)

        dx_kj = rion(1,kk) - rion(1,jj)   
        dy_kj = rion(2,kk) - rion(2,jj)
        dz_kj = rion(3,kk) - rion(3,jj)

        dx_kl = rion(1,kk) - rion(1,ll)
        dy_kl = rion(2,kk) - rion(2,ll)
        dz_kl = rion(3,kk) - rion(3,ll)

        call lattice_min_difference(dx_ij,dy_ij,dz_ij)
        call lattice_min_difference(dx_ik,dy_ik,dz_ik)
        call lattice_min_difference(dx_jl,dy_jl,dz_jl)
        call lattice_min_difference(dx_kj,dy_kj,dz_kj)
        call lattice_min_difference(dx_kl,dy_kl,dz_kl)

        dx_ijk = dy_ij * dz_kj - dy_kj * dz_ij
        dy_ijk = dz_ij * dx_kj - dz_kj * dx_ij
        dz_ijk = dx_ij * dy_kj - dx_kj * dy_ij

        dx_kjl = dy_kj * dz_kl - dy_kl * dz_kj
        dy_kjl = dz_kj * dx_kl - dz_kl * dx_kj
        dz_kjl = dx_kj * dy_kl - dx_kl * dy_kj

        dx_ijkl = dy_ijk * dz_kjl - dy_kjl * dz_ijk
        dy_ijkl = dz_ijk * dx_kjl - dz_kjl * dx_ijk
        dz_ijkl = dx_ijk * dy_kjl - dx_kjl * dy_ijk

        rr_ijk = 1.0d0 / dsqrt( dx_ijk * dx_ijk  +
     >                          dy_ijk * dy_ijk  +
     >                          dz_ijk * dz_ijk )
        rr_kjl = 1.0d0 / dsqrt( dx_kjl * dx_kjl  +
     >                          dy_kjl * dy_kjl  +
     >                          dz_kjl * dz_kjl )
        cosphi =  rr_ijk * rr_kjl * ( dx_ijk * dx_kjl  +
     >                                dy_ijk * dy_kjl  +
     >                                dz_ijk * dz_kjl )

        if (cosphi.gt.1.0d0)  cosphi =  1.0d0
        if (cosphi.lt.-1.0d0) cosphi = -1.0d0
        phi    = dacos(cosphi)
        
        E = E + Kr0(1,b)*(1.0d0+dcos(Kr0(2,b)*phi-Kr0(3,b)))
      end do

      pspw_qmmm_dihedral_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *      pspw_qmmm_dihedral_frag_fion      *
*     *                                        *
*     ******************************************
      subroutine pspw_qmmm_dihedral_frag_fion(nangle,indx,Kr0,
     >                                        rion,fion)
      implicit none
      integer nangle,indx(4,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,kk,ll,b
      real*8  dx_ij,dy_ij,dz_ij
      real*8  dx_ik,dy_ik,dz_ik
      real*8  dx_jl,dy_jl,dz_jl
      real*8  dx_kj,dy_kj,dz_kj
      real*8  dx_kl,dy_kl,dz_kl
      real*8  dx_ijk,dy_ijk,dz_ijk
      real*8  dx_kjl,dy_kjl,dz_kjl
      real*8  dx_ijkl,dy_ijkl,dz_ijkl
      real*8  rr_ijk,rr_kjl
      real*8  dE,cosphi,sinphi,sigphi,phi
      real*8  fxi,fyi,fzi
      real*8  fxj,fyj,fzj
      real*8  fxk,fyk,fzk
      real*8  fxl,fyl,fzl
      real*8  dx_a,dy_a,dz_a
      real*8  dx_b,dy_b,dz_b
c      real*8   dmax
c      external dmax

      do b=1,nangle
        ii = indx(1,b)
        jj = indx(2,b)
        kk = indx(3,b)
        ll = indx(4,b)
        dx_ij = rion(1,ii) - rion(1,jj)
        dy_ij = rion(2,ii) - rion(2,jj)
        dz_ij = rion(3,ii) - rion(3,jj)

        dx_ik = rion(1,ii) - rion(1,kk)
        dy_ik = rion(2,ii) - rion(2,kk)
        dz_ik = rion(3,ii) - rion(3,kk)

        dx_jl = rion(1,jj) - rion(1,ll)
        dy_jl = rion(2,jj) - rion(2,ll)
        dz_jl = rion(3,jj) - rion(3,ll)

        dx_kj = rion(1,kk) - rion(1,jj)
        dy_kj = rion(2,kk) - rion(2,jj)
        dz_kj = rion(3,kk) - rion(3,jj)

        dx_kl = rion(1,kk) - rion(1,ll)
        dy_kl = rion(2,kk) - rion(2,ll)
        dz_kl = rion(3,kk) - rion(3,ll)

        call lattice_min_difference(dx_ij,dy_ij,dz_ij)
        call lattice_min_difference(dx_ik,dy_ik,dz_ik)
        call lattice_min_difference(dx_jl,dy_jl,dz_jl)
        call lattice_min_difference(dx_kj,dy_kj,dz_kj)
        call lattice_min_difference(dx_kl,dy_kl,dz_kl)

        dx_ijk = dy_ij * dz_kj - dy_kj * dz_ij
        dy_ijk = dz_ij * dx_kj - dz_kj * dx_ij
        dz_ijk = dx_ij * dy_kj - dx_kj * dy_ij

        dx_kjl = dy_kj * dz_kl - dy_kl * dz_kj
        dy_kjl = dz_kj * dx_kl - dz_kl * dx_kj
        dz_kjl = dx_kj * dy_kl - dx_kl * dy_kj

        dx_ijkl = dy_ijk * dz_kjl - dy_kjl * dz_ijk
        dy_ijkl = dz_ijk * dx_kjl - dz_kjl * dx_ijk
        dz_ijkl = dx_ijk * dy_kjl - dx_kjl * dy_ijk

        rr_ijk = 1.0d0 / dsqrt( dx_ijk * dx_ijk  +
     >                          dy_ijk * dy_ijk  +
     >                          dz_ijk * dz_ijk )

        rr_kjl = 1.0d0 / dsqrt( dx_kjl * dx_kjl  +
     >                          dy_kjl * dy_kjl  +
     >                          dz_kjl * dz_kjl )
        cosphi =  rr_ijk * rr_kjl * ( dx_ijk * dx_kjl  +
     >                                dy_ijk * dy_kjl  +
     >                                dz_ijk * dz_kjl )


        sigphi = dsign ( 1.0d0 , (dx_kj * dx_ijkl  +
     >                             dy_kj * dy_ijkl  +
     >                             dz_kj * dz_ijkl) )

        if (cosphi.gt.1.0d0)  cosphi =  1.0d0
        if (cosphi.lt.-1.0d0) cosphi = -1.0d0

        phi    = dacos(cosphi)
        sinphi = dsin(phi)

cccc  E =  Kr0(1,b)*(1.0d0+dcos(Kr0(2,b)*phi-Kr0(3,b))) cccc

        dE =  -Kr0(1,b)*(Kr0(2,b)*dsin(Kr0(2,b)*phi-Kr0(3,b)))
        dE = dE/max(1.0d-8,sinphi)

        dx_ijk = dx_ijk * rr_ijk
        dy_ijk = dy_ijk * rr_ijk
        dz_ijk = dz_ijk * rr_ijk

        dx_kjl = dx_kjl * rr_kjl
        dy_kjl = dy_kjl * rr_kjl
        dz_kjl = dz_kjl * rr_kjl

        dx_a = dx_kjl - cosphi * dx_ijk
        dy_a = dy_kjl - cosphi * dy_ijk
        dz_a = dz_kjl - cosphi * dz_ijk

        dx_b = dx_ijk - cosphi * dx_kjl
        dy_b = dy_ijk - cosphi * dy_kjl
        dz_b = dz_ijk - cosphi * dz_kjl

        fxi = rr_ijk * ( dy_kj * dz_a - dy_a * dz_kj )
        fyi = rr_ijk * ( dz_kj * dx_a - dz_a * dx_kj )
        fzi = rr_ijk * ( dx_kj * dy_a - dx_a * dy_kj )

        fxl = rr_kjl * ( dy_kj * dz_b - dy_b * dz_kj )
        fyl = rr_kjl * ( dz_kj * dx_b - dz_b * dx_kj )
        fzl = rr_kjl * ( dx_kj * dy_b - dx_b * dy_kj )

        fxj = (rr_ijk * ( dy_ik * dz_a - dy_a * dz_ik ) -
     >         rr_kjl * ( dy_kl * dz_b - dy_b * dz_kl ) )
        fyj = (rr_ijk * ( dz_ik * dx_a - dz_a * dx_ik ) -
     >         rr_kjl * ( dz_kl * dx_b - dz_b * dx_kl ) )
        fzj = (rr_ijk * ( dx_ik * dy_a - dx_a * dy_ik ) -
     >         rr_kjl * ( dx_kl * dy_b - dx_b * dy_kl ) )

        fxk = (rr_kjl * ( dy_jl * dz_b - dy_b * dz_jl ) -
     >         rr_ijk * ( dy_ij * dz_a - dy_a * dz_ij ) )
        fyk = (rr_kjl * ( dz_jl * dx_b - dz_b * dx_jl ) -
     >         rr_ijk * ( dz_ij * dx_a - dz_a * dx_ij ) )
        fzk = (rr_kjl * ( dx_jl * dy_b - dx_b * dy_jl ) -
     >         rr_ijk * ( dx_ij * dy_a - dx_a * dy_ij ) )
      end do

      return
      end



*     *************************************
*     *                                   *
*     *   pspw_qmmm_spring_coordnum_frag  *
*     *                                   *
*     *************************************

*   This routine returns the coordination number spring energy.

      real*8 function pspw_qmmm_spring_coordnum_frag(ncoords,
     >                                           nsize1,indx1,
     >                                           nsize2,indx2,
     >                                            Kr0,rion)
      implicit none
      integer ncoords
      integer nsize1(*),indx1(*)
      integer nsize2(*),indx2(*)
      real*8  Kr0(6,*)
      real*8  rion(3,*)

*     **** local variables ****
      logical sprik
      integer i,j,ii,jj,a1,a2,c,natom1,natom2
      real*8  r,E,x,y,z,r0,f,n,m

      E = 0.0d0
      ii = 0 
      jj = 0 
      do c=1,ncoords
         natom1 = nsize1(c)
         natom2 = nsize2(c)
         n      = Kr0(1,c)
         m      = Kr0(2,c)
         r0     = Kr0(3,c)
         sprik  = (Kr0(6,c).gt.0.0d0)
      
         f = 0.0d0
         do a1=1,natom1
            do a2=1,natom2
              i = indx1(ii+a1) 
              j = indx2(jj+a2)
              x = rion(1,i)-rion(1,j)
              y = rion(2,i)-rion(2,j)
              z = rion(3,i)-rion(3,j)
              call lattice_min_difference(x,y,z)
              r = dsqrt( x**2 + y**2 + z**2)
              if (sprik) then
                 f = f + 1.d0/(1.d0+dexp(n*(r-r0)))
              else
                 f = f + (1.0d0-(r/r0)**n)/(1.0d0-(r/r0)**m)
              end if
            end do
         end do
         E = E + Kr0(4,c)*(f-Kr0(5,c))**2
         ii = ii + natom1
         jj = jj + natom2
      end do

      pspw_qmmm_spring_coordnum_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *  pspw_qmmm_spring_coordnum_frag_fion   *
*     *                                        *
*     ******************************************

*   This routine calculates the coordination number spring force 
*
*   !!!!This routine can be sped up.
*
      subroutine pspw_qmmm_spring_coordnum_frag_fion(ncoords,
     >                                           nsize1,indx1,
     >                                           nsize2,indx2,
     >                                           Kr0,rion,
     >                                           fion)
      implicit none
      integer ncoords
      integer nsize1(*),indx1(*)
      integer nsize2(*),indx2(*)
      real*8  Kr0(6,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      logical sprik
      integer i,j,ii,jj,c,natom1,natom2,a1,a2
      real*8  r,dE,x,y,z,n,m,f,df,r0,rn,rm

      ii = 0 
      jj = 0 
      do c=1,ncoords
         natom1 = nsize1(c)
         natom2 = nsize2(c)
         n      = Kr0(1,c)
         m      = Kr0(2,c)
         r0     = Kr0(3,c)
         sprik  = (Kr0(6,c).gt.0.0d0)

         !*** calculate f ****
         f = 0.0d0
         do a1=1,natom1
            do a2=1,natom2
              i = indx1(ii+a1) 
              j = indx2(jj+a2)
              x = rion(1,i)-rion(1,j)
              y = rion(2,i)-rion(2,j)
              z = rion(3,i)-rion(3,j)
              call lattice_min_difference(x,y,z)
              r = dsqrt( x**2 + y**2 + z**2)
              if (sprik) then
                 f = f + 1.0d0/(1.0d0+dexp(n*(r-r0)))
              else
                 f = f + (1.0d0-(r/r0)**n)/(1.0d0-(r/r0)**m)
              end if
            end do
         end do
         dE =  2.0d0*Kr0(4,c)*(f-Kr0(5,c))

         do a1=1,natom1
            do a2=1,natom2
               i = indx1(ii+a1) 
               j = indx2(jj+a2)
               x = rion(1,i)-rion(1,j)
               y = rion(2,i)-rion(2,j)
               z = rion(3,i)-rion(3,j)
               call lattice_min_difference(x,y,z)
               r = dsqrt( x**2 + y**2 + z**2)
               if (sprik) then
                  rn  = 1.0d0+dexp(n*(r-r0))
                  df = -n*(rn-1.0d0)/(rn*rn)
               else
                  rn  = (1.0d0-(r/r0)**n)
                  rm  = (1.0d0-(r/r0)**m)
                  df = (-n*rm/r0*(r/r0)**(n-1) + m*rn/r0*(r/r0)**(m-1))
     >                 / (rm**2)
               end if
 
               fion(1,i) = fion(1,i) - (x/r)*df*dE
               fion(2,i) = fion(2,i) - (y/r)*df*dE
               fion(3,i) = fion(3,i) - (z/r)*df*dE

               fion(1,j) = fion(1,j) + (x/r)*df*dE
               fion(2,j) = fion(2,j) + (y/r)*df*dE
               fion(3,j) = fion(3,j) + (z/r)*df*dE
            end do
         end do
         ii = ii + natom1
         jj = jj + natom2
      end do

      return
      end


*     *************************************
*     *                                   *
*     *   pspw_qmmm_coordnum              *
*     *                                   *
*     *************************************

*   This routine returns the coordination number spring energy.

      real*8 function pspw_qmmm_coordnum(nsize1,indx1,
     >                                   nsize2,indx2,
     >                                   Kr0,rion)
      implicit none
      integer ncoords
      integer nsize1,indx1(*)
      integer nsize2,indx2(*)
      real*8  Kr0(6)
      real*8  rion(3,*)

*     **** local variables ****
      logical sprik
      integer i,j,a1,a2
      real*8  r,x,y,z,r0,f,n,m

      n      = Kr0(1)
      m      = Kr0(2)
      r0     = Kr0(3)
      sprik  = (Kr0(6).gt.0.0d0)
      
      f = 0.0d0
      do a1=1,nsize1
         do a2=1,nsize2
           i = indx1(a1) 
           j = indx2(a2)
           x = rion(1,i)-rion(1,j)
           y = rion(2,i)-rion(2,j)
           z = rion(3,i)-rion(3,j)
           call lattice_min_difference(x,y,z)
           r = dsqrt( x**2 + y**2 + z**2)
           if (sprik) then
              f = f + 1.d0/(1.d0+dexp(n*(r-r0)))
           else
              f = f + (1.0d0-(r/r0)**n)/(1.0d0-(r/r0)**m)
           end if
         end do
      end do

      pspw_qmmm_coordnum = f
      return
      end



*     *************************************
*     *                                   *
*     *   pspw_qmmm_spring_gpenalty_frag  *
*     *                                   *
*     *************************************

*   This routine returns the gaussian penalty spring energy.
*   Epenalty = sum(i) A*exp(-B*((xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2))

      real*8 function pspw_qmmm_spring_gpenalty_frag(npenalty,
     >                                           nsize1,indx1,
     >                                           Kr0,rion)
      implicit none
      integer npenalty
      integer nsize1(*),indx1(*)
      real*8  Kr0(5,*)
      real*8  rion(3,*)

      integer i,ii,a1,c,natom1
      real*8  r2,E,x0,y0,z0,x,y,z,A,B

      E = 0.0d0
      ii = 0
      do c=1,npenalty
         natom1 = nsize1(c)
         A      = Kr0(1,c)
         B      = Kr0(2,c)
         x0     = Kr0(3,c)
         y0     = Kr0(4,c)
         z0     = Kr0(5,c)

         do a1=1,natom1
           i = indx1(ii+a1)
           x = rion(1,i)-x0
           y = rion(2,i)-y0
           z = rion(3,i)-z0
           call lattice_min_difference(x,y,z)
           r2 =  x**2 + y**2 + z**2
           E = E + A*exp(-B*r2)
         end do
         ii = ii + natom1
      end do

      pspw_qmmm_spring_gpenalty_frag = E

      return
      end



*     *************************************
*     *                                   *
*     *   pspw_qmmm_spring_bondings_frag  *
*     *                                   *
*     *************************************

*   This routine returns the gaussian penalty spring energy.
*   Ebondings = sum(i) 

      real*8 function pspw_qmmm_spring_bondings_frag(nbondings,
     >                                           nsize1,indx1,
     >                                           Kr0,rion)
      implicit none
      integer nbondings
      integer nsize1(*),indx1(*)
      real*8  Kr0(*)
      real*8  rion(3,*)

*     **** local variables ****
      integer i,ii,jj,n1,i1,i2,c
      real*8 x1,y1,z1,x2,y2,z2,dx12,dy12,dz12
      real*8 gamm,gamm0,K0,r12,coef,E


      E = 0.0d0
      ii = 0
      jj = 0
      do c=1,nbondings
         n1    = nsize1(c)
         K0    = Kr0(ii+1)
         gamm0 = Kr0(ii+2)

         gamm = 0.0d0
         do i=1,n1
            i1 = indx1(ii+2*i-1)
            i2 = indx1(ii+2*i)
            coef = Kr0(jj+2+i)
            x1 = rion(1,i1)
            y1 = rion(2,i1)
            z1 = rion(3,i1)
            x2 = rion(1,i2)
            y2 = rion(2,i2)
            z2 = rion(3,i2)
            dx12 = x1-x2
            dy12 = y1-y2
            dz12 = z1-z2
            call lattice_min_difference(dx12,dy12,dz12)
            r12  = dsqrt(dx12**2 + dy12**2 + dz12**2)
            gamm = gamm + Kr0(ii+2+i)*r12
         end do
         ii = ii + 2*n1
         jj = jj + n1 + 2 
         E = E + K0*(gamm-gamm0)**2
      end do

      pspw_qmmm_spring_bondings_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *   pspw_qmmm_spring_bondings_frag_fion  *
*     *                                        *
*     ******************************************

*   This routine returns the gaussian penalty spring energy.
*   Ebondings = sum(i)

      subroutine pspw_qmmm_spring_bondings_frag_fion(nbondings,
     >                                           nsize1,indx1,
     >                                           Kr0,rion,fion)
      implicit none
      integer nbondings
      integer nsize1(*),indx1(*)
      real*8  Kr0(*)
      real*8  rion(3,*),fion(3,*)

*     **** local variables ****
      integer i,ii,jj,n1,c,i1,i2
      real*8 x1,y1,z1,x2,y2,z2,dx12,dy12,dz12
      real*8 gamm,gamm0,K0,r12
      real*8 E,energy,denergy,coef

      E = 0.0d0
      ii = 0
      jj = 0
      do c=1,nbondings
         n1    = nsize1(c)
         K0    = Kr0(jj+1)
         gamm0 = Kr0(jj+2)

         gamm = 0.0d0
         do i=1,n1
            i1 = indx1(ii+2*i-1)
            i2 = indx1(ii+2*i)
            coef = Kr0(jj+2+i)
            x1 = rion(1,i1)
            y1 = rion(2,i1)
            z1 = rion(3,i1)
            x2 = rion(1,i2)
            y2 = rion(2,i2)
            z2 = rion(3,i2)
            dx12 = x1-x2
            dy12 = y1-y2
            dz12 = z1-z2
            call lattice_min_difference(dx12,dy12,dz12)
            r12  = dsqrt(dx12**2 + dy12**2 + dz12**2)
            gamm = gamm + coef*r12
         end do
         energy  = K0*(gamm-gamm0)**2
         denergy = -2.0d0*K0*(gamm-gamm0)

         gamm = 0.0d0
         do i=1,n1
            i1 = indx1(ii+2*i-1)
            i2 = indx1(ii+2*i)
            coef = Kr0(jj+2+i)
            x1 = rion(1,i1)
            y1 = rion(2,i1)
            z1 = rion(3,i1)
            x2 = rion(1,i2)
            y2 = rion(2,i2)
            z2 = rion(3,i2)
            dx12 = x1-x2
            dy12 = y1-y2
            dz12 = z1-z2
            call lattice_min_difference(dx12,dy12,dz12)
            r12  = dsqrt(dx12**2 + dy12**2 + dz12**2)
            fion(1,i1) = fion(1,i1) + denergy*coef*(dx12/r12)
            fion(2,i1) = fion(2,i1) + denergy*coef*(dy12/r12)
            fion(3,i1) = fion(3,i1) + denergy*coef*(dz12/r12)
            fion(1,i2) = fion(1,i2) - denergy*coef*(dx12/r12)
            fion(2,i2) = fion(2,i2) - denergy*coef*(dy12/r12)
            fion(3,i2) = fion(3,i2) - denergy*coef*(dz12/r12)
         end do
         ii = ii + 2*n1
         jj = jj + n1 + 2
         E = E + energy
      end do

      return
      end
