c
c $Id: pspw_springs.F,v 1.1 2008-09-22 22:55:16 bylaska Exp $
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

      real*8 function pspw_qmmm_morse_bond_frag(nbond,indx,Kr0,rion)
      implicit none
      integer nbond,indx(2,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)

*     **** local variables ****
      integer ii,jj,b
      real*8  r,x,y,z
      real*8  E

      E = 0.0d0
      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)
 
        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)

        E = E + Kr0(1,b)*( 1.0d0 - dexp(-Kr0(2,b)*(r-Kr0(3,b))) )**2
      end do

      pspw_qmmm_morse_bond_frag = E
      return
      end


*     ******************************************
*     *                                        *
*     *    pspw_qmmm_morse_bond_frag_fion      *
*     *                                        *
*     ******************************************
      subroutine pspw_qmmm_morse_bond_frag_fion(nbond,indx,Kr0,rion,
     >                                           fion)
      implicit none
      integer nbond,indx(2,*)
      real*8  Kr0(3,*)
      real*8  rion(3,*)
      real*8  fion(3,*)

*     **** local variables ****
      integer ii,jj,b
      real*8  x,y,z,r
      real*8  dE,ee

      do b=1,nbond
        ii = indx(1,b)
        jj = indx(2,b)

        x = rion(1,ii)-rion(1,jj)
        y = rion(2,ii)-rion(2,jj)
        z = rion(3,ii)-rion(3,jj)
        call lattice_min_difference(x,y,z)
        r = dsqrt( x**2 + y**2 + z**2)

        ee = dexp(-Kr0(2,b)*(r-Kr0(3,b)))
        dE = 2.0d0*Kr0(1,b)*Kr0(2,b)*(1.0d0-ee)*ee/r

        fion(1,ii) = fion(1,ii) - x*dE 
        fion(2,ii) = fion(2,ii) - y*dE 
        fion(3,ii) = fion(3,ii) - z*dE 

        fion(1,jj) = fion(1,jj) + x*dE 
        fion(2,jj) = fion(2,jj) + y*dE 
        fion(3,jj) = fion(3,jj) + z*dE 

      end do

      return
      end

