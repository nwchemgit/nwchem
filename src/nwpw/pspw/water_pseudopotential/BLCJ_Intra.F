**********************************************************
*                                                        *
*              BLCJ_Intra module                         *
*                                                        *
*          BLCJ Electron-Water pseudopotential           *
*                                                        *
*          used by BLCJ_Water module                     *
*                                                        *
*    -- developed by Eric J. Bylaska on February 6, 2001 *
*                                                        *
**********************************************************



*     **********************************
*     *	                               *
*     *        BLCJ_Intramolecular     *
*     *                                *
*     **********************************

      real*8 function BLCJ_Intramolecular(ro,rh1,rh2)
      implicit none
      real*8 ro(3),rh1(3),rh2(3)

*     **** BLCJ Intramolecular parameters ****
      real*8 theta_eq,R_eq
      real*8 alpha1,alpha2,alpha3
      real*8 D1,D2,D3
      parameter (theta_eq = 1.824218134d0, R_eq=1.8088d0)
      parameter (alpha1=1.1298d0,alpha2=1.1298d0,alpha3=0.70337d0)
      parameter (    D1=0.20912d0,   D2=0.20912d0,   D3=0.15791d0)

*     **** local variables ****
      real*8 xo1,yo1,zo1,ro1
      real*8 xo2,yo2,zo2,ro2
      real*8 s1,s2,s3
      real*8 theta,cos_theta
      real*8 ewater      
   

*     **** cacluate distances between O and H's ****
      xo1 = rh1(1) - ro(1)
      yo1 = rh1(2) - ro(2)
      zo1 = rh1(3) - ro(3)
      ro1 = dsqrt(xo1*xo1 + yo1*yo1 + zo1*zo1)

      xo2 = rh2(1) - ro(1)
      yo2 = rh2(2) - ro(2)
      zo2 = rh2(3) - ro(3)
      ro2 = dsqrt(xo2*xo2 + yo2*yo2 + zo2*zo2)
      

*     **** calculate angle and then scaled distances ****
      cos_theta = (xo1*xo2 + yo1*yo2 + zo1*zo2)/(ro1*ro2)
      theta = dacos(cos_theta)
      theta = 0.5d0*(theta-theta_eq)

      s1 = ro1*dcos(theta) - R_eq
      s2 = ro2*dcos(theta) - R_eq
      s3 = (ro1+ro2)*dsin(theta)/R_eq

*     **** calculate Energy ****
      ewater =  D1 * (1.0d0 - dexp(-alpha1*s1))**2
     >       +  D2 * (1.0d0 - dexp(-alpha2*s2))**2
     >       +  D3 * (1.0d0 - dexp(-alpha3*s3))**2


      BLCJ_Intramolecular = ewater
      return
      end


*     **********************************
*     *	                               *
*     *   BLCJ_Intramolecular_Fwater   *
*     *                                *
*     **********************************

      subroutine BLCJ_Intramolecular_Fwater(ro,rh1,rh2,
     >                                      fo,f1, f2)
      implicit none
      real*8 ro(3),rh1(3),rh2(3)
      real*8 fo(3),f1(3), f2(3)    

*     **** BLCJ Intramolecular parameters ****
      real*8 theta_eq,R_eq
      real*8 alpha1,alpha2,alpha3
      real*8 D1,D2,D3
      parameter (theta_eq = 1.824218134d0, R_eq=1.8088d0)
      parameter (alpha1=1.1298d0,alpha2=1.1298d0,alpha3=0.70337d0)
      parameter (    D1=0.20912d0,   D2=0.20912d0,   D3=0.15791d0)

*     **** local variables ****
      real*8 xo1,yo1,zo1,ro1
      real*8 xo2,yo2,zo2,ro2
      real*8 s1,s2,s3
      real*8 theta,cos_theta
      real*8 x,ct,cx,st,sx
      real*8 t_ct
      real*8 ro1_xo,ro1_yo,ro1_zo
      real*8 ro1_x1,ro1_y1,ro1_z1
      real*8 ro2_xo,ro2_yo,ro2_zo
      real*8 ro2_x2,ro2_y2,ro2_z2
      real*8 ct_xo,ct_yo,ct_zo
      real*8 ct_x1,ct_y1,ct_z1
      real*8 ct_x2,ct_y2,ct_z2
      real*8 s1_ro1,s1_t
      real*8 s2_ro2,s2_t
      real*8 s3_ro1,s3_ro2,s3_t     
      real*8 E_s1,E_s2,E_s3
      real*8 E_xo,E_yo,E_zo
      real*8 E_x1,E_y1,E_z1
      real*8 E_x2,E_y2,E_z2

*     **** cacluate distances between O and H's ****
      xo1 = rh1(1) - ro(1)
      yo1 = rh1(2) - ro(2)
      zo1 = rh1(3) - ro(3)
      ro1 = dsqrt(xo1*xo1 + yo1*yo1 + zo1*zo1)
      ro1_xo = -xo1/ro1
      ro1_yo = -yo1/ro1
      ro1_zo = -zo1/ro1
      ro1_x1 =  xo1/ro1
      ro1_y1 =  yo1/ro1
      ro1_z1 =  zo1/ro1

      xo2 = rh2(1) - ro(1)
      yo2 = rh2(2) - ro(2)
      zo2 = rh2(3) - ro(3)
      ro2 = dsqrt(xo2*xo2 + yo2*yo2 + zo2*zo2)
      ro2_xo = -xo2/ro2
      ro2_yo = -yo2/ro2
      ro2_zo = -zo2/ro2
      ro2_x2 =  xo2/ro2
      ro2_y2 =  yo2/ro2
      ro2_z2 =  zo2/ro2
      

*     **** calculate angle ****
      cos_theta = (xo1*xo2 + yo1*yo2 + zo1*zo2)/(ro1*ro2)
      theta = dacos(cos_theta)
      x = 0.5d0*(theta-theta_eq)
      cx    = dcos(x)
      ct    = dcos(theta)
      sx    = dsin(x)
      st    = dsin(theta)
      
      t_ct = -1.0d0/st

      ct_xo = -xo2/(ro1*ro2) - xo1/(ro1*ro2)
     >        - (ct/ro1)*ro1_xo
     >        - (ct/ro2)*ro2_xo

      ct_yo = -yo2/(ro1*ro2) - yo1/(ro1*ro2)
     >        - (ct/ro1)*ro1_yo
     >        - (ct/ro2)*ro2_yo

      ct_zo = -zo2/(ro1*ro2) - zo1/(ro1*ro2)
     >        - (ct/ro1)*ro1_zo
     >        - (ct/ro2)*ro2_zo


      ct_x1 = xo2/(ro1*ro2) - (ct/ro1)*ro1_x1
      ct_y1 = yo2/(ro1*ro2) - (ct/ro1)*ro1_y1
      ct_z1 = zo2/(ro1*ro2) - (ct/ro1)*ro1_z1
 
      ct_x2 = xo1/(ro1*ro2) - (ct/ro2)*ro2_x2
      ct_y2 = yo1/(ro1*ro2) - (ct/ro2)*ro2_y2
      ct_z2 = zo1/(ro1*ro2) - (ct/ro2)*ro2_z2

      s1 = ro1*cx - R_eq
      s2 = ro2*cx - R_eq
      s3 = (ro1+ro2)*sx/R_eq

      s1_ro1 = cx
      s1_t   = -0.5d0*ro1*sx
 
      s2_ro2 = cx
      s2_t   = -0.5d0*ro2*sx
  
      s3_ro1 = sx/R_eq
      s3_ro2 = sx/R_eq

      s3_t   = 0.5d0*(ro1+ro2)*cx/R_eq
       
      E_s1 = 2.0d0*alpha1*D1*dexp(-alpha1*s1)
     >     * (1.0d0 - dexp(-alpha1*s1))

      E_s2 = 2.0d0*alpha2*D2*dexp(-alpha2*s2)
     >     * (1.0d0 - dexp(-alpha2*s2)) 

      E_s3 = 2.0d0*alpha3*D3*dexp(-alpha3*s3)
     >     * (1.0d0 - dexp(-alpha3*s3)) 
     


      E_xo = E_s1*(s1_ro1*ro1_xo                 + s1_t*t_ct*ct_xo)
     >     + E_s2*(s2_ro2*ro2_xo                 + s2_t*t_ct*ct_xo)
     >     + E_s3*(s3_ro1*ro1_xo + s3_ro2*ro2_xo + s3_t*t_ct*ct_xo)

      E_yo = E_s1*(s1_ro1*ro1_yo                 + s1_t*t_ct*ct_yo)
     >     + E_s2*(s2_ro2*ro2_yo                 + s2_t*t_ct*ct_yo)
     >     + E_s3*(s3_ro1*ro1_yo + s3_ro2*ro2_yo + s3_t*t_ct*ct_yo)

      E_zo = E_s1*(s1_ro1*ro1_zo                 + s1_t*t_ct*ct_zo)
     >     + E_s2*(s2_ro2*ro2_zo                 + s2_t*t_ct*ct_zo)
     >     + E_s3*(s3_ro1*ro1_zo + s3_ro2*ro2_zo + s3_t*t_ct*ct_zo)


      E_x1 = E_s1*(s1_ro1*ro1_x1                 + s1_t*t_ct*ct_x1)
     >     + E_s2*(                              + s2_t*t_ct*ct_x1)
     >     + E_s3*(s3_ro1*ro1_x1                 + s3_t*t_ct*ct_x1)

      E_y1 = E_s1*(s1_ro1*ro1_y1                 + s1_t*t_ct*ct_y1)
     >     + E_s2*(                              + s2_t*t_ct*ct_y1)
     >     + E_s3*(s3_ro1*ro1_y1                 + s3_t*t_ct*ct_y1)

      E_z1 = E_s1*(s1_ro1*ro1_z1                 + s1_t*t_ct*ct_z1)
     >     + E_s2*(                              + s2_t*t_ct*ct_z1)
     >     + E_s3*(s3_ro1*ro1_z1                 + s3_t*t_ct*ct_z1)



      E_x2 = E_s1*(                              + s1_t*t_ct*ct_x2)
     >     + E_s2*(s2_ro2*ro2_x2                 + s2_t*t_ct*ct_x2)
     >     + E_s3*(                s3_ro2*ro2_x2 + s3_t*t_ct*ct_x2)

      E_y2 = E_s1*(                              + s1_t*t_ct*ct_y2)
     >     + E_s2*(s2_ro2*ro2_y2                 + s2_t*t_ct*ct_y2)
     >     + E_s3*(                s3_ro2*ro2_y2 + s3_t*t_ct*ct_y2)

      E_z2 = E_s1*(                              + s1_t*t_ct*ct_z2)
     >     + E_s2*(s2_ro2*ro2_z2                 + s2_t*t_ct*ct_z2)
     >     + E_s3*(                s3_ro2*ro2_z2 + s3_t*t_ct*ct_z2)

      fo(1) = -E_xo
      fo(2) = -E_yo
      fo(3) = -E_zo

      f1(1) = -E_x1
      f1(2) = -E_y1
      f1(3) = -E_z1

      f2(1) = -E_x2
      f2(2) = -E_y2
      f2(3) = -E_z2
  


     
      return
      end


c $Id$
