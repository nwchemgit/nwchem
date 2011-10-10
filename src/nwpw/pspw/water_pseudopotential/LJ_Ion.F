*     **********************************
*     *	                               *
*     *          LJ_Ion_Repulsion      *
*     *                                *
*     **********************************

*   This routine calculates the LJ Repulsion between two centers 
*
      subroutine LJ_Ion_Repulsion(r1,e1,s1,r2,e2,s2,er)
      implicit none
      real*8 r1(3),e1,s1
      real*8 r2(3),e2,s2
      real*8 er

*     **** local variables ****
      real*8 xo,yo,zo,d2
      real*8 e,s
     
      e = dsqrt(e1*e2)
      s = dsqrt(s1*s2)

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)
        
      er = 4.0d0*e * ((s*s/d2)**6 - (s*s/d2)**3)    
      return 
      end

*     **********************************
*     *	                               *
*     *      LJ_Ion_Repulsion_Force    *
*     *                                *
*     **********************************

*   This routine calculates the LJ Repulsion force between two centers
*
      subroutine LJ_Ion_Repulsion_Force(r1,e1,s1,r2,e2,s2,f2)
      implicit none
      real*8 r1(3),e1,s1
      real*8 r2(3),e2,s2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,do,d2
      real*8 der
      real*8 e,s
     
      e = dsqrt(e1*e2)
      s = dsqrt(s1*s2)

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)
      do = dsqrt(d2)

      der = (4.0d0*e/do) * ( -12.0d0*(s*s/d2)**6 
     >                       +  6.0d0*(s*s/d2)**3 ) 
         
      f2(1) = -(xo/do)*der
      f2(2) = -(yo/do)*der
      f2(3) = -(zo/do)*der  

      return 
      end

c $Id$
