c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.1  2003/02/19 23:41:16  marat
c added double factorial, hypergeometric, and gamma function
c
c MV
c
c

!**************************************************
!
!   Name: paw_double_factorial
!
!   Purpose: calculates double factorial n!!
!                for n from -1 to 16      
!
!  Created: Marat 8/05/2002
!**************************************************
      function paw_double_factorial(n)
      
      implicit none
      integer n
      integer paw_double_factorial
      integer ndf(18)


      data ndf /    1,  
     >               1,
     >               1,
     >               2,
     >               3,
     >               8,
     >               15,
     >               48,
     >               105,
     >               384,
     >               945,
     >               3840,
     >               10395,
     >               46080,
     >               135135,
     >               645120,
     >               2027025,
     >               10321920/


      
      paw_double_factorial = ndf(n+2)
      return
      end
