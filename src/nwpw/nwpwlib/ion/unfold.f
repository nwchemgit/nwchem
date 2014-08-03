*
* $Id$
*

      subroutine nwpw_unfold(ni,unita,r1,r2)
*     ======================================================= 
*	This routine unfolds r2
*     ======================================================= 
      implicit none
      integer ni
      real*8  unita(3,3)
      real*8  r1(3,*)
      real*8  r2(3,*)

*     **** Local variables defined ****
      real*8  f1(3),f2(3)
      real*8  a(3,3),b(3,3),volume
      integer i,j


*     ***** Determine the unit lattice vectors and distances ******
      do j=1,3
      do i=1,3
        a(i,j) = unita(i,j)
      end do
      end do

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      volume = a(1,1)*b(1,1)
     >       + a(2,1)*b(2,1)
     >       + a(3,1)*b(3,1)

*     **** no unfolding for small unit cell ****
      if (dabs(volume).lt.2.0) return

      volume = 1.0d0/volume
      call dscal(9,volume,b,1)


      do i =1,ni
	
*        *** Break the Ion positions into the a1, a2, and a3 components ***
         f1(1) =  b(1,1) * r1(1,i)
     >         +  b(2,1) * r1(2,i)
     >         +  b(3,1) * r1(3,i)

         f1(2) =  b(1,2) * r1(1,i)
     >         +  b(2,2) * r1(2,i)
     >         +  b(3,2) * r1(3,i)

         f1(3) =  b(1,3) * r1(1,i)
     >         +  b(2,3) * r1(2,i)
     >         +  b(3,3) * r1(3,i)

         f2(1) =  b(1,1) * r2(1,i)
     >         +  b(2,1) * r2(2,i)
     >         +  b(3,1) * r2(3,i)

         f2(2) =  b(1,2) * r2(1,i)
     >         +  b(2,2) * r2(2,i)
     >         +  b(3,2) * r2(3,i)

         f2(3) =  b(1,3) * r2(1,i)
     >         +  b(2,3) * r2(2,i)
     >         +  b(3,3) * r2(3,i)


*	**** Change the a1, a2 and a3  components to ****

   23   IF ((f2(1)-f1(1)) .GT. (0.5d0)) THEN
*          WRITE (*,*) 'a1>', I, R2A1, DA(1)/2.0d0
           r2(1,i) = r2(1,i) - unita(1,1)
           r2(2,i) = r2(2,i) - unita(2,1)
           r2(3,i) = r2(3,i) - unita(3,1)

            f2(1) =  b(1,1) * r2(1,i)
     >            +  b(2,1) * r2(2,i)
     >            +  b(3,1) * r2(3,i)
           GO TO 23
        ENDIF
	   
   24   IF ((f2(1)-f1(1)) .LE. (-0.5d0)) THEN
*          WRITE (*,*) 'a1<', I, R2A1, DA(1)/2.0d0
           r2(1,i) = r2(1,i) + unita(1,1)
           r2(2,i) = r2(2,i) + unita(2,1)
           r2(3,i) = r2(3,i) + unita(3,1)

            f2(1) =  b(1,1) * r2(1,i)
     >            +  b(2,1) * r2(2,i)
     >            +  b(3,1) * r2(3,i)
            GO TO 24
        ENDIF

   25   IF ((f2(2)-f1(2)) .GT. (0.5d0)) THEN
*          WRITE (*,*) 'a2>', I, R2A2, DA(2)/2.0d0
           r2(1,i) = r2(1,i) - unita(1,2)
           r2(2,i) = r2(2,i) - unita(2,2)
           r2(3,i) = r2(3,i) - unita(3,2)

            f2(2) =  b(1,2) * r2(1,i)
     >            +  b(2,2) * r2(2,i)
     >            +  b(3,2) * r2(3,i)
          GO TO 25
        ENDIF
	   
   26   IF ((f2(2)-f1(2)) .LE. (-0.5d0)) THEN
*          WRITE (*,*) 'a2<', I, R2A2, DA(2)/2.0d0
           r2(1,i) = r2(1,i) + unita(1,2)
           r2(2,i) = r2(2,i) + unita(2,2)
           r2(3,i) = r2(3,i) + unita(3,2)

            f2(2) =  b(1,2) * r2(1,i)
     >            +  b(2,2) * r2(2,i)
     >            +  b(3,2) * r2(3,i)
           GO TO 26
        ENDIF


   27   IF ((f2(3)-f1(3)) .GT. (0.5d0)) THEN
*         WRITE (*,*) 'a3>', i, R2A3, DA(3)/2.0d0
          r2(1,i) = r2(1,i) - unita(1,3)
          r2(2,i) = r2(2,i) - unita(2,3)
          r2(3,i) = r2(3,i) - unita(3,3)

          f2(3) =  b(1,3) * r2(1,i)
     >          +  b(2,3) * r2(2,i)
     >          +  b(3,3) * r2(3,i)
           GO TO 27
        ENDIF
	   
   28   IF ((f2(3)-f1(3)) .LE. (-0.5d0)) THEN
*         WRITE (*,*) 'a3<', I, R2A3, DA(3)/2.0d0
          r2(1,i) = r2(1,i) + unita(1,3)
          r2(2,i) = r2(2,i) + unita(2,3)
          r2(3,i) = r2(3,i) + unita(3,3)

          f2(3) =  b(1,3) * r2(1,i)
     >          +  b(2,3) * r2(2,i)
     >          +  b(3,3) * r2(3,i)
           GO TO 28
        ENDIF

      end do 

      return
      end
