*
* $Id$
*
      subroutine incell2(ni,r2,r1)

*     ======================================================= 
*	This routine places the ions inside the cell defined 
*	by the lattice vectors centered at zero.  
*     ======================================================= 
      implicit none
      integer ni
      real*8  r2(3,*),r1(3,*)

*     **** Local variables defined ****
      real*8  fa1,fa2,fa3
      real*8  a(3,3),b(3,3),volume
      integer i,j

*      **** external functions ****
       real*8   lattice_unita
       external lattice_unita

*     ***** Determine the unit lattice vectors and distances ******
      do j=1,3
      do i=1,3
        a(i,j) = lattice_unita(i,j)
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

      volume = 1.0d0/volume
      call dscal(9,volume,b,1)

!$OMP DO
      do i =1,ni

*        *** Break the Ion positions into the a1, a2, and a3 components ***
         fa1 =  b(1,1) * r2(1,i)
     >       +  b(2,1) * r2(2,i)
     >       +  b(3,1) * r2(3,i)

         fa2 =  b(1,2) * r2(1,i)
     >       +  b(2,2) * r2(2,i)
     >       +  b(3,2) * r2(3,i)

         fa3 =  b(1,3) * r2(1,i)
     >       +  b(2,3) * r2(2,i)
     >       +  b(3,3) * r2(3,i)


*	**** Change the a1, a2 and a3  components to ****
*	**** make the ion be in the cell             ****

   23   IF (fa1 .GT. (0.5d0)) THEN
*          WRITE (*,*) 'a1>', I, R2A1, DA(1)/2.0d0
           r2(1,i) = r2(1,i) - lattice_unita(1,1)
           r2(2,i) = r2(2,i) - lattice_unita(2,1)
           r2(3,i) = r2(3,i) - lattice_unita(3,1)

           r1(1,i) = r1(1,i) - lattice_unita(1,1)
           r1(2,i) = r1(2,i) - lattice_unita(2,1)
           r1(3,i) = r1(3,i) - lattice_unita(3,1)


           fa1 =  b(1,1) * r2(1,i)
     >         +  b(2,1) * r2(2,i)
     >         +  b(3,1) * r2(3,i)
           GO TO 23
        ENDIF
   
   24   IF (fa1 .LE. (-0.5d0)) THEN
*          WRITE (*,*) 'a1<', I, R2A1, DA(1)/2.0d0
           r2(1,i) = r2(1,i) + lattice_unita(1,1)
           r2(2,i) = r2(2,i) + lattice_unita(2,1)
           r2(3,i) = r2(3,i) + lattice_unita(3,1)

           r1(1,i) = r1(1,i) + lattice_unita(1,1)
           r1(2,i) = r1(2,i) + lattice_unita(2,1)
           r1(3,i) = r1(3,i) + lattice_unita(3,1)


           fa1 =  b(1,1) * r2(1,i)
     >         +  b(2,1) * r2(2,i)
     >         +  b(3,1) * r2(3,i)
            GO TO 24
        ENDIF

   25   IF (fa2 .GT. (0.5d0)) THEN
*          WRITE (*,*) 'a2>', I, R2A2, DA(2)/2.0d0
           r2(1,i) = r2(1,i) - lattice_unita(1,2)
           r2(2,i) = r2(2,i) - lattice_unita(2,2)
           r2(3,i) = r2(3,i) - lattice_unita(3,2)

           r1(1,i) = r1(1,i) - lattice_unita(1,2)
           r1(2,i) = r1(2,i) - lattice_unita(2,2)
           r1(3,i) = r1(3,i) - lattice_unita(3,2)


           fa2 =  b(1,2) * r2(1,i)
     >         +  b(2,2) * r2(2,i)
     >         +  b(3,2) * r2(3,i)
          GO TO 25
        ENDIF
   
   26   IF (fa2 .LE. (-0.5d0)) THEN
*          WRITE (*,*) 'a2<', I, R2A2, DA(2)/2.0d0
           r2(1,i) = r2(1,i) + lattice_unita(1,2)
           r2(2,i) = r2(2,i) + lattice_unita(2,2)
           r2(3,i) = r2(3,i) + lattice_unita(3,2)

           r1(1,i) = r1(1,i) + lattice_unita(1,2)
           r1(2,i) = r1(2,i) + lattice_unita(2,2)
           r1(3,i) = r1(3,i) + lattice_unita(3,2)


           fa2 =  b(1,2) * r2(1,i)
     >         +  b(2,2) * r2(2,i)
     >         +  b(3,2) * r2(3,i)
           GO TO 26
        ENDIF


   27   IF (fa3 .GT. (0.5d0)) THEN
*         WRITE (*,*) 'a3>', i, R2A3, DA(3)/2.0d0
          r2(1,i) = r2(1,i) - lattice_unita(1,3)
          r2(2,i) = r2(2,i) - lattice_unita(2,3)
          r2(3,i) = r2(3,i) - lattice_unita(3,3)

          r1(1,i) = r1(1,i) - lattice_unita(1,3)
          r1(2,i) = r1(2,i) - lattice_unita(2,3)
          r1(3,i) = r1(3,i) - lattice_unita(3,3)

          fa3 =  b(1,3) * r2(1,i)
     >        +  b(2,3) * r2(2,i)
     >        +  b(3,3) * r2(3,i)
           GO TO 27
        ENDIF
  
   28   IF (fa3 .LE. (-0.5d0)) THEN
*         WRITE (*,*) 'a3<', I, R2A3, DA(3)/2.0d0
          r2(1,i) = r2(1,i) + lattice_unita(1,3)
          r2(2,i) = r2(2,i) + lattice_unita(2,3)
          r2(3,i) = r2(3,i) + lattice_unita(3,3)

          r1(1,i) = r1(1,i) + lattice_unita(1,3)
          r1(2,i) = r1(2,i) + lattice_unita(2,3)
          r1(3,i) = r1(3,i) + lattice_unita(3,3)

          fa3 =  b(1,3) * r2(1,i)
     >        +  b(2,3) * r2(2,i)
     >        +  b(3,3) * r2(3,i)
           GO TO 28
        ENDIF

      end do 
!$OMP END DO

      return
      end
