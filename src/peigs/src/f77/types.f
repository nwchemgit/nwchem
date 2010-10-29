*
* $Id$
*
      subroutine type1(a,b,alpha,beta,q,e,NIN,NDIM,record)
C
C-----------------------------------------------------------------------
C
C     This routine reads Fann's matrices, i.e., alpha's/beta's,
C     and computes the representation a's/b's, and q's/e's. 
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
      character        record*48
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          i,j
      double precision shift
      character        matrix*24,string*24
C
      intrinsic        sqrt
C
C.... read the alpha's and beta's ......................................
C
cgif      call parser (record,string)
cgif      read (string,fmt='(a)') matrix
cgif      call parser (record,string)
cgif      read (string,fmt='(e24.0)') shift
C
cgif      open (9,err=1,file=matrix,form='formatted')
cgif      read (9,err=1,fmt=*) NIN 
cgif      if ( NIN.GT.NDIM ) stop 'TYPE1: not enough space in Z'
cgif      read (9,err=1,fmt=*) (j,alpha(i),beta(i),i=1,NIN)
C
cgif      do j = 1,NIN-1
cgif         alpha(j) = alpha(j) + shift
cgif         beta(j) = beta(j+1)
cgif      end do
cgif      alpha(NIN) = alpha(NIN) + shift
C
C.... compute the q's and e's ..........................................
C
      q(1) = alpha(1)
      do j = 1,NIN-1
         e(j) = (beta(j)/q(j)) * beta(j)
         q(j+1) = alpha(j+1) - e(j) 
         if ( q(j+1).lt.zero ) 
     &      STOP 'Matrix is not positive definite, change shift'
      end do 
C
C.... compute the a's and b's ..........................................
C
      a(1) = sqrt(q(1))
      do i=1,NIN-1
         b(i) = sqrt(e(i))
         a(i+1) = sqrt(q(i+1))
      end do
C
      return
    1 stop '* I/O error, check file for Fann''s matrix'
      end
C
C***********************************************************************
C
      subroutine TYPE2 (a,b,alpha,beta,q,e,NIN,NDIM,record)
C
C-----------------------------------------------------------------------
C
C     This routine reads parameters for Gamma glued matrices, i.e.,
C     ncopies and gamma, and computes the representation alpha's/
C     beta's, a's/b's, and q's/e's.      
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
      character        record*48
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          j,k,ncopies,N
      double precision gamma
      character        string*24
C
      intrinsic        dble
C
      call parser (record,string)
      read (string,fmt='(i24)') ncopies
      call parser (record,string)
      read (string,fmt='(e24.0)') gamma
C
      N = 11
      NIN = N*ncopies 
      if ( NIN.GT.NDIM ) stop 'TYPE2: not enough space in Z'
C
C.... compute the a's and b's ..........................................
C
      do j = 1,NIN-1
         b(j) = dble(1)
      end do
      do k = 0,ncopies-1
         a( 1+k*11) = dble( 1)
         a( 2+k*11) = dble(11)
         a( 3+k*11) = dble(21)
         a( 4+k*11) = dble(31)
         a( 5+k*11) = dble(41)
         a( 6+k*11) = dble(51)
         a( 7+k*11) = dble(41)
         a( 8+k*11) = dble(31)
         a( 9+k*11) = dble(21)
         a(10+k*11) = dble(11)
         a(11+k*11) = dble( 1)
         b(10+k*11) = gamma
      end do
C
C.... compute the alpha's and beta's ...................................
C
      alpha(1) = a(1)**2
      do j = 1,NIN-1
         beta(j) = a(j)*b(j)
         alpha(j+1) = a(j+1)**2 + b(j)**2
      end do
C
C.... compute the q's and e's ..........................................
C
      do j = 1,NIN-1
         q(j) = a(j)**2
         e(j) = b(j)**2 
      end do
      q(NIN) = a(NIN)**2
      e(NIN) = zero
C
      return
      end
C
C***********************************************************************
C
      subroutine TYPE3 (a,b,alpha,beta,q,e,NIN,NDIM,record)
C
C-----------------------------------------------------------------------
C
C     This routine reads Inder's matrices, i.e., a's and b's, and
C     computes the representation alpha's/beta's and q's/e's. 
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
      character        record*48
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          i
      character        matrix*24,string*24
C
C.... read the a's and b's .............................................
C
      call parser (record,string)
      read (string,fmt='(a)') matrix
C
      open (9,err=1,file=matrix,form='formatted')
      read (9,err=1,fmt=*) NIN 
      if ( NIN.GT.NDIM ) stop 'TYPE3: not enough space in Z'
C
      read (9,err=1,fmt=*) (a(i),i=1,NIN  )
      read (9,err=1,fmt=*) (b(i),i=1,NIN-1)
C
C.... compute the alpha's and beta's ...................................
C
      alpha(1) = a(1)**2
      do i = 1,NIN-1
         beta(i) = a(i)*b(i)
         alpha(i+1) = a(i+1)**2 + b(i)**2
      end do
C
C.... compute the q's and e's ..........................................
C
      do i = 1,NIN-1
         q(i) = a(i)**2
         e(i) = b(i)**2 
      end do
      q(NIN) = a(NIN)**2
      e(NIN) = zero
C
      return
    1 stop '* I/O error, check file for Inder''s matrix'
      end
C
C***********************************************************************
C
      subroutine TYPE4 (a,b,alpha,beta,q,e,NIN,NDIM,record)
C
C-----------------------------------------------------------------------
C
C     This routine generates a one_two_one matrix of dimension NIN
C     in the representation alpha's/beta's, a's/b's, and q's/e's.
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
      character        record*48
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          i,j
      character        string*24
C
      intrinsic        sqrt
C
C.... read the alpha's and beta's ......................................
C
      call parser (record,string)
      read (string,fmt='(i24)') NIN

      if ( NIN.GT.NDIM ) stop 'TYPE4: not enough space in Z'
C
      do j = 1,NIN-1
         alpha(j) = 2.0d0
         beta(j) = 1.0d0
      end do
      alpha(NIN) = 2.0d0
C
C.... compute the q's and e's ..........................................
C
      q(1) = alpha(1)
      do j = 1,NIN-1
         e(j) = (beta(j)/q(j)) * beta(j)
         q(j+1) = alpha(j+1) - e(j) 
         if ( q(j+1).lt.zero ) 
     &      STOP 'Matrix is not positive definite, change shift'
      end do 
C
C.... compute the a's and b's ..........................................
C
      a(1) = sqrt(q(1))
      do i=1,NIN-1
         b(i) = sqrt(e(i))
         a(i+1) = sqrt(q(i+1))
      end do
C
      return
      end
C
C***********************************************************************
C
      subroutine TYPE5 (a,b,alpha,beta,q,e,NIN,NDIM,record)
C
C-----------------------------------------------------------------------
C
C     This routine reads parameters for generating a random matrix
C     and computes the representation alpha's/beta's, a's/b's, and 
C     q's/e's.      
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
      character        record*48
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          iseed(4),j
      character        string*24
C
      double precision dnrm2
      intrinsic        dble
C
      call parser (record,string)
      read (string,fmt='(i24)') NIN
C
      if ( NIN.GT.NDIM ) stop 'TYPE5: not enough space in Z'
C
C.... compute the a's and b's ..........................................
C
      do j = 1,NIN
         iseed(1) = j
         iseed(2) = j + 1
         iseed(3) = j + 2
         iseed(4) = 2*j + 3
         call dlarnv(3,iseed,NIN+1-j,alpha)
         a(j) = dnrm2(NIN+1-j,alpha,1)
         iseed(1) = j
         iseed(2) = j + 1
         iseed(3) = j + 2
         iseed(4) = 2*j + 11
         call dlarnv(3,iseed,NIN-j,alpha)
         b(j) = dnrm2(NIN-j,alpha,1)
      end do
C
C.... compute the alpha's and beta's ...................................
C
      alpha(1) = a(1)**2
      do j = 1,NIN-1
         beta(j) = a(j)*b(j)
         alpha(j+1) = a(j+1)**2 + b(j)**2
      end do
C
C.... compute the q's and e's ..........................................
C
      do j = 1,NIN-1
         q(j) = a(j)**2
         e(j) = b(j)**2 
      end do
      q(NIN) = a(NIN)**2
      e(NIN) = zero
C
      return
      end
C
C***********************************************************************
C
      subroutine parser (record,string)
C
      integer   i,j
      character char*1,record*40,string*24
      logical   advance
C
      i = 0
      j = 0
      string = ' '
      advance = .true.
C
      do while ( advance )
         i = i + 1
         char = record(i:i)
         if ( char.ne.' ' .and. char.ne.',' ) then
            j = j + 1
            string(j:j) = char
         else 
            record(i:i) = ' '
            advance = j.eq.0 
         end if
      end do
C
      record = record(i+1:40)
C
      return
      end
