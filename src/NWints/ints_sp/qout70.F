      subroutine qout70(gout, eri, isht, jsht, ksht, lsht,
     $     ib1, jb1, kb1, lb1)
c $Id$
      implicit none
      double precision gout(*)  ! [input]
      double precision eri(*)   ! [output]
      integer isht,jsht,ksht,lsht ! [input]
      integer ib1, jb1, kb1, lb1 ! [input]
c
c     Take the integrals from the internal buffer of the g70
c     integral routines and, permuting as necessary, shove
c     them into the final output buffer, taking into account
c     that if a shell was a p, rather than an sp, some integrals
c     need to be discarded.
c
c     !! NOTE the meaning of ib is different to that of GAMESS
c     
      integer i, j, k, l
      integer n1, n2
      integer kmin(-1:1), kmax(-1:1) ! Index range for shell types
      integer ib(4,4)           ! Offsets in gout for permuted indices
      integer mini, maxi, minj, maxj, mink, maxk, minl, maxl, ijkl
      save ib,kmin,kmax         ! So that init. is done once only
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data kmin/1,1,2/          ! sp = 1,4; s = 1,1;  p = 2,4;
      data kmax/4,1,4/
c     
      mini = kmin(isht)
      minj = kmin(jsht)
      mink = kmin(ksht)
      minl = kmin(lsht)
      maxi = kmax(isht)
      maxj = kmax(jsht)
      maxk = kmax(ksht)
      maxl = kmax(lsht)
      ijkl = 0
      do i = mini,maxi
         do j = minj,maxj
            n1 = ib(ib1,i) + ib(jb1,j) + 1
            do k = mink,maxk
               n2 = n1 + ib(kb1,k)
               do l = minl,maxl
                  ijkl = ijkl + 1
                  eri(ijkl) = gout(n2 + ib(lb1,l))
               enddo
            enddo
         enddo
      enddo
      end
