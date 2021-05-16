#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
      subroutine dradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
*      
* $Id$
*      
      double precision cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(1), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c      
      integer j1, j3, j4, j5, j6, j7, ik1, ik2, ik3, i1, i2
      doubleprecision d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12,
     . d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, 
     .d26, d27, d28, d29, d30
      arg = tpi/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
      i2 = iand(max0(ido,0),3)
      do i = 1, i2
         ch(i,k,1) = cc(i,1,k)
      end do
      do i = i2 + 1, ido, 4
         d27 = cc(i,1,k)
         d28 = cc(i+1,1,k)
         d29 = cc(i+2,1,k)
         d30 = cc(i+3,1,k)
         ch(i,k,1) = d27
         ch(i+1,k,1) = d28
         ch(i+2,k,1) = d29
         ch(i+3,k,1) = d30
      end do
  102 continue
      go to 106
c      
  103 continue
      do k = 1, l1
         i1 = iand(max0(ido,0),3)
         do i = 1, i1
            ch(i,k,1) = cc(i,1,k)
         end do
         do i = i1 + 1, ido, 4
            d23 = cc(i,1,k)
            d24 = cc(i+1,1,k)
            d25 = cc(i+2,1,k)
            d26 = cc(i+3,1,k)
            ch(i,k,1) = d23
            ch(i+1,k,1) = d24
            ch(i+2,k,1) = d25
            ch(i+3,k,1) = d26
         end do
      end do
c      
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
      do k = 1, l1
         ch(1,k,j) = cc(ido,j2-2,k) + cc(ido,j2-2,k)
         ch(1,k,jc) = cc(1,j2-1,k) + cc(1,j2-1,k)
      end do
  108 continue
c      
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
      do i = 1, ((ido - 1)/2)
         ch(i*2,k,j) = cc(i*2,2*j-1,k) + cc(idp2-(i+1)*2,2*j-2,k)
         ch(i*2,k,jc) = cc(i*2,2*j-1,k) - cc(idp2-(i+1)*2,2*j-2,k)
         ch(1+i*2,k,j) = cc(1+i*2,2*j-1,k) - cc(idp2-1-i*2,2*j-2,k)
         ch(1+i*2,k,jc) = cc(1+i*2,2*j-1,k) + cc(idp2-1-i*2,2*j-2,k)
      end do
  110    continue
  111 continue
      go to 116
c      
  112 do 115 j=2,ipph
         jc = ipp2-j
      do k = 1, l1
         do i = 1, ((ido - 1)/2)
            ch(i*2,k,j) = cc(i*2,2*j-1,k) + cc(idp2-(i+1)*2,2*j-2,k)
            ch(i*2,k,jc) = cc(i*2,2*j-1,k) - cc(idp2-(i+1)*2,2*j-2,k)
            ch(1+i*2,k,j) = cc(1+i*2,2*j-1,k) - cc(idp2-1-i*2,2*j-2,k)
            ch(1+i*2,k,jc) = cc(1+i*2,2*j-1,k) + cc(idp2-1-i*2,2*j-2,k)
         end do
      end do
  115 continue
c      
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
      do ik = 1, idl1
         c2(ik,l) = ch2(ik,1) + ar1*ch2(ik,2)
         c2(ik,lc) = ai1*ch2(ik,ip)
      end do
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
      do ik = 1, idl1
         c2(ik,l) = c2(ik,l) + ar2*ch2(ik,j)
         c2(ik,lc) = c2(ik,lc) + ai2*ch2(ik,jc)
      end do
  119    continue
  120 continue
c      
      if (ipph - 1 .gt. 0) then
         j1 = iand(ipph - 1,3)
c
         do j = 1, j1
            ik3 = iand(max0(idl1,0),3)
            do ik = 1, ik3
               ch2(ik,1) = ch2(ik,1) + ch2(ik,j+1)
            end do
            do ik = ik3 + 1, idl1, 4
               d19 = ch2(ik,1) + ch2(ik,j+1)
               d20 = ch2(ik+1,1) + ch2(ik+1,j+1)
               d21 = ch2(ik+2,1) + ch2(ik+2,j+1)
               d22 = ch2(ik+3,1) + ch2(ik+3,j+1)
               ch2(ik,1) = d19
               ch2(ik+1,1) = d20
               ch2(ik+2,1) = d21
               ch2(ik+3,1) = d22
            end do
         end do
c
         do j = j1 + 1, ipph - 1, 4
            ik2 = iand(max0(idl1,0),3)
            do ik = 1, ik2
               ch2(ik,1) = ch2(ik,1) + ch2(ik,j+1) + ch2(ik,j+2) + ch2(
     1            ik,j+3) + ch2(ik,j+4)
            end do
            do ik = ik2 + 1, idl1, 4
               d15 = ch2(ik,1) + ch2(ik,j+1) + ch2(ik,j+2) + ch2(ik,j+3)
     1             + ch2(ik,j+4)
               d16 = ch2(ik+1,1) + ch2(ik+1,j+1) + ch2(ik+1,j+2) + ch2(
     1            ik+1,j+3) + ch2(ik+1,j+4)
               d17 = ch2(ik+2,1) + ch2(ik+2,j+1) + ch2(ik+2,j+2) + ch2(
     1            ik+2,j+3) + ch2(ik+2,j+4)
               d18 = ch2(ik+3,1) + ch2(ik+3,j+1) + ch2(ik+3,j+2) + ch2(
     1            ik+3,j+3) + ch2(ik+3,j+4)
               ch2(ik,1) = d15
               ch2(ik+1,1) = d16
               ch2(ik+2,1) = d17
               ch2(ik+3,1) = d18
            end do
         end do
      endif
c      
      do 124 j=2,ipph
         jc = ipp2-j
      do k = 1, l1
         ch(1,k,j) = c1(1,k,j) - c1(1,k,jc)
         ch(1,k,jc) = c1(1,k,j) + c1(1,k,jc)
      end do
  124 continue
c      
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
      do i = 1, ((ido - 1)/2)
         ch(i*2,k,j) = c1(i*2,k,j) - c1(1+i*2,k,jc)
         ch(i*2,k,jc) = c1(i*2,k,j) + c1(1+i*2,k,jc)
         ch(1+i*2,k,j) = c1(1+i*2,k,j) + c1(i*2,k,jc)
         ch(1+i*2,k,jc) = c1(1+i*2,k,j) - c1(i*2,k,jc)
      end do
  126    continue
  127 continue
      go to 132
c      
  128 do 131 j=2,ipph
         jc = ipp2-j
      do k = 1, l1
         do i = 1, ((ido - 1)/2)
            ch(i*2,k,j) = c1(i*2,k,j) - c1(1+i*2,k,jc)
            ch(i*2,k,jc) = c1(i*2,k,j) + c1(1+i*2,k,jc)
            ch(1+i*2,k,j) = c1(1+i*2,k,j) + c1(i*2,k,jc)
            ch(1+i*2,k,jc) = c1(1+i*2,k,j) - c1(i*2,k,jc)
         end do
      end do
  131 continue
  132 continue
c      
      if (ido .eq. 1) return
      ik1 = iand(max0(idl1,0),3)
      do ik = 1, ik1
         c2(ik,1) = ch2(ik,1)
      end do
      do 133 ik = ik1 + 1, idl1, 4
         d11 = ch2(ik,1)
         d12 = ch2(ik+1,1)
         d13 = ch2(ik+2,1)
         d14 = ch2(ik+3,1)
         c2(ik,1) = d11
         c2(ik+1,1) = d12
         c2(ik+2,1) = d13
         c2(ik+3,1) = d14
  133 continue
c
      j7 = iand(max0((ip - 1)*l1,0),3)
      do j = 1, j7
         c1(1,j,2) = ch(1,j,2)
      end do
      do j = j7 + 1, (ip - 1)*l1, 4
         d7 = ch(1,j,2)
         d8 = ch(1,j+1,2)
         d9 = ch(1,j+2,2)
         d10 = ch(1,j+3,2)
         c1(1,j,2) = d7
         c1(1,j+1,2) = d8
         c1(1,j+2,2) = d9
         c1(1,j+3,2) = d10
      end do
c      
      if (nbd .gt. l1) go to 139
      is = -ido
      do i = 1, ((ido - 1)/2)
         if (l1 .gt. 0) then
            j3 = iand(l1,1)
            do k = 1, j3
               j4 = 0
               do j = 1, ip - 1
                  c1(2*i,k,1+j) = wa(ido+2*i-1+is+j4)*ch(2*i,k,1+j) - wa
     1               (ido+2*i+is+j4)*ch(2*i+1,k,1+j)
                  c1(2*i+1,k,1+j) = wa(ido+2*i-1+is+j4)*ch(2*i+1,k,1+j)
     1                + wa(ido+2*i+is+j4)*ch(2*i,k,1+j)
                  j4 = j4 + ido
               end do
            end do
            do k = j3 + 1, l1, 2
               j5 = 0
               do j = 1, ip - 1
                  d2 = wa(ido+2*i+is+j5)
                  d1 = wa(ido+2*i-1+is+j5)
                  d3 = d1*ch(2*i,k,1+j) - d2*ch(2*i+1,k,1+j)
                  d4 = d1*ch(2*i,k+1,1+j) - d2*ch(2*i+1,k+1,1+j)
                  c1(2*i,k,1+j) = d3
                  c1(2*i,k+1,1+j) = d4
                  d5 = d1*ch(2*i+1,k,1+j) + d2*ch(2*i,k,1+j)
                  d6 = d1*ch(2*i+1,k+1,1+j) + d2*ch(2*i,k+1,1+j)
                  c1(2*i+1,k,1+j) = d5
                  c1(2*i+1,k+1,1+j) = d6
                  j5 = j5 + ido
               end do
            end do
         endif
      end do
      go to 143
c      
  139 is = -ido
      do k = 1, l1
         do i = 1, ((ido - 1)/2)
            j6 = 0
            do j = 1, ip - 1
               c1(2*i,k,1+j) = wa(ido+2*i-1+is+j6)*ch(2*i,k,1+j) - wa(
     1            ido+2*i+is+j6)*ch(2*i+1,k,1+j)
               c1(2*i+1,k,1+j) = wa(ido+2*i-1+is+j6)*ch(2*i+1,k,1+j) + 
     1            wa(ido+2*i+is+j6)*ch(2*i,k,1+j)
               j6 = j6 + ido
            end do
         end do
      end do
c      
  143 return
      end
