#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:36:36
C****************************************************************
C
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2
C   Version 6.1C1 on 12/ 8/97 at 23:35:25
C****************************************************************
C      
      subroutine dradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(1), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c      
      integer j1, j3, j4, j5, j6, k1, k2, j7, j8, i1, i2, ik1, ik2, ik3,
     . j9, ik4
      doubleprecision 
     . d16, d17, d18, d20, d21, d22
      integer j10, j11, k3, k4, j12, j13, j14, j15, j16, i3, i4, i5, i6,
     . ik5, ik6, ik7, ik8, k5, j17, j18, ik9, ik10
      doubleprecision d35, d36, d37, d38, d39, d40, d41, d42, d43, d44, 
     .d45, d46, d47, d48, d49, d50, d51, d52, d53, d54, d55, d56, d57, d
     .58, d59, d60, d61, d62, d63, d64, d65, d66, d67, d68, d69, d70, d7
     .1, d72, d73, d74, d75, d76, d77, d78, d79, d80, d81, d82, d83, d84
     ., d85, d86, d87, d88
      arg = tpi/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      ik4 = iand(max0(idl1,0),3)
      ik10 = iand(max0(ik4,0),3)
      do ik = 1, ik10
         ch2(ik,1) = c2(ik,1)
      end do
      do ik = ik10 + 1, ik4, 4
         d85 = c2(ik,1)
         d86 = c2(ik+1,1)
         d87 = c2(ik+2,1)
         d88 = c2(ik+3,1)
         ch2(ik,1) = d85
         ch2(ik+1,1) = d86
         ch2(ik+2,1) = d87
         ch2(ik+3,1) = d88
      end do
      ik9 = iand(max0(((idl1 + 3 - ik4)/4)*4,0),3)
      do ik = 1, ik9
         ch2(ik4+ik,1) = c2(ik4+ik,1)
      end do
      do ik = ik9 + 1, ((idl1 + 3 - ik4)/4)*4, 4
         d81 = c2(ik4+ik,1)
         d82 = c2(ik4+ik+1,1)
         d83 = c2(ik4+ik+2,1)
         d84 = c2(ik4+ik+3,1)
         ch2(ik4+ik,1) = d81
         ch2(ik4+ik+1,1) = d82
         ch2(ik4+ik+2,1) = d83
         ch2(ik4+ik+3,1) = d84
      end do
      j9 = iand(max0((ip - 1)*l1,0),3)
      j18 = iand(max0(j9,0),3)
      do j = 1, j18
         ch(1,j,2) = c1(1,j,2)
      end do
      do j = j18 + 1, j9, 4
         d77 = c1(1,j,2)
         d78 = c1(1,j+1,2)
         d79 = c1(1,j+2,2)
         d80 = c1(1,j+3,2)
         ch(1,j,2) = d77
         ch(1,j+1,2) = d78
         ch(1,j+2,2) = d79
         ch(1,j+3,2) = d80
      end do
      j17 = iand(max0((((ip - 1)*l1 + 3 - j9)/4)*4,0),3)
      do j = 1, j17
         ch(1,j9+j,2) = c1(1,j9+j,2)
      end do
      do j = j17 + 1, (((ip - 1)*l1 + 3 - j9)/4)*4, 4
         d73 = c1(1,j9+j,2)
         d74 = c1(1,j9+j+1,2)
         d75 = c1(1,j9+j+2,2)
         d76 = c1(1,j9+j+3,2)
         ch(1,j9+j,2) = d73
         ch(1,j9+j+1,2) = d74
         ch(1,j9+j+2,2) = d75
         ch(1,j9+j+3,2) = d76
      end do
c      
      if (nbd .gt. l1) go to 107
      is = -ido
      do i = 1, (ido - 1)/2
         if (l1 .gt. 0) then
            j1 = iand(l1,1)
      if (j1 .gt. 0) then
         j12 = iand(j1,1)
         do k = 1, j12
            j3 = 0
            j13 = 0
            do j = 1, ip - 1
               ch(2*i,k,1+j) = wa(ido+2*i-1+is+j3+j13)*c1(2*i,k,1+j) + 
     1            wa(ido+2*i+is+j3+j13)*c1(2*i+1,k,1+j)
               ch(2*i+1,k,1+j) = wa(ido+2*i-1+is+j3+j13)*c1(2*i+1,k,1+j)
     1             - wa(ido+2*i+is+j3+j13)*c1(2*i,k,1+j)
               j13 = j13 + ido
            end do
         end do
         do k = j12 + 1, j1, 2
            j3 = 0
            j14 = 0
            do j = 1, ip - 1
               d36 = wa(ido+2*i+is+j3+j14)
               d35 = wa(ido+2*i-1+is+j3+j14)
               d37 = d35*c1(2*i,k,1+j) + d36*c1(2*i+1,k,1+j)
               d38 = d35*c1(2*i,k+1,1+j) + d36*c1(2*i+1,k+1,1+j)
               ch(2*i,k,1+j) = d37
               ch(2*i,k+1,1+j) = d38
               d39 = d35*c1(2*i+1,k,1+j) - d36*c1(2*i,k,1+j)
               d40 = d35*c1(2*i+1,k+1,1+j) - d36*c1(2*i,k+1,1+j)
               ch(2*i+1,k,1+j) = d39
               ch(2*i+1,k+1,1+j) = d40
               j14 = j14 + ido
            end do
         end do
      endif
      k5 = iand(max0(((l1 + 1 - j1)/2),0),3)
      j4 = 0        ! take care of compiler warnings
      do k = 1, k5
         j4 = 0
      end do
      do k = k5 + 1, ((l1 + 1 - j1)/2), 4
         j4 = 0
         j4 = 0
         j4 = 0
         j4 = 0
      end do
      do k3 = 0, ((l1 + 1 - j1)/2) - 1, 50
         k4 = min0(((l1 + 1 - j1)/2) - k3,50)
         do j10 = 0, ip - 2, 50
            j11 = min0(ip - 1 - j10,50)
            do k = 1, k4
               j15 = 0
               do j = 1, j11
                  ch(2*i,2*(k+k3)-1+j1,j10+1+j) = wa(ido*(j10+1)+2*i-1+
     1               is+j4+j15)*c1(2*i,2*(k+k3)-1+j1,j10+1+j) + wa(ido*(
     2               j10+1)+2*i+is+j4+j15)*c1(2*i+1,2*(k+k3)-1+j1,j10+1+
     3               j)
                  ch(2*i,2*(k+k3)+j1,j10+1+j) = wa(ido*(j10+1)+2*i-1+is+
     1               j4+j15)*c1(2*i,2*(k+k3)+j1,j10+1+j) + wa(ido*(j10+1
     2               )+2*i+is+j4+j15)*c1(2*i+1,2*(k+k3)+j1,j10+1+j)
                  ch(2*i+1,2*(k+k3)-1+j1,j10+1+j) = wa(ido*(j10+1)+2*i-1
     1               +is+j4+j15)*c1(2*i+1,2*(k+k3)-1+j1,j10+1+j) - wa(
     2               ido*(j10+1)+2*i+is+j4+j15)*c1(2*i,2*(k+k3)-1+j1,j10
     3               +1+j)
                  ch(2*i+1,2*(k+k3)+j1,j10+1+j) = wa(ido*(j10+1)+2*i-1+
     1               is+j4+j15)*c1(2*i+1,2*(k+k3)+j1,j10+1+j) - wa(ido*(
     2               j10+1)+2*i+is+j4+j15)*c1(2*i,2*(k+k3)+j1,j10+1+j)
                  j15 = j15 + ido
               end do
            end do
         end do
      end do
         endif
      end do
      go to 111
c      
  107 is = -ido
      do k = 1, l1
         do i = 1, (ido - 1)/2
            j5 = 0
      j16 = 0
      do j = 1, ip - 1
         ch(2*i,k,1+j) = wa(ido+2*i-1+is+j5+j16)*c1(2*i,k,1+j) + wa(ido+
     1      2*i+is+j5+j16)*c1(2*i+1,k,1+j)
         ch(2*i+1,k,1+j) = wa(ido+2*i-1+is+j5+j16)*c1(2*i+1,k,1+j) - wa(
     1      ido+2*i+is+j5+j16)*c1(2*i,k,1+j)
         j16 = j16 + ido
      end do
         end do
      end do
c      
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
      do i = 1, (ido - 1)/2
         c1(i*2,k,j) = ch(i*2,k,j) + ch(i*2,k,jc)
         c1(i*2,k,jc) = ch(1+i*2,k,j) - ch(1+i*2,k,jc)
         c1(1+i*2,k,j) = ch(1+i*2,k,j) + ch(1+i*2,k,jc)
         c1(1+i*2,k,jc) = ch(i*2,k,jc) - ch(i*2,k,j)
      end do
  113    continue
  114 continue
      go to 121
c      
  115 do 118 j=2,ipph
         jc = ipp2-j
      do k = 1, l1
      do i = 1, (ido - 1)/2
         c1(i*2,k,j) = ch(i*2,k,j) + ch(i*2,k,jc)
         c1(i*2,k,jc) = ch(1+i*2,k,j) - ch(1+i*2,k,jc)
         c1(1+i*2,k,j) = ch(1+i*2,k,j) + ch(1+i*2,k,jc)
         c1(1+i*2,k,jc) = ch(i*2,k,jc) - ch(i*2,k,j)
      end do
      end do
  118 continue
      go to 121
c      
  119 continue
      ik3 = iand(max0(idl1,0),3)
      ik8 = iand(max0(ik3,0),3)
      do ik = 1, ik8
         c2(ik,1) = ch2(ik,1)
      end do
      do ik = ik8 + 1, ik3, 4
         d69 = ch2(ik,1)
         d70 = ch2(ik+1,1)
         d71 = ch2(ik+2,1)
         d72 = ch2(ik+3,1)
         c2(ik,1) = d69
         c2(ik+1,1) = d70
         c2(ik+2,1) = d71
         c2(ik+3,1) = d72
      end do
      ik7 = iand(max0(((idl1 + 3 - ik3)/4)*4,0),3)
      do ik = 1, ik7
         c2(ik3+ik,1) = ch2(ik3+ik,1)
      end do
      do ik = ik7 + 1, ((idl1 + 3 - ik3)/4)*4, 4
         d65 = ch2(ik3+ik,1)
         d66 = ch2(ik3+ik+1,1)
         d67 = ch2(ik3+ik+2,1)
         d68 = ch2(ik3+ik+3,1)
         c2(ik3+ik,1) = d65
         c2(ik3+ik+1,1) = d66
         c2(ik3+ik+2,1) = d67
         c2(ik3+ik+3,1) = d68
      end do
c      
  121 do 123 j=2,ipph
         jc = ipp2-j
      do k = 1, l1
         c1(1,k,j) = ch(1,k,j) + ch(1,k,jc)
         c1(1,k,jc) = ch(1,k,jc) - ch(1,k,j)
      end do
  123 continue
c      
      ar1 = 1.d0
      ai1 = 0.d0
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
      do ik = 1, idl1
         ch2(ik,l) = c2(ik,1) + ar1*c2(ik,2)
         ch2(ik,lc) = ai1*c2(ik,ip)
      end do
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
      do ik = 1, idl1
         ch2(ik,l) = ch2(ik,l) + ar2*c2(ik,j)
         ch2(ik,lc) = ch2(ik,lc) + ai2*c2(ik,jc)
      end do
  126    continue
  127 continue
c      
      if (ipph - 1 .gt. 0) then
         j6 = iand(ipph - 1,3)
c      
         do j = 1, j6
            ik2 = iand(max0(idl1,0),3)
      ik6 = iand(max0(ik2,0),3)
      do ik = 1, ik6
         ch2(ik,1) = ch2(ik,1) + c2(ik,j+1)
      end do
      do ik = ik6 + 1, ik2, 4
         d61 = ch2(ik,1) + c2(ik,j+1)
         d62 = ch2(ik+1,1) + c2(ik+1,j+1)
         d63 = ch2(ik+2,1) + c2(ik+2,j+1)
         d64 = ch2(ik+3,1) + c2(ik+3,j+1)
         ch2(ik,1) = d61
         ch2(ik+1,1) = d62
         ch2(ik+2,1) = d63
         ch2(ik+3,1) = d64
      end do
      do ik = 1, ((idl1 + 3 - ik2)/4)
         d20 = ch2(ik2-2+ik*4,1) + c2(ik2-2+ik*4,j+1)
         d21 = ch2(ik2-1+ik*4,1) + c2(ik2-1+ik*4,j+1)
         d22 = ch2(ik2+ik*4,1) + c2(ik2+ik*4,j+1)
         ch2(ik2-3+ik*4,1) = ch2(ik2-3+ik*4,1) + c2(ik2-3+ik*4,j+1)
         ch2(ik2-2+ik*4,1) = d20
         ch2(ik2-1+ik*4,1) = d21
         ch2(ik2+ik*4,1) = d22
      end do
         end do
c      
         do j = j6 + 1, ipph - 1, 4
            ik1 = iand(max0(idl1,0),3)
      ik5 = iand(max0(ik1,0),3)
      do ik = 1, ik5
         ch2(ik,1) = ch2(ik,1) + c2(ik,j+1) + c2(ik,j+2) + c2(ik,j+3) + 
     1      c2(ik,j+4)
      end do
      do ik = ik5 + 1, ik1, 4
         d57=ch2(ik,1)+c2(ik,j+1)+c2(ik,j+2)+c2(ik,j+3)+c2(ik,j+4)
         d58 = ch2(ik+1,1) + c2(ik+1,j+1) + c2(ik+1,j+2) + c2(ik+1,j+3)
     1       + c2(ik+1,j+4)
         d59 = ch2(ik+2,1) + c2(ik+2,j+1) + c2(ik+2,j+2) + c2(ik+2,j+3)
     1       + c2(ik+2,j+4)
         d60 = ch2(ik+3,1) + c2(ik+3,j+1) + c2(ik+3,j+2) + c2(ik+3,j+3)
     1       + c2(ik+3,j+4)
         ch2(ik,1) = d57
         ch2(ik+1,1) = d58
         ch2(ik+2,1) = d59
         ch2(ik+3,1) = d60
      end do
      do ik = 1, ((idl1 + 3 - ik1)/4)
         d16 = ch2(ik1-2+ik*4,1) + c2(ik1-2+ik*4,j+1) + c2(ik1-2+ik*4,j+
     1      2) + c2(ik1-2+ik*4,j+3) + c2(ik1-2+ik*4,j+4)
         d17 = ch2(ik1-1+ik*4,1) + c2(ik1-1+ik*4,j+1) + c2(ik1-1+ik*4,j+
     1      2) + c2(ik1-1+ik*4,j+3) + c2(ik1-1+ik*4,j+4)
         d18 = ch2(ik1+ik*4,1) + c2(ik1+ik*4,j+1) + c2(ik1+ik*4,j+2) + 
     1      c2(ik1+ik*4,j+3) + c2(ik1+ik*4,j+4)
         ch2(ik1-3+ik*4,1) = ch2(ik1-3+ik*4,1) + c2(ik1-3+ik*4,j+1) + c2
     1      (ik1-3+ik*4,j+2) + c2(ik1-3+ik*4,j+3) + c2(ik1-3+ik*4,j+4)
         ch2(ik1-2+ik*4,1) = d16
         ch2(ik1-1+ik*4,1) = d17
         ch2(ik1+ik*4,1) = d18
      end do
         end do
      endif
c      
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
      i2 = iand(max0(ido,0),3)
      i6 = iand(max0(i2,0),3)
      do i = 1, i6
         cc(i,1,k) = ch(i,k,1)
      end do
      do i = i6 + 1, i2, 4
         d53 = ch(i,k,1)
         d54 = ch(i+1,k,1)
         d55 = ch(i+2,k,1)
         d56 = ch(i+3,k,1)
         cc(i,1,k) = d53
         cc(i+1,1,k) = d54
         cc(i+2,1,k) = d55
         cc(i+3,1,k) = d56
      end do
      i5 = iand(max0(((ido + 3 - i2)/4)*4,0),3)
      do i = 1, i5
         cc(i2+i,1,k) = ch(i2+i,k,1)
      end do
      do i = i5 + 1, ((ido + 3 - i2)/4)*4, 4
         d49 = ch(i2+i,k,1)
         d50 = ch(i2+i+1,k,1)
         d51 = ch(i2+i+2,k,1)
         d52 = ch(i2+i+3,k,1)
         cc(i2+i,1,k) = d49
         cc(i2+i+1,1,k) = d50
         cc(i2+i+2,1,k) = d51
         cc(i2+i+3,1,k) = d52
      end do
  131 continue
      go to 135
c      
  132 continue
      do k = 1, l1
         i1 = iand(max0(ido,0),3)
      i4 = iand(max0(i1,0),3)
      do i = 1, i4
         cc(i,1,k) = ch(i,k,1)
      end do
      do i = i4 + 1, i1, 4
         d45 = ch(i,k,1)
         d46 = ch(i+1,k,1)
         d47 = ch(i+2,k,1)
         d48 = ch(i+3,k,1)
         cc(i,1,k) = d45
         cc(i+1,1,k) = d46
         cc(i+2,1,k) = d47
         cc(i+3,1,k) = d48
      end do
      i3 = iand(max0(((ido + 3 - i1)/4)*4,0),3)
      do i = 1, i3
         cc(i1+i,1,k) = ch(i1+i,k,1)
      end do
      do i = i3 + 1, ((ido + 3 - i1)/4)*4, 4
         d41 = ch(i1+i,k,1)
         d42 = ch(i1+i+1,k,1)
         d43 = ch(i1+i+2,k,1)
         d44 = ch(i1+i+3,k,1)
         cc(i1+i,1,k) = d41
         cc(i1+i+1,1,k) = d42
         cc(i1+i+2,1,k) = d43
         cc(i1+i+3,1,k) = d44
      end do
      end do
c      
  135 continue
      do j7 = 0, ipph - 2, 50
         j8 = min0(ipph - 1 - j7,50)
         do k1 = 0, l1 - 1, 50
            k2 = min0(l1 - k1,50)
            do j = 1, j8
      do k = 1, k2
         cc(ido,2*(j+j7),k1+k) = ch(1,k1+k,j+1+j7)
         cc(1,2*(j+j7)+1,k1+k) = ch(1,k1+k,ipp2-1-j7-j)
      end do
            end do
         end do
      end do
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
      do i = 1, (ido - 1)/2
         cc(i*2,j2-1,k) = ch(i*2,k,j) + ch(i*2,k,jc)
         cc(idp2-(i+1)*2,j2-2,k) = ch(i*2,k,j) - ch(i*2,k,jc)
         cc(1+i*2,j2-1,k) = ch(1+i*2,k,j) + ch(1+i*2,k,jc)
         cc(idp2-1-i*2,j2-2,k) = ch(1+i*2,k,jc) - ch(1+i*2,k,j)
      end do
  139    continue
  140 continue
      return
c      
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
      do k = 1, l1
      do i = 1, (ido - 1)/2
         cc(i*2,j2-1,k) = ch(i*2,k,j) + ch(i*2,k,jc)
         cc(idp2-(i+1)*2,j2-2,k) = ch(i*2,k,j) - ch(i*2,k,jc)
         cc(1+i*2,j2-1,k) = ch(1+i*2,k,j) + ch(1+i*2,k,jc)
         cc(idp2-1-i*2,j2-2,k) = ch(1+i*2,k,jc) - ch(1+i*2,k,j)
      end do
      end do
  144 continue
c      
      return
      end
