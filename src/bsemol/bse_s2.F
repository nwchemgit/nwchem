      subroutine bse_s2(nocc, nmo, npoles, xup, xdw, yup, ydw, 
     $                  ovlp, s2)
      implicit none

      integer nocc(2), npoles(2), nmo
      double precision xup(npoles(1)), xdw(npoles(2))
      double precision yup(npoles(1)), ydw(npoles(2))
      double precision s2
      double precision ovlp(nmo,nmo)

      double precision ov,xx,yy
      integer i,j,a,b,ipole,jpole,nvir(2)
      integer llipole,ulipole,lljpole,uljpole

      double precision, external :: ddot

      s2 = 0d0
      nvir(1) = nmo - nocc(1)
      nvir(2) = nmo - nocc(2)

      ! Up-Up
      do i=1,nocc(1)
        do j=1,nocc(1)
          llipole = (i-1)*nvir(1) + 1
          ulipole = i*nvir(1)
          lljpole = (j-1)*nvir(1) + 1
          uljpole = j*nvir(1)
          ov = dot_product(ovlp(i,1:nocc(2)),ovlp(j,1:nocc(2)))
          xx = dot_product(xup(llipole:ulipole),xup(lljpole:uljpole))

          yy = dot_product(yup(llipole:ulipole),yup(lljpole:uljpole))
          s2 = s2 + ov*(xx+yy)
        enddo
      enddo
      do a=nocc(1)+1,nmo
        do b=nocc(1)+1,nmo
          ov = dot_product(ovlp(a,1:nocc(2)),ovlp(b,1:nocc(2)))
          llipole = a - nocc(1)
          lljpole = b - nocc(1)
          s2 = s2 - ov*ddot(nocc(1),xup(llipole),nvir(1),
     $                              xup(lljpole),nvir(1))
          s2 = s2 - ov*ddot(nocc(1),yup(llipole),nvir(1),
     $                              yup(lljpole),nvir(1))
        enddo
      enddo

      ! Dw-Dw
      do i=1,nocc(2)
        do j=1,nocc(2)
          llipole = (i-1)*nvir(2) + 1
          ulipole = i*nvir(2)
          lljpole = (j-1)*nvir(2) + 1
          uljpole = j*nvir(2)
          ov = dot_product(ovlp(1:nocc(1),i),ovlp(1:nocc(1),j))
          xx = dot_product(xdw(llipole:ulipole),xdw(lljpole:uljpole))
          yy = dot_product(ydw(llipole:ulipole),ydw(lljpole:uljpole))
          s2 = s2 + ov*(xx+yy)
        enddo
      enddo
      do a=nocc(2)+1,nmo
        do b=nocc(2)+1,nmo
          ov = dot_product(ovlp(1:nocc(1),a),ovlp(1:nocc(1),b))
          llipole = a - nocc(2)
          lljpole = b - nocc(2)
          s2 = s2 - ov*ddot(nocc(2),xdw(llipole),nvir(2),
     $                              xdw(lljpole),nvir(2))
          s2 = s2 - ov*ddot(nocc(2),ydw(llipole),nvir(2),
     $                              ydw(lljpole),nvir(2))
        enddo
      enddo

      ! Up-Dw
      do ipole=1,npoles(1)
        i = ipole/nvir(1) + 1
        a = ipole - (i-1)*nvir(1) + nocc(1)
        do jpole=1,npoles(2)
          j = jpole/nvir(2) + 1
          b = jpole - (j-1)*nvir(2) + nocc(2)
          s2 = s2 - 
     $     2d0*xup(ipole)*xdw(jpole)*ovlp(i,j)*ovlp(a,b) -
     $     2d0*yup(ipole)*ydw(jpole)*ovlp(i,j)*ovlp(a,b) +
     $     2d0*xup(ipole)*ydw(jpole)*ovlp(a,j)*ovlp(i,b)
        enddo
      enddo

      ! Dw-Up
      do ipole=1,npoles(2)
        i = ipole/nvir(2) + 1
        a = ipole - (i-1)*nvir(2) + nocc(2)
        do jpole=1,npoles(1)
          j = jpole/nvir(1) + 1
          b = jpole - (j-1)*nvir(1) + nocc(1)
          s2 = s2 + 
     $      2d0*xdw(ipole)*yup(jpole)*ovlp(j,a)*ovlp(b,i)
        enddo
      enddo

      end subroutine

