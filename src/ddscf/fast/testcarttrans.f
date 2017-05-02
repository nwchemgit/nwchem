      implicit none
*
* $Id$
*
      integer lmax, npoly, numl, lmax2
      double precision util_random
      external util_random
      parameter (lmax = 3, npoly = 7)
      parameter (lmax2 = lmax+lmax)
      parameter (numl = (lmax+1)*(lmax+2)*(lmax+3)/6)
      integer i, j, k, l, m, ijk, ind, nlm, n
      integer i1, i2, j1, j2, k1, k2, l1, l2, l1p, l2p
      double precision x1, x2, y1, y2, z1, z2
      double precision test1, test2, test3, r, factor
      double precision  coeff(npoly,numl)
      double precision scoeff(npoly,numl)
      double precision ncoeff(npoly,numl)
      double precision dens(numl*numl),ndens(numl*numl),work(numl*numl)
      double precision junk, a, b, c, x, y, z, test, testn, err
      double precision d(((lmax2+1)*(lmax2+2))/2, -lmax2:lmax2, 0:lmax2)
      double precision 
     $     dinv(((lmax2+1)*(lmax2+2)*(lmax2+3))/6, -lmax2:lmax2,0:lmax2)
      double precision q(-lmax2:lmax2,0:lmax2)
c
      junk = util_random(55512121)
      call xlm_init
      call xlm_coeff_inv(lmax2,d,dinv)
c
      do l = 0, lmax
         do m = 1, npoly
            ijk = (l*(l+1)*(l+2))/6
            do i = l,0,-1
               do j = l-i,0,-1
                  k = l-i-j
                  ijk = ijk + 1
                  scoeff(m,ijk) = util_random(0) - 0.5d0
                  coeff(m,ijk) = scoeff(m,ijk)
               enddo
            enddo
         enddo
      enddo
c
      a = util_random(0)-0.5d0
      b = util_random(0)-0.5d0
      c = util_random(0)-0.5d0
      x = util_random(0)-0.5d0
      y = util_random(0)-0.5d0
      z = util_random(0)-0.5d0
      write(6,1) ' New center ', a, b, c
      write(6,1) ' Target     ', x, y, z
 1    format(a,3f12.6)
c
      call cart_poly_translate(lmax,npoly,coeff,ncoeff,a,b,c)
c
      err = 0.0d0
      do m = 1, npoly
         test = 0.0d0
         testn= 0.0d0
         ijk = 0
         do l = 0, lmax
            do i = l,0,-1
               do j = l-i,0,-1
                  k = l-i-j
                  ijk = ijk + 1
                  test = test + scoeff(m,ijk)*x**i*y**j*z**k
                  testn = testn + ncoeff(m,ijk)*
     $                 (x-a)**i*(y-b)**j*(z-c)**k
               enddo
            enddo
         enddo
         err = max(err,abs(test-testn))
**         write(6,*) m, test, testn, test-testn
      enddo
      write(6,*) ' Err from translation of poly ', err
c
      do i = 1, numl**2
         dens(i) = util_random(0) - 0.5d0
      enddo
      x1 = util_random(0)-0.5d0
      y1 = util_random(0)-0.5d0
      z1 = util_random(0)-0.5d0
      x2 = util_random(0)-0.5d0
      y2 = util_random(0)-0.5d0
      z2 = util_random(0)-0.5d0
c
      write(6,1) ' New center ', a, b, c
      write(6,1) ' Center1    ', x1, y1, z1
      write(6,1) ' Center2    ', x2, y2, z2
      err = 0.0d0
      do l1 = 0, lmax
         do l2 = 0, lmax
c
            test1 = 0.0d0
            ind = 0
            do i2 = l2,0,-1
               do j2 = l2-i2,0,-1
                  k2 = l2-i2-j2
                  do i1 = l1,0,-1
                     do j1 = l1-i1,0,-1
                        k1 = l1-i1-j1
                        ind = ind + 1
                        test1 = test1 + dens(ind)*
     $                       (x-x1)**i1*(y-y1)**j1*(z-z1)**k1*
     $                       (x-x2)**i2*(y-y2)**j2*(z-z2)**k2
                     enddo
                  enddo
               enddo
            enddo
c                        
            call cart_dens_translate(l1,x1,y1,z1,l2,x2,y2,z2,
     $           a,b,c,dens,ndens,work)
c
            test2 = 0.0d0
            ind = 0
            do i2 = l2,0,-1
               do j2 = l2-i2,0,-1
                  k2 = l2-i2-j2
                  do i1 = l1,0,-1
                     do j1 = l1-i1,0,-1
                        k1 = l1-i1-j1
                        ind = ind + 1
                        test2 = test2 + dens(ind)*
     $                       (x-x1)**i1*(y-y1)**j1*(z-z1)**k1*
     $                       (x-x2)**i2*(y-y2)**j2*(z-z2)**k2
                     enddo
                  enddo
               enddo
            enddo
c                        
            test3 = 0.0d0
            ind = 0
            do l2p = 0, l2
               do i2 = l2p,0,-1
                  do j2 = l2p-i2,0,-1
                     k2 = l2p-i2-j2
                     do l1p = 0, l1
                        do i1 = l1p,0,-1
                           do j1 = l1p-i1,0,-1
                              k1 = l1p-i1-j1
                              ind = ind + 1
                              test3 = test3 + ndens(ind)*
     $                             (x-a)**i1*(y-b)**j1*(z-c)**k1*
     $                             (x-a)**i2*(y-b)**j2*(z-c)**k2
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
c
            err = max(err,abs(test3-test1))
**            write(6,*) l1, l2, test1, test3, test3-test1
c
            call cart_dens_product(l1,l2,ndens,work)
            ijk = 0
            test = 0.0d0
            do l = 0, l1+l2
               do i = l,0,-1
                  do j = l-i,0,-1
                     k = l-i-j
                     ijk = ijk + 1
                     test = test + work(ijk)*
     $                    (x-a)**i*(y-b)**j*(z-c)**k
                  enddo
               enddo
            enddo
c
            err = max(err,abs(test-test1))
**            write(6,*) test1-test
c
            call cart_dens_to_sph(l1+l2,work,ndens,dinv,lmax2)
            test = 0.0d0
            r = sqrt((x-a)**2 + (y-b)**2 + (z-c)**2)
            call xlm(l1+l2,x-a,y-b,z-c,q,lmax2)
            nlm = 0
            do n = 0, l1+l2
               do l = n, 0, -2
                  factor = r**(n-l)
                  do m = -l, l
                     nlm = nlm + 1
                     test = test + factor*q(m,l)*ndens(nlm)
                  enddo
               enddo
            enddo
            err = max(err,abs(test-test1))
**            write(6,*) test1-test
c
            call cart_dens_trans_prod_sph(
     $           l1, x1, y1, z1,
     $           l2, x2, y2, z2, 
     $           a, b, c, work, dinv, lmax2, dens, ndens)

            test = 0.0d0
            r = sqrt((x-a)**2 + (y-b)**2 + (z-c)**2)
            call xlm(l1+l2,x-a,y-b,z-c,q,lmax2)
            nlm = 0
            do n = 0, l1+l2
               do l = n, 0, -2
                  factor = r**(n-l)
                  do m = -l, l
                     nlm = nlm + 1
                     test = test + factor*q(m,l)*ndens(nlm)
                  enddo
               enddo
            enddo
            err = max(err,abs(test-test1))
**            write(6,*) test1-test
         enddo
      enddo
      write(6,*) ' Error from translation of density ', err
c
      end
