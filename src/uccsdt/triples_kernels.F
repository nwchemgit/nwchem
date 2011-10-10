      subroutine w_zcji_zkab
     &  (w, zcji, zkab, 
     &   alo, ahi, blo, bhi, clo, chi, 
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zcji(zlo:zhi,clo:chi,jlo:jhi,ilo:ihi)
      double precision zkab(zlo:zhi,klo:khi,alo:ahi,blo:bhi)
      double precision r(zlo:zhi,alo:ahi,blo:bhi)
      double precision x(alo:ahi,blo:bhi,clo:chi,jlo:jhi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, jdim, jbegin, jend, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do k = klo, khi
c
         do b = blo, bhi
            do a = alo, ahi
               do z = zlo, zhi
                  r(z,a,b) = zkab(z,k,a,b)
               enddo
            enddo
         enddo
c
         ibegin = ilo
         if (spin.eq.spinck)ibegin = max(ilo,k+1)
         do i = ibegin, ihi
c
            jbegin = jlo
            if (spin.eq.spinck)jbegin = max(jlo,k+1)
            jend = min(jhi,i-1)
            jdim = jend - jbegin + 1
c
            if (jdim.gt.0)then
c
               call dgemm('T', 'N', adim*bdim, cdim*jdim, zdim, alpha, 
     &              r, zdim, zcji(zlo,clo,jbegin,i), zdim, 
     &              0.0d0, x(alo,blo,clo,jbegin), adim*bdim)
c
               do j = jbegin, jend
                  do a = alo, ahi
                     do b = blo, min(bhi,a-1)
                        cend = chi
                        if (spin.eq.spinck)cend = min(chi,b-1)
                        do c = clo, cend
                           w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(a,b,c,j)
                        end do
                     end do
                  end do
               end do
c     
               flopcount = flopcount + 
     &              2d0*dble(adim*bdim*cdim)*zdim*jdim
c
            endif
         end do
      end do
c
      call uccsdt_flops_exit(1, flopcount)
c
      end
      subroutine w_zcjk_ziab
     &  (w, zcjk, ziab,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zcjk(zlo:zhi,clo:chi,jlo:jhi,klo:khi)
      double precision ziab(zlo:zhi,ilo:ihi,alo:ahi,blo:bhi)
      double precision r(zlo:zhi,alo:ahi,blo:bhi)
      double precision x(alo:ahi,blo:bhi,clo:chi,jlo:jhi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, jdim, jbegin, jend
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do i = ilo, ihi
c
         do b = blo, bhi
            do a = alo, ahi
               do z = zlo, zhi
                  r(z,a,b) = ziab(z,i,a,b)
               enddo
            enddo
         enddo
c
         kend = khi
         if (spin.eq.spinck) kend = min(khi,i-1)
         do k = klo, kend
c
            jbegin = jlo
            if (spin.eq.spinck)jbegin = max(jlo,k+1)
            jend = min(jhi,i-1)
c
            jdim = jend - jbegin + 1
c
            if (jdim.gt.0)then
c
               call dgemm('T', 'N', adim*bdim, cdim*jdim, zdim, alpha, 
     &              r, zdim, zcjk(zlo,clo,jbegin,k), zdim, 
     &              0.0d0, x(alo,blo,clo,jbegin), adim*bdim)
c
               do j = jbegin, jend
                  do a = alo, ahi
                     do b = blo, min(bhi,a-1)
                        cend = chi
                        if (spin.eq.spinck)cend = min(chi,b-1)
                        do c = clo, cend
                           w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(a,b,c,j)
                        end do
                     end do
                  end do
               end do
c     
               flopcount = flopcount + 
     &              2d0*dble(adim*bdim*cdim)*zdim*jdim
c
            endif
         end do
      end do
c
      call uccsdt_flops_exit(2,flopcount)
c
      end
      subroutine w_zcik_zjab
     &  (w, zcik, zjab,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zcik(zlo:zhi,clo:chi,ilo:ihi,klo:khi)
      double precision zjab(zlo:zhi,jlo:jhi,alo:ahi,blo:bhi)
      double precision r(zlo:zhi,alo:ahi,blo:bhi)
      double precision x(alo:ahi,blo:bhi,clo:chi,ilo:ihi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, idim, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do j = jlo, jhi
c
         do b = blo, bhi
            do a = alo, ahi
               do z = zlo, zhi
                  r(z,a,b) = zjab(z,j,a,b)
               enddo
            enddo
         enddo
c
         kend = khi
         if (spin.eq.spinck) kend = min(khi,j-1)
         do k = klo, kend
c
            ibegin = max(ilo,j+1)
c
            idim = ihi - ibegin + 1
c
            if (idim.gt.0)then
c
               call dgemm('T', 'N', adim*bdim, cdim*idim, zdim, alpha, 
     &              r, zdim, zcik(zlo,clo,ibegin,k), zdim, 
     &              0.0d0, x(alo,blo,clo,ibegin), adim*bdim)
c
               do i = ibegin, ihi
                  do a = alo, ahi
                     do b = blo, min(bhi,a-1)
                        cend = chi
                        if (spin.eq.spinck)cend = min(chi,b-1)
                        do c = clo, cend
                           w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(a,b,c,i)
                        end do
                     end do
                  end do
               end do
c     
               flopcount = flopcount + 
     &              2d0*dble(adim*bdim*cdim)*zdim*idim
c
            endif
         end do
      end do
c
      call uccsdt_flops_exit(3, flopcount)
c
      end
      subroutine w_zaji_zkbc
     &  (w, zaji, zkbc,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zaji(zlo:zhi,alo:ahi,jlo:jhi,ilo:ihi)
      double precision zkbc(zlo:zhi,klo:khi,blo:bhi,clo:chi)
      double precision r(zlo:zhi,blo:bhi,clo:chi)
      double precision x(blo:bhi,clo:chi,alo:ahi,jlo:jhi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, jdim, jbegin, jend, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do k = klo, khi
c
         do c = clo, chi
            do b = blo, bhi
               do z = zlo, zhi
                  r(z,b,c) = zkbc(z,k,b,c)
               enddo
            enddo
         enddo
c
         ibegin = ilo
         if (spin.eq.spinck)ibegin = max(ilo,k+1)
         do i = ibegin, ihi
c
            jbegin = jlo
            if (spin.eq.spinck)jbegin = max(jlo,k+1)
            jend = min(jhi,i-1)
            jdim = jend - jbegin + 1
c
            if (jdim.gt.0)then
c
               call dgemm('T', 'N', bdim*cdim, adim*jdim, zdim, alpha, 
     &              r, zdim, zaji(zlo,alo,jbegin,i), zdim, 
     &              0.0d0, x(blo,clo,alo,jbegin), bdim*cdim)
c
               do j = jbegin, jend
                  do a = alo, ahi
                     do b = blo, min(bhi,a-1)
                        cend = chi
                        if (spin.eq.spinck)cend = min(chi,b-1)
                        do c = clo, cend
                           w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(b,c,a,j)
                        end do
                     end do
                  end do
               end do
c     
               flopcount = flopcount + 
     &              2d0*dble(adim*bdim*cdim)*zdim*jdim
c
            endif
         end do
      end do
c
      call uccsdt_flops_exit(4, flopcount)
c
      end
      subroutine w_zbji_zkac
     &  (w, zbji, zkac,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zbji(zlo:zhi,blo:bhi,jlo:jhi,ilo:ihi)
      double precision zkac(zlo:zhi,klo:khi,alo:ahi,clo:chi)
      double precision r(zlo:zhi,alo:ahi,clo:chi)
      double precision x(alo:ahi,clo:chi,blo:bhi,jlo:jhi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, jdim, jbegin, jend, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()

      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do k = klo, khi
c
         do c = clo, chi
            do a = alo, ahi
               do z = zlo, zhi
                  r(z,a,c) = zkac(z,k,a,c)
               enddo
            enddo
         enddo
c
         ibegin = ilo
         if (spin.eq.spinck)ibegin = max(ilo,k+1)
         do i = ibegin, ihi
c
            jbegin = jlo
            if (spin.eq.spinck)jbegin = max(jlo,k+1)
            jend = min(jhi,i-1)
            jdim = jend - jbegin + 1
c
            if (jdim.gt.0)then
c
               call dgemm('T', 'N', adim*cdim, bdim*jdim, zdim, alpha, 
     &              r, zdim, zbji(zlo,blo,jbegin,i), zdim, 
     &              0.0d0, x(alo,clo,blo,jbegin), adim*cdim)
c
               do j = jbegin, jend
                  do a = alo, ahi
                     do b = blo, min(bhi,a-1)
                        cend = chi
                        if (spin.eq.spinck)cend = min(chi,b-1)
                        do c = clo, cend
                           w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(a,c,b,j)
                        end do
                     end do
                  end do
               end do
c     
               flopcount = flopcount + 
     &              2d0*dble(adim*bdim*cdim)*zdim*jdim
c
            endif
         end do
      end do
c
      call uccsdt_flops_exit(5, flopcount)
c      
      end
      subroutine w_zakj_zicb
     &  (w, zakj, zicb,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zakj(zlo:zhi,alo:ahi,klo:khi,jlo:jhi)
      double precision zicb(zlo:zhi,ilo:ihi,clo:chi,blo:bhi)
      double precision r(zlo:zhi,clo:chi,blo:bhi)
      double precision x(clo:chi,blo:bhi,alo:ahi,klo:khi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, kdim
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do i = ilo, ihi
c
         do b = blo, bhi
            do c = clo, chi
               do z = zlo, zhi
                  r(z,c,b) = zicb(z,i,c,b)
               enddo
            enddo
         enddo
         do j = jlo, min(jhi,i-1)
            kend = khi
            if (spin.eq.spinck) kend = min(khi,j-1)
            kdim = kend - klo + 1
c
            call dgemm('T', 'N', cdim*bdim, adim*kdim, zdim, alpha, 
     &           r, zdim, zakj(zlo,alo,klo,j), zdim, 
     &           0.0d0, x, cdim*bdim)
c
            do k = klo, kend
               do a = alo, ahi
                  do b = blo, min(bhi,a-1)
                     cend = chi
                     if (spin.eq.spinck)cend = min(chi,b-1)
                     do c = clo, cend
                        w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(c,b,a,k)
                     end do
                  end do
               end do
            end do
c
            flopcount = flopcount + 2d0*dble(adim*bdim*cdim)*zdim*kdim
c
         end do
      end do
      call uccsdt_flops_exit(6, flopcount)
c
      end
      subroutine w_zbkj_zica
     &  (w, zbkj, zica,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zbkj(zlo:zhi,blo:bhi,klo:khi,jlo:jhi)
      double precision zica(zlo:zhi,ilo:ihi,clo:chi,alo:ahi)
      double precision r(zlo:zhi,clo:chi,alo:ahi)
      double precision x(clo:chi,alo:ahi,blo:bhi,klo:khi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, kdim
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do i = ilo, ihi
c
         do a = alo, ahi
            do c = clo, chi
               do z = zlo, zhi
                  r(z,c,a) = zica(z,i,c,a)
               enddo
            enddo
         enddo
         do j = jlo, min(jhi,i-1)
            kend = khi
            if (spin.eq.spinck) kend = min(khi,j-1)
            kdim = kend - klo + 1
c
            call dgemm('T', 'N', cdim*adim, bdim*kdim, zdim, alpha, 
     &           r, zdim, zbkj(zlo,blo,klo,j), zdim, 
     &           0.0d0, x, cdim*adim)
c
            do k = klo, kend
               do a = alo, ahi
                  do b = blo, min(bhi,a-1)
                     cend = chi
                     if (spin.eq.spinck)cend = min(chi,b-1)
                     do c = clo, cend
                        w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(c,a,b,k)
                     end do
                  end do
               end do
            end do
c
            flopcount = flopcount + 2d0*dble(adim*bdim*cdim)*zdim*kdim
c
         end do
      end do
c
      call uccsdt_flops_exit(7, flopcount)
c
      end
      subroutine w_zaki_zjcb
     &  (w, zaki, zjcb,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zaki(zlo:zhi,alo:ahi,klo:khi,ilo:ihi)
      double precision zjcb(zlo:zhi,jlo:jhi,clo:chi,blo:bhi)
      double precision r(zlo:zhi,clo:chi,blo:bhi)
      double precision x(clo:chi,blo:bhi,alo:ahi,klo:khi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, kdim, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do j = jlo, jhi
c
         do b = blo, bhi
            do c = clo, chi
               do z = zlo, zhi
                  r(z,c,b) = zjcb(z,j,c,b)
               enddo
            enddo
         enddo
c
         ibegin = max(ilo,j+1)
         do i = ibegin, ihi
c
            kend = khi
            if (spin.eq.spinck) kend = min(khi,j-1)
            kdim = kend - klo + 1
c
            call dgemm('T', 'N', cdim*bdim, adim*kdim, zdim, alpha, 
     &           r, zdim, zaki(zlo,alo,klo,i), zdim, 
     &           0.0d0, x, cdim*bdim)
c
            do k = klo, kend
               do a = alo, ahi
                  do b = blo, min(bhi,a-1)
                     cend = chi
                     if (spin.eq.spinck)cend = min(chi,b-1)
                     do c = clo, cend
                        w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(c,b,a,k)
                     end do
                  end do
               end do
            end do
c
            flopcount = flopcount + 2d0*dble(adim*bdim*cdim)*zdim*kdim
c
         end do
      end do
c
      call uccsdt_flops_exit(8, flopcount)
c
      end
      subroutine w_zbki_zjca
     &  (w, zbki, zjca,
     &   alo, ahi, blo, bhi, clo, chi,
     &   ilo, ihi, jlo, jhi, klo, khi,
     &   zlo, zhi, spin, spinck, r, x, alpha)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     &        ilo, ihi, jlo, jhi, klo, khi,
     &        zlo, zhi, spin, spinck
c
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision zbki(zlo:zhi,blo:bhi,klo:khi,ilo:ihi)
      double precision zjca(zlo:zhi,jlo:jhi,clo:chi,alo:ahi)
      double precision r(zlo:zhi,clo:chi,alo:ahi)
      double precision x(clo:chi,alo:ahi,blo:bhi,klo:khi)
      double precision alpha
      integer i, j, k, a, b, c, z, cend, kend
      integer adim, bdim, cdim, zdim, kdim, ibegin
      double precision flopcount
c
      flopcount = 0d0
      call uccsdt_flops_enter()
c
      adim = ahi - alo + 1
      bdim = bhi - blo + 1
      cdim = chi - clo + 1
      zdim = zhi - zlo + 1
c
      do j = jlo, jhi
c     
         do a = alo, ahi
            do c = clo, chi
               do z = zlo, zhi
                  r(z,c,a) = zjca(z,j,c,a)
               enddo
            enddo
         enddo
c
         ibegin = max(ilo,j+1)
         do i = ibegin, ihi
c
            kend = khi
            if (spin.eq.spinck) kend = min(khi,j-1)
            kdim = kend - klo + 1
c
            call dgemm('T', 'N', cdim*adim, bdim*kdim, zdim, alpha, 
     &           r, zdim, zbki(zlo,blo,klo,i), zdim, 
     &           0.0d0, x, cdim*adim)
c
            do k = klo, kend
               do a = alo, ahi
                  do b = blo, min(bhi,a-1)
                     cend = chi
                     if (spin.eq.spinck)cend = min(chi,b-1)
                     do c = clo, cend
                        w(c,b,a,k,j,i) = w(c,b,a,k,j,i) + x(c,a,b,k)
                     end do
                  end do
               end do
            end do
c
            flopcount = flopcount + 2d0*dble(adim*bdim*cdim)*zdim*kdim
c
         end do
      end do
c
      call uccsdt_flops_exit(9, flopcount)
c
      end
      subroutine w_energy
     $  (w, t1, t1_ck, eps, energy,
     $   int_abij, int_abjk, int_abik, int_bcij,
     $   int_acij, int_bcjk, int_acjk, int_bcik,
     $   int_acik,
     $   alo, ahi, blo, bhi, clo, chi,
     $   ilo, ihi, jlo, jhi, klo, khi,
     $   vlo, vhi, olo, ohi, 
     $   vlo_ck, vhi_ck, olo_ck, ohi_ck, 
     $   spin, spinck, nmo)
c
      implicit none
      integer alo, ahi, blo, bhi, clo, chi,
     $        ilo, ihi, jlo, jhi, klo, khi,
     $        vlo, vhi, olo, ohi, 
     $        vlo_ck, vhi_ck, olo_ck, ohi_ck, 
     $        spin, spinck, nmo
c
      double precision t1(vlo:vhi,olo:ohi)
      double precision t1_ck(vlo_ck:vhi_ck,olo_ck:ohi_ck)
      double precision w(clo:chi,blo:bhi,alo:ahi,
     &                   klo:khi,jlo:jhi,ilo:ihi)
      double precision int_abij(alo:ahi,blo:bhi,ilo:ihi,jlo:jhi)
      double precision int_abjk(alo:ahi,blo:bhi,jlo:jhi,klo:khi)
      double precision int_abik(alo:ahi,blo:bhi,ilo:ihi,klo:khi)
      double precision int_bcij(blo:bhi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_acij(alo:ahi,clo:chi,ilo:ihi,jlo:jhi)
      double precision int_bcjk(blo:bhi,clo:chi,jlo:jhi,klo:khi)
      double precision int_acjk(alo:ahi,clo:chi,jlo:jhi,klo:khi)
      double precision int_bcik(blo:bhi,clo:chi,ilo:ihi,klo:khi)
      double precision int_acik(alo:ahi,clo:chi,ilo:ihi,klo:khi)
      double precision eps(nmo,2), energy, v, d
      integer i, j, k, a, b, c, ctop, ktop
      do a = alo, ahi
         do b = blo, min(bhi,a-1)
           ctop = chi
           if (spin.eq.spinck) ctop = min(chi,b-1)
            do c = clo, ctop
               do i = ilo, ihi
                  do j = jlo, min(jhi,i-1)
                    ktop = khi
                    if (spin.eq.spinck) ktop = min(khi,j-1)
                     do k = klo, ktop
c
                        v = w(c,b,a,k,j,i)
     $                        + t1_ck(c,k)*
     $                          int_abij(a,b,i,j)
     $                        + t1(a,i)*
     $                          int_bcjk(b,c,j,k)
     $                        - t1(b,i)*
     $                          int_acjk(a,c,j,k)
     $                        - t1(a,j)*
     $                          int_bcik(b,c,i,k)
     $                        + t1(b,j)*
     $                          int_acik(a,c,i,k)
c
                        if (spin .eq. spinck) then
                          v = v 
     $                        + t1(c,i)*
     $                          int_abjk(a,b,j,k)
     $                        - t1(c,j)*
     $                          int_abik(a,b,i,k)
     $                        + t1(a,k)*
     $                          int_bcij(b,c,i,j)
     $                        - t1(b,k)*
     $                          int_acij(a,c,i,j)
                        end if
c
                        d = eps(a,spin)+
     $                      eps(b,spin)+
     $                      eps(c,spinck)-
     $                      eps(i,spin)-
     $                      eps(j,spin)-
     $                      eps(k,spinck)
c
                        energy = energy - 
     $                           v*w(c,b,a,k,j,i)/d
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      end
c $Id$
