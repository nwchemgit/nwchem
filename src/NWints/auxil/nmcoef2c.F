      subroutine nmcoef2c(z,cc,l,np,scr,lscr)
c $Id$
c.. routine to compute normalized contraction based on (cgtf|cgtf) 
c..  normalization factor and do normalization
c
      implicit none
#include "errquit.fh"
c::passed
      integer l                   ! [input] angular momentum type
      integer np                  ! [input] number of primitives in contraction
      double precision z(np)      ! [input] exponents
      double precision cc(np)     ! [input/output] contraction coefficients
      integer lscr                ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
c::local
      double precision factor, sum, almost_one
      double precision xyz(3), dumexp, dumcof
      double precision PI
      integer nint, i, nbf, l_lim
      integer ip, j, index
      logical warning
      logical oldway
c
      oldway = .false.
c
      PI = 2.0d00*acos(0.0d00)
c
c... dummy coordinates
      xyz(1) = 0.0d00
      xyz(2) = 0.0d00
      xyz(3) = 0.0d00
c... dummy coef/exponent
      dumexp = 0.0d00
      dumcof = 1.0d00
c... number of integrals
      if (l.lt.0) call errquit
     &    ('nmcoef2c: how did I get an sp or spd shell here?',911,
     &       BASIS_ERR)
      nbf = (l+1)*(l+2)/2
      nint = nbf*nbf
      if (.not.oldway) then
*
*...
*... normalize for overlap of basis function in all space
*... e.g., normalize with respect to primitive coefficients
*...
*
         factor = 1.0d00
         l_lim = 2*l - 1
         do i = l_lim,2,-2
            factor = factor*dble(i)
         enddo
         if(l.eq.0) then
            do ip = 1,np
               cc(ip) = cc(ip) *
     &              ((2.D0*z(ip)/PI)**0.75D0)/ sqrt(factor)
            enddo
         else
            do ip = 1,np
               cc(ip) = cc(ip) *
     &              ((2.D0*z(ip)/PI)**0.75D0)*(4.D0*z(ip))**(0.5D0*l)
               cc(ip) = cc(ip) / sqrt(factor)
            enddo
         endif
*     
*...
* normalize contraction coefficients to make (a|a) = 1 
*...
*
        sum = 0.0d00
        do i = 1,np
          do j = i,np
*debug:            write(6,*)' '
*debug:            write(6,*)'---------------------------------------'
*debug:            write(6,*)' xyz ',xyz
*debug:            write(6,*)' i,exp,coeff: ',i,z(i),cc(i)
*debug:            write(6,*)' j,exp,coeff: ',j,z(j),cc(j)
            call hf2(
     &          xyz,z(i),cc(i),1,1,l,
     &          xyz,dumexp,dumcof,1,1,0,
     &          xyz,z(j),cc(j),1,1,l,
     &          xyz,dumexp,dumcof,1,1,0,
     &          scr,nint,.false.,.false.,.false.,.false.,
     &          scr(1+nint),(lscr-nint))
* find minimum diagonal element
            call find_min_diag_nmcoef2c(scr,nbf,factor)
*debug:            write(6,*)' sum 1     ', sum, factor
            if(i.ne.j) factor = factor*2d0
            sum = sum + factor
*debug:            write(6,*)' sum 2     ', sum
          enddo
        enddo
*debug:        write(6,*)' sum final ', sum
        factor = 1.0d00/sqrt(sum)
*debug:        write(6,*)' factor    ',factor
        do i = 1,np
          cc(i) = cc(i)*factor
        enddo
c
      else
        call hf2 (
     &        xyz,z,cc,np,1,l,
     &        xyz,dumexp,dumcof,1,1,0,
     &        xyz,z,cc,np,1,l,
     &        xyz,dumexp,dumcof,1,1,0,
     &        scr,nint,.false.,.false.,.false.,.false.,
     &        scr(1+nint),(lscr-nint))
*        
crak:--
c.. find minimum diagonal element
c
        factor = abs(scr(1))
        index = 0
        do i = 1,nbf
          do j = 1,nbf
            index = index + 1
            if (i.eq.j) then
              if (abs(scr(index)).gt.1.0d-05) then
                factor = min(factor,abs(scr(index)))
              endif
*debug:              write(6,*)'factor ',factor,' val ',scr(index)
            endif
          enddo
        enddo
        factor = 1.0d00/sqrt(factor)
c
c... normalize contraction coefficients
        do i = 1, np
          cc(i) = cc(i)*factor
        enddo
      endif
*debug:* compute sum of squares
*debug:      factor = 0.0d00
*debug:      do ip = 1,np
*debug:        factor = factor + cc(ip)*cc(ip)
*debug:      enddo
*debug:      write(6,*) ' sum of squares :',factor
*
*...
*...  compute integrals and check diagonal elements
*...  diagonal elements must be >= 1.0
*
      call hf2 (
     &        xyz,z,cc,np,1,l,
     &        xyz,dumexp,dumcof,1,1,0,
     &        xyz,z,cc,np,1,l,
     &        xyz,dumexp,dumcof,1,1,0,
     &        scr,nint,.false.,.false.,.false.,.false.,
     &        scr(1+nint),(lscr-nint))
      warning = .false.
      almost_one = 1.0d00 - 1.0d-10  ! allow a little slop
00001 continue
      index = 0
      do i = 1,nbf
        do j = 1,nbf
          index = index + 1
          if (i.eq.j) then
            factor = scr(index)
            if ((.not.warning).and.(factor.lt.almost_one)) then
              warning = .true.
              write(6,*)' nmcoef2c: problem normalizing basis set ',
     &            '::factor       ::', factor
              factor = 1.0d00 - factor
              write(6,*)' nmcoef2c: problem normalizing basis set ',
     &            '::diff from 1.0::', factor
              goto 00001
            endif
            if (warning) then
              write(6,10000)l,i,j,factor
            endif
          endif
        enddo
      enddo
      if (warning)
     &    call errquit(' nmcoef2c: problem normalizing basis set',911,
     &       BASIS_ERR)
10000 format(1x,'ints after normalization(',i3,
     &    ') <',i5,'|',i5,'> =',1pd20.10)
      end
*
      subroutine find_min_diag_nmcoef2c(buf,nr,factor)
      implicit none
      integer nr
      double precision buf(nr,nr)
      double precision factor
      logical negative
c
      integer i,j
*debug:      write(6,*)' nr = ',nr
      factor = abs(buf(1,1))
      negative = buf(1,1).lt.0.0d00
      do i = 1,nr
        do j = 1,nr
          if (i.eq.j) then
*debug:            write(6,10000)i,j,buf(i,j)
            if (abs(buf(i,j)).lt.factor) then
              factor = abs(buf(i,j))
              negative = buf(i,j) .lt. 0.0d00
            endif
          endif
        enddo
      enddo
      if (negative) factor = -1.0d00*factor
10000 format(' findbuffer(',i3,',',i3,') =',1pd20.10)
      end
