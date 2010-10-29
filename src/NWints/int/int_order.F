* $Id$
      program int_order
      implicit none
c
      integer icase
c
      external int_order_BD
c
      call defNxyz(6)
c
00001 continue
      write(6,*)' '
      write(6,*)' '
      write(6,*)' enter a choice '
      write(6,*)' 1 = (overlap, kinetic, potential) 1e integrals '
      write(6,*)' 2 = 4 center 2 electron integrals ' 
      write(6,*)' 3 = 3 center overlap integrals '
      write(6,*)' 4 = dipole integrals     '
      write(6,*)' 5 = 2 center 2 electron integrals '
      write(6,*)' 6 = 3 center 2 electron integrals '
      write(6,*)' 0 = exit'
      read(5,*)icase
      if (icase.eq.0) then
        goto 99999
      elseif (icase.eq.1) then
        call int_case1_order()
        goto 00001
      elseif (icase.eq.2) then
        write(6,*)' not implemented yet'
        goto 00001
      elseif (icase.eq.3) then
        write(6,*)' not implemented yet'
        goto 00001
      elseif (icase.eq.4) then
        call int_case4_order()
        goto 00001
      elseif (icase.eq.5) then
        write(6,*)' not implemented yet'
        goto 00001
      elseif (icase.eq.6) then
        write(6,*)' not implemented yet'
        goto 00001
      else
        write(6,*)' invalid case ', icase
        goto 00001
      endif
99999 continue
*rak:      write(6,*)' testing now'
*rak:      do lval=0,6
*rak:        l2 = (lval+1)*(lval+2)/2
*rak:        write(6,*)' cartesan lval=',lval
*rak:        do indd = 1,l2
*rak:          call int_str_ang(lval,indd,.false.,string,lstr)
*rak:          write(6,*)' len of str ',lstr,' <',string(1:lstr),'>'
*rak:        enddo
*rak:        write(6,*)' spherical lval=',lval
*rak:        l2 = 2*lval+1
*rak:        do indd = 1,l2
*rak:          call int_str_ang(lval,indd,.true.,string,lstr)
*rak:          write(6,*)' len of str ',lstr,' <',string(1:lstr),'>'
*rak:        enddo
*rak:      enddo
c
      end
      subroutine int_str_ang(lval,indx,ospherical,str,lstr)
      implicit none
*
* routine returns printable string for lval and index of lval
* lval is the angular momentum and index is the particular element
* of the array.
* for example: lval = 1 is a p function therefore there are 3 elements
*              so valid indx values are 1, 2, or 3
*
#include "inp.fh"
#include "int_order.fh"
#include "errquit.fh"
*::passed
      integer lval    ! [input] lvalue
      integer indx    ! [input] index value of angular block
      logical ospherical ! [input] is the integral block spherical
      character*(*) str  ! [output] string associated with this
      integer lstr       ! [output] actual length of the string from inp_strlen
*::local
      integer l2, strval, i
c
      if (lval.lt.-1.or.lval.gt.6) call errquit
     &    ('int_str_ang: invalid l-value of ',lval,
     &       INT_ERR)
      str = ' '
c
      if (lval.eq.-1) then
        if (indx.eq.1) then
          str = 's'
        else if (indx.ge.2.and.indx.le.4) then
          call int_str_xyz(lval,(indx-1),str,lstr)
        else
          write(6,*)' invalid index for lval of -1'
          write(6,*)' indx ',indx
          call errquit(' ',911, INT_ERR)
        endif
        lstr = inp_strlen(str)
        return
      endif
      if (ospherical) then
        l2 = 2*lval+1
        if (.not.(indx.ge.1.and.indx.le.l2)) then
          write(6,*)' lval   = ',lval
          write(6,*)' indx   = ',indx
          call errquit('int_str_ang: indx invalid',911, INT_ERR)
        endif
        strval = 0 - lval - 1  + indx
        if (.not.(strval.ge.(-lval).and.strval.le.lval)) then
          write(6,*)' lval   = ',lval
          write(6,*)' indx   = ',indx
          write(6,*)' strval = ',strval
          call errquit ('int_str_ang: strval invalid ',911, INT_ERR)
        endif
        str = shell_types(lval) 
        write(str(3:4),'(i2)')strval
        do i = 5,len(str)
          str(i:i) = ' '
        enddo
        lstr = inp_strlen(str)
      else
        l2 = (lval+1)*(lval+2)/2
        if (.not.(indx.ge.1.and.indx.le.l2)) then
          write(6,*)' lval   = ',lval
          write(6,*)' indx   = ',indx
          call errquit('int_str_ang: indx invalid',911, INT_ERR)
        endif
        call int_str_xyz(lval,indx,str,lstr)
      endif
c
      end
      subroutine int_str_xyz(lval,indx,str,lstr)
      implicit none
#include "inp.fh"
      integer lval
      integer indx
      character*(*) str
      integer lstr
c
      integer nxyz(3)
c
      character*1 pxyz(3)
      integer i,j, icnt
c
      if (lval.eq.0) then
        str = 's'
        lstr = inp_strlen(str)
        return
      endif
      call getNxyz(lval,indx,nxyz)
      pxyz(1) = 'x'
      pxyz(2) = 'y'
      pxyz(3) = 'z'
      str=' '
      icnt = 1
      do i=1,3
        do j=1,nxyz(i)
          str(icnt:icnt) = pxyz(i)
          icnt = icnt + 1
        enddo
      enddo
      do i = icnt,len(str)
        str(i:i) = ' '
      enddo
      lstr = inp_strlen(str)
c
      end
      block data int_order_BD
#include "int_order.fh"      
*            
*                        -1   0    1    2    3    4    5    6    7   8
      data shell_types /'sp','s ','p ','d ','f ','g ','h ','i ','k','l'/
      end
      subroutine int_case1_order()
      implicit none
c
      integer int_bra_sph, bra_type, l2_bra, ibra, l_bra
      logical obrasph
      character*256 str_bra
      integer int_ket_sph, ket_type, l2_ket, iket, l_ket
      logical oketsph
      character*256 str_ket
c
      write(6,*)' integrals of type <ishell| Operator | jshell>'
      write(6,*)' enter a 1 for spherical bra basis; else 0 '
      read(5,*)int_bra_sph
      write(6,*)' enter a 1 for spherical ket basis; else 0 '
      read(5,*)int_ket_sph
      write(6,*)' enter the type of the bra function (ishell)'
      write(6,*)' 0=s, 1=p, 2=d, etc.'
      read(5,*)bra_type
      write(6,*)' enter the type of the ket function (ishell)'
      write(6,*)' 0=s, 1=p, 2=d, etc.'
      read(5,*)ket_type
c
      if (int_bra_sph.eq.1) then
        obrasph = .true.
        l2_bra = 2*bra_type + 1
      else
        obrasph = .false.
        l2_bra = (bra_type+1)*(bra_type+2)/2
      endif
      if (int_ket_sph.eq.1) then
        oketsph = .true.
        l2_ket = 2*ket_type + 1
      else
        oketsph = .false.
        l2_ket = (ket_type+1)*(ket_type+2)/2
      endif
      do ibra = 1,l2_bra
        call int_str_ang(bra_type,ibra,obrasph,str_bra,l_bra)
        do iket = 1,l2_ket
          call int_str_ang(ket_type,iket,oketsph,str_ket,l_ket)
          write(6,*)'<',str_bra(1:l_bra),'|',str_ket(1:l_ket),'>'
        enddo
      enddo
      end
      subroutine int_case4_order()
      implicit none
c
      integer int_bra_sph, bra_type, l2_bra, ibra, l_bra
      logical obrasph
      character*256 str_bra
      integer int_ket_sph, ket_type, l2_ket, iket, l_ket
      logical oketsph
      character*256 str_ket
      integer lval_max, dip_val
      character*256 str_dip
      integer l2_dip, idip, l_dip
      
c
      write(6,*)' integrals of type <ishell| dipole | jshell>'
      write(6,*)' enter a 1 for spherical bra basis; else 0 '
      read(5,*)int_bra_sph
      write(6,*)' enter a 1 for spherical ket basis; else 0 '
      read(5,*)int_ket_sph
      write(6,*)' enter the type of the bra function (ishell)'
      write(6,*)' 0=s, 1=p, 2=d, etc.'
      read(5,*)bra_type
      write(6,*)' enter the type of the ket function (ishell)'
      write(6,*)' 0=s, 1=p, 2=d, etc.'
      read(5,*)ket_type
      write(6,*)' the dipole code generates the dipole integrals for '
      write(6,*)' ALL lval of 0 up to lval_max'
      write(6,*)' enter the lval_max for the dipole '
      read (5,*)lval_max 
c
      if (int_bra_sph.eq.1) then
        obrasph = .true.
        l2_bra = 2*bra_type + 1
      else
        obrasph = .false.
        l2_bra = (bra_type+1)*(bra_type+2)/2
      endif
      if (int_ket_sph.eq.1) then
        oketsph = .true.
        l2_ket = 2*ket_type + 1
      else
        oketsph = .false.
        l2_ket = (ket_type+1)*(ket_type+2)/2
      endif
* do the zero case
      do ibra = 1,l2_bra
        call int_str_ang(bra_type,ibra,obrasph,str_bra,l_bra)
        do iket = 1,l2_ket
          call int_str_ang(ket_type,iket,oketsph,str_ket,l_ket)
          write(6,*)'<',str_bra(1:l_bra),'| |',str_ket(1:l_ket),'>'
        enddo
      enddo
      do dip_val = 1,lval_max
        l2_dip = (dip_val+1)*(dip_val+2)/2
        do ibra = 1,l2_bra
          call int_str_ang(bra_type,ibra,obrasph,str_bra,l_bra)
          do idip = 1,l2_dip
            call int_str_xyz(dip_val,idip,str_dip,l_dip)
            do iket = 1,l2_ket
              call int_str_ang(ket_type,iket,oketsph,str_ket,l_ket)
              write(6,*)'<',str_bra(1:l_bra),'|',str_dip(1:l_dip),
     &            '|',str_ket(1:l_ket),'>'
            enddo
          enddo
        enddo
      enddo
      end
