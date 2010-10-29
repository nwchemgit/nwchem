      program compare
* $Id$
      implicit none
#include "nwc_const.fh"
#include "stdio.fh"
c
c     This code compares the basis sets and ecp's in the new basis set library with
c     the current one. It only compares the basis sets that are available in both 
c     libraries. 
c
c     The code assumes that in the running directory there is a file library (which is
c     the current library) and a file library.new (which is the newely generated basis
c     set).
c
c     Code created by Bert de Jong
c
      character*40 name
      integer pmax, cmax
      parameter (pmax=nw_max_prim*nw_max_angular,cmax=nw_max_gen_con)
      double precision ee(pmax),ff(pmax),cc(pmax,cmax),dd(pmax,cmax)
      integer rex1a(pmax),rex1b(pmax)
      double precision rex2a(pmax,2),rex2b(pmax,2)
      logical okay, find
      integer ip, ig, nprim, ng, nbasis, nwrong, ncomp, nnew
c
      open(unit=7,file='library',status='old',
     &    access='sequential',form='formatted',err=90000)
      open(unit=8,file='library.new',status='old',
     &    access='sequential',form='formatted',err=90000)
c
c     First do basis sets
c
      write(*,*) 'Checking basis sets !!!!'
      write(*,*) ' '
      ncomp = 0
      nnew  = 0
      nwrong = 0
      nbasis = 0
911   continue
      find = .false.
      name = " "
c
c     Read basis from new library
c
      call read_basis(nbasis,name,find,8,cc,ee)
      if (.not.find) then
         write(*,'(A,I4,A)') 'Found ',nbasis,' basis sets in library.'
         write(*,'(A,I4,A)') 'Of those ',nnew,' are new ones.'
         write(*,'(A,I4,A,I4,A)') 'Of the ',ncomp,' dupicates ',nwrong,
     $        ' differ going from the old to the new library.'
         write(*,*) ' '
         goto 0611
      endif
      nbasis = nbasis + 1
c
c     Find corresponding basis in old library
c
      call read_basis(1,name,find,7,dd,ff)
      if (.not.find) then 
         nnew = nnew + 1
         goto 911
      endif
c
c     Compare exponents and coefficients
c
      okay = .true.
      call compare_exp(ee,ff,okay)
      call compare_coef(cc,dd,okay)
      if (.not.okay) then
         nwrong = nwrong + 1
         write(*,'(A,A)') "Error in basis set ",name
      else
         write(*,'(I4,A,A)') nbasis,". Successfull comparision of ",name
         ncomp = ncomp + 1
      endif
      goto 911
0611  continue
c
c     Now do ECP's
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'Checking ECP sets !!!!'
      write(*,*) ' '
      ncomp = 0
      nnew  = 0
      nwrong = 0
      nbasis = 0
1011  continue
      find = .false.
      name = " "
c
c     Read ecp from new library
c
      call read_ecp(nbasis,name,find,8,rex1a,rex2a)
      if (.not.find) then
         write(*,'(A,I4,A)') 'Found ',nbasis,' ecp sets in library.'
         write(*,'(A,I4,A)') 'Of those ',nnew,' are new ones.'
         write(*,'(A,I4,A,I4,A)') 'Of the ',ncomp,' dupicates ',nwrong,
     $        ' differ going from the old to the new library.'
         write(*,*) ' '
         goto 1611
      endif
      nbasis = nbasis + 1
c
c     Find corresponding ecp in old library
c
      call read_ecp(1,name,find,7,rex1b,rex2b)
      if (.not.find) then 
         nnew = nnew + 1
         goto 1011
      endif
c
c     Compare power and exponents+coefficients
c
      okay = .true.
      call compare_power(rex1a,rex1b,okay)
      call compare_expcoef(rex2a,rex2b,okay)
      if (.not.okay) then
         nwrong = nwrong + 1
         write(*,'(A,A)') "Error in ecp set ",name
      else
         write(*,'(I4,A,A)') nbasis,". Successfull comparision of ",name
         ncomp = ncomp + 1
      endif
      goto 1011
1611  continue
c
      stop 'Done checking!!'
90000 stop 'Error opening library file(s)'
      end
c
      subroutine read_basis(nbas,basname,dofind,unit,coeff,expon)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "stdio.fh"
      character*5 btype
      character*180 field
      character*40 basname,name
      integer pmax,cmax,ip,ig,nprim,ng,nbasis
      parameter (pmax=nw_max_prim*nw_max_angular,cmax=nw_max_gen_con)
      double precision expon(pmax),coeff(pmax,cmax)
      logical dofind         ! Search for basis in basname if .true.
                             ! else read basis # nbas and put name in basname
                             ! If dofind is true and basname is not found dofind will be set to false
      integer nbas           ! Start at basis # nbas
      integer unit           ! Look for basis on file unit
c
      nbasis = 0
      name = " "
      rewind(unit)
      call inp_init(unit,LuOut)
00001 if (.not.inp_read()) then
c
c     Did not find a basis set, or reached end of file, return
c
         dofind = .false.
         goto 00009
      endif
      call inp_set_field(0) 
      if (.not.inp_a(btype)) goto 80000
      if (.not.inp_a(name)) goto 80000
      if (btype.eq."basis") then
         nbasis = nbasis + 1
         if (((nbasis.gt.nbas).and.(.not.dofind)) .or.  
     $      ((basname.eq.name).and.dofind)) then
            basname = name
c
c     Read exponents and coefficients for all l values until end is found
c
            do ip = 1,pmax
              expon(ip)  = 0.0d00
              do ig = 1,cmax
                coeff(ip,ig) = 0.0d00
              enddo
            enddo
            nprim = 0
c
c     Read "atom l-value" or find "end"
c
            if (.not.inp_read()) goto 80000
00003       call inp_set_field(0)
            if (.not.inp_a(field)) goto 80000
            if (field(1:3).eq.'end') goto 00006
c
c     Read next line: 
c       If exponent + contraction put it in array
c       Else new atom l-value or end -> goto 00003
c
00004       if (.not.inp_read()) goto 80000
            call inp_set_field(0) 
            if (.not.inp_f(expon(nprim+1))) goto 00003
            ng = inp_n_field() - 1
            do ig = 1,ng
              if (.not.inp_f(coeff(nprim+1,ig))) goto 80000
            enddo
            nprim = nprim + 1
            goto 00004
         else
            goto 00008
         endif
      else
         goto 00008
      endif
c
c     Not a basis, or wrong basis -> read until you find "end"
c
00008 if (.not.inp_read()) goto 80000
      if (.not.inp_a(field)) goto 80000
      if (field(1:3).eq.'end') goto 00001
      goto 00008
c
00006 dofind = .true.
00009 return
c
80000 continue
*      call inp_outrec()
      write(LuOut,*)'end of file or error reading',name
      stop ' done with error'
      end
c
      subroutine read_ecp(nbas,basname,dofind,unit,coeff,expon)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "stdio.fh"
      character*5 btype
      character*180 field
      character*40 basname,name
      integer pmax,cmax,ip,ig,nprim,ng,nbasis
      parameter (pmax=nw_max_prim*nw_max_angular,cmax=nw_max_gen_con)
      double precision expon(pmax,2)
      integer coeff(pmax)
      logical dofind         ! Search for basis in basname if .true.
                             ! else read basis # nbas and put name in basname
                             ! If dofind is true and basname is not found dofind will be set to false
      integer nbas           ! Start at basis # nbas
      integer unit           ! Look for basis on file unit
c
      nbasis = 0
      name = " "
      rewind(unit)
      call inp_init(unit,LuOut)
10001 if (.not.inp_read()) then
c
c     Did not find a ecp set, or reached end of file, return
c
         dofind = .false.
         goto 10009
      endif
      call inp_set_field(0) 
      if (.not.inp_a(btype)) goto 80001
      if (.not.inp_a(name)) goto 80001
      if (btype.eq."ecp") then
         nbasis = nbasis + 1
         if (((nbasis.gt.nbas).and.(.not.dofind)) .or.  
     $      ((basname.eq.name).and.dofind)) then
            basname = name
c
c     Read exponents and coefficients for all l values until end is found
c
            do ip = 1,pmax
              expon(ip,1)  = 0.0d00
              expon(ip,2)  = 0.0d00
              coeff(ip) = 0
            enddo
            nprim = 0
c
c           Read line nelec
c
            if (.not.inp_read()) goto 80001
            call inp_set_field(0)
            if (.not.inp_a(field)) goto 80001
            if (field(1:3).eq.'end') goto 10006
c
c     Read "atom l-value" or find "end"
c
            if (.not.inp_read()) goto 80001
10003       call inp_set_field(0)
            if (.not.inp_a(field)) goto 80001
            if (field(1:3).eq.'end') goto 10006
c
c     Read next line: 
c       If power + exponent + contraction put it in array
c       Else new atom l-value or end -> goto 10003
c
10004       if (.not.inp_read()) goto 80001
            call inp_set_field(0) 
            if (.not.inp_i(coeff(nprim+1))) goto 10003
            if (.not.inp_f(expon(nprim+1,1))) goto 80001
            if (.not.inp_f(expon(nprim+1,2))) goto 80001
            nprim = nprim + 1
            goto 10004
         else
            goto 10008
         endif
      else
         goto 10008
      endif
c
c     Not a ecp, or wrong ecp -> read until you find "end"
c
10008 if (.not.inp_read()) goto 80001
      if (.not.inp_a(field)) goto 80001
      if (field(1:3).eq.'end') goto 10001
      goto 10008
c
10006 dofind = .true.
10009 return
c
80001 continue
*      call inp_outrec()
      write(LuOut,*)'end of file or error reading',name
      stop ' done with error'
      end
c
      subroutine compare_expcoef(exp1,exp2,nodiff)
      implicit none
#include "nwc_const.fh"
      integer pmax,i
      parameter (pmax=nw_max_prim*nw_max_angular)
      double precision exp1(pmax,2),exp2(pmax,2)
      logical nodiff
c
      do i = 1, pmax
         if ((exp1(i,1)-exp2(i,1)).gt.1.0d-6) nodiff = .false.
         if ((exp1(i,2)-exp2(i,2)).gt.1.0d-6) nodiff = .false.
      enddo
      return
      end
c  
      subroutine compare_power(exp1,exp2,nodiff)
      implicit none
#include "nwc_const.fh"
      integer pmax,i
      parameter (pmax=nw_max_prim*nw_max_angular)
      integer exp1(pmax),exp2(pmax)
      logical nodiff
c
      do i = 1, pmax
         if ((exp1(i)-exp2(i)).gt.0) nodiff = .false.
      enddo
      return
      end
c  
      subroutine compare_exp(exp1,exp2,nodiff)
      implicit none
#include "nwc_const.fh"
      integer pmax,i
      parameter (pmax=nw_max_prim*nw_max_angular)
      double precision exp1(pmax),exp2(pmax)
      logical nodiff
c
      do i = 1, pmax
         if ((exp1(i)-exp2(i)).gt.1.0d-6) nodiff = .false.
      enddo
      return
      end
c  
      subroutine compare_coef(coef1,coef2,nodiff)
      implicit none
#include "nwc_const.fh"
      integer pmax,cmax,i,j
      parameter (pmax=nw_max_prim*nw_max_angular,cmax=nw_max_gen_con)
      double precision coef1(pmax,cmax),coef2(pmax,cmax)
      logical nodiff
c
      do i = 1, pmax
         do j = 1, cmax
         if (abs(coef1(i,j)-coef2(i,j)).gt.1.0d-6) nodiff = .false.
         enddo
      enddo
      return
      end
