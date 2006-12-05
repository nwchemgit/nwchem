      subroutine transp_input(rtdb)
*
* $Id: transp_input.F,v 1.3 2006-12-05 13:55:50 edo Exp $
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "transport.fh"

      integer rtdb
c
      integer num_dirs,ind,ipos,plno,atom0,atom1,
     N     nlen,nfield,i,j,
     A     istart,iend
      integer atoi,l_a,k_a,n_input
      external atoi
      character*1 plno_char
      parameter (num_dirs=5)

      character*255 test
      character*4 catom0,catom1
      character*15 dirs(num_dirs)
      double precision vecs(3)
      data dirs /'pl', 'atoms', 'lattice', 'pl_thrs',
     E     'end'/
c

C
 100  if (.not. inp_read())call errquit
     $     ('transp_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test)) 
     $     call errquit('tra_input: failed to read keyword', 0,
     &       INPUT_ERR)
      write(0,*) ' got test ',test(1:inp_strlen(test))
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, 'transport')
c
      goto 100
C player definition
C   PL0 atoms 81-120
C   PL0 lattice 0.0 -4.92 0.0
      else if (inp_contains(.false., 'pl', test,ipos)) then
         plno=atoi(test(3:3))
         write(0,*) ' player no. ',plno
      if (.not. inp_a(test)) 
     $     call errquit('tra_input: needs something after PL?', 0,
     &       INPUT_ERR)
      if (inp_contains(.false.,'atoms', test,ipos)) then
c
c     do a ma_alloc of size natoms
c     we can the copy and reduce to actual size
      if (.not.MA_Push_Get(MT_int,nw_max_atom,'atlist',l_a,k_a))
     &   call errquit('trainp: cannot allocate atlist',0, 0)
c
         nfield = inp_n_field()
         n_input=0
         do i=1,nfield-2
            if (inp_i(atom0) ) then
               write(LuOut,*)'trinp: atoms list found ',atom0
                n_input=n_input+1
                int_mb(k_a+n_input-1)=atom0
            else
               if (.not. inp_a(test)) 
     $              call errquit('tra_input: sthing after atomss', 0,
     &              INPUT_ERR)
               istart=0
               if(.not.inp_strtok(test(1:inp_strlen(test)),
     C              ':', istart, iend)) call errquit(
     A              ' column char not found' ,0,0)
               catom0=test(1:iend)
c     
               catom1=test((iend+2):inp_strlen(test))
               atom0=atoi(catom0(1:iend))
               atom1=atoi(catom1(1:inp_strlen(test)-iend-1))
               write(0,*) ' atoms from ',atom0, ' to ',atom1
               do j=atom0,atom1
                n_input=n_input+1
                int_mb(k_a+n_input-1)=atom0
               enddo
            endif
         enddo
         write(0,*) ' PL no.',plno,' made of ',n_input,' atoms'
         
      elseif(inp_contains(.false.,'lattice', test,ipos)) then
c
c     read the lattice 3 vecs
c
         do i=1,3
            write(0,*) i,plno
            if(.not.inp_f(vecs(i)))call errquit(
     .           ' vec not read ',0,0)
            write(0,*) i,vecs(i)
            latvec(i,plno)=vecs(i)
            write(0,*) latvec(i,plno)
         enddo
      else
         write(0,*) ' keyword ',test(1:inp_strlen(test)),
     C        ' not recognized'
         call errquit(' INPUT error ',0,0)
      endif
      
      goto 100
      else if (inp_compare(.false.,'end', test)) then
         write(0,*) ' got end'
       continue
      endif
      call errquit(' aiaiaia ',0,0)
      return
      end
#if 1
      subroutine player
      call errquit(' stub ',0,0)
      return
      end
#endif
