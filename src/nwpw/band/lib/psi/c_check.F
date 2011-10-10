*     ***********************************
*     *             			*
*     *          c_OrthoCheck		*
*     *             			*
*     ***********************************

      subroutine c_OrthoCheck(indx1,ispin,ne,nbrill,psi1)
      implicit none 
      integer indx1
      integer ispin,ne(2),nbrill
      complex*16 psi1(*)

*    *** local variables ***
      integer npack1

      call Cram_max_npack(npack1)
      call Pneb_orthoCheck(indx1,0,0,npack1,psi1)
      return
      end

c*     ***********************************
c*     *             			*
c*     *       Stfl_g_OrthoCheckMake	*
c*     *             			*
c*     ***********************************
c
c      subroutine Stfl_g_OrthoCheckMake(coutput,
c     >                                 ispin,ne,npack1,neall,nbrill,A)
c      implicit none
c      logical coutput
c      integer ispin,ne(2)
c      integer npack1,neall,nbrill
c      complex*16 A(npack1,neall,nbrill)
c      integer mbb,nbb
c      mbb=0
c      if (ispin.eq.1) mbb=1
c      nbb=0
c      if (nbrill.eq.1) nbb=1
c      call Pneb_OrthoCheckMake(coutput,mbb,nbb,npack1,A)
c      return
c      end
c*     ***********************************
c*     *             			*
c*     *          Stfl_MakeOrtho		*
c*     *             			*
c*     ***********************************
c
c      subroutine Stfl_g_MakeOrtho(nb,nfft3d,ne,psi)
c      implicit none 
c      integer    nb,nfft3d,ne
c      complex*16 psi(nfft3d,ne)
c      integer mbb,nbb
c      mbb=0
c      nbb=0
c      if (nb.eq.0) nbb=1
c      call Pneb_makeOrtho(mbb,nbb,nfft3d,psi)
c      return
c      end

ccccccccc
      character*7 function c_index_name(i)
      integer i

      integer itmp,j0,j1,j2,j3,j4,j5
      character*7 name

      itmp = i

      j5 = itmp/100000
      itmp = itmp - j5*100000
      j4 = itmp/10000
      itmp = itmp - j4*10000
      j3 = itmp/1000
      itmp = itmp - j3*1000
      j2 = itmp/100
      itmp = itmp - j2*100
      j1 = itmp/10
      itmp = itmp - j1*10
      j0 = itmp/1

      name(1:1) = '_'
      name(2:2) = CHAR(j5+ICHAR('0'))
      name(3:3) = CHAR(j4+ICHAR('0'))
      name(4:4) = CHAR(j3+ICHAR('0'))
      name(5:5) = CHAR(j2+ICHAR('0'))
      name(6:6) = CHAR(j1+ICHAR('0'))
      name(7:7) = CHAR(j0+ICHAR('0'))
      c_index_name = name
      return
      end
c $Id$
