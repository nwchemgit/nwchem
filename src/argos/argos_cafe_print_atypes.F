      subroutine argos_cafe_print_atypes(lfnout,nam,num,wgt,vdw)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfnout
      character*6 nam(mat,3)
      integer num(mat,3)
      real*8 wgt(mat,6),vdw(mat,mat,map,6)
c
      integer i,j,k,l
c
      write(lfnout,1000)
 1000 format(/,' ATOM TYPES',//,t9,'Set',t13,'Type',
     + t24,'Atomic',t32,'Mass',/,t24,'number')
c
      do 1 i=1,mat
      write(lfnout,1001) i,(j,nam(i,j),num(i,j),wgt(i,j),j=1,nparms)
 1001 format(/,2i5,2x,a6,4x,i5,f12.6,/,(5x,i5,2x,a6,4x,i5,f12.6))
    1 continue
c
      write(lfnout,1002)
 1002 format(/,' VAN DER WAALS PARAMETERS',//,
     + t14,'Set',t19,'Atom types',t30,4(' v. d. Waals'),/,
     + t30,2('  dispersion'),2('   repulsion'),/,
     + t30,2(12x,'  3rd neighb'))
c
      do 2 i=1,mat
      do 3 j=i,mat
      write(lfnout,1003) i,j,(k,nam(i,k)(1:5),nam(j,k)(1:5),
     + (vdw(i,j,l,k),l=1,4),k=1,nparms)
 1003 format(/,3i5,3x,a5,1x,a5,4(1pe12.5),/,
     + (10x,i5,3x,a5,1x,a5,4(1pe12.5)))
    3 continue
    2 continue
c
      return
      end
c $Id$
