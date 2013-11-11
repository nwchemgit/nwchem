C> \ingroup selci
C> @{
      subroutine selci_anal(irange, roots, ept, ptnorm, enew, ncinew,
     &           ptnorm_mp, ept_mp)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
      dimension irange(21), roots(nroot), ept(nroot), ptnorm(nroot),
     $     enew(nroot)
      dimension ept_mp(nroot), ptnorm_mp(nroot)
      integer selci_isum
c
      write(ifllog,101) (10.0d0**(dble(1-i)*0.5d0), irange(i),
     $     i=1,21)
 101  format(/
     $     ' Analysis of PT Interactions by Size '/
     $     ' ----------------------------------- '//
     $     '      Interaction         Count '/
     $     '     -------------        ----- '/
     $     21 (:5x,f13.10, 5x,i8/))
      write(ifllog,1011) selci_isum(21,irange,1)
 1011 format(5x,13x,4x,'---------'/5x,13x,4x,i9)
c
c     print out final energies
c
      if (iwpt.ne.1) then
         write(ifllog,91) (roots(i), ept(i), roots(i)+ept(i),
     $      dsqrt(ptnorm(i)), i=1,nroot)
      endif
      if (iwpt.ne.0) then
         write(ifllog,92) (roots(i), ept_mp(i), roots(i)+ept_mp(i),
     $      dsqrt(ptnorm_mp(i)), i=1,nroot)
      endif
      if (iwpt.eq.0 .or. iwpt.eq.2) then
         write(ifllog,*)' Selection based on EN pt '
      else
         write(ifllog,*)' Selection based on MP pt '
      endif
      
91    format(/
     $ '       Total Perturbation Corrections for CI energies '/
     $ '       ---------------------------------------------- '//
     $ ' EN|   ci energy      pt correction     ci+pt energy',
     $     '    PT norm '/
     $ ' EN| ---------------  ---------------  ---------------',
     $     '  ---------'/
     $ 20(:' EN:',f15.8,2x,f15.8,2x,f15.8,2x,f9.6/))
92    format(/
     $ '      Total Perturbation Corrections for CI energies '/
     $ '      ---------------------------------------------- '//
     $ ' MP|   ci energy      pt correction     ci+pt energy',
     $     '    PT norm '/
     $ ' MP| ---------------  ---------------  ---------------',
     $     '  ---------'/
     $ 20(:' MP:',f15.8,2x,f15.8,2x,f15.8,2x,f9.6/))
c
      write(ifllog,31) ncinew, noconf
 31   format(/
     $     ' No. of CSF in new ci           ',i8/
     $     ' No. of orbital conf in new ci  ',i8//)
c
      write(ifllog,12) (roots(i), enew(i), roots(i)+enew(i),i=1,nroot)
 12   format(/
     $     '           Predicted New CI Energies '/
     $     '           ------------------------- '//
     $ '   ci energy      pt correction     ci+pt energy'/
     $ ' ---------------  ---------------  ---------------'/
     $ 20(:f15.8,2x,f15.8,2x,f15.8/))
c
      end
C> @}
