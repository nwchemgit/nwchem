	subroutine printbasis
	implicit none
#include "cguess.fh"
        integer isym, ibas, lobas, iprim
c
*       write(6,*)
*       write(6,*) ' zn=', zn, 'nbct=', nbct, ' nsym=', nsym, ' nsht=', nsht
C$Id$
        lobas =  1
	do isym = 1, nsym
           write(6,*) ' isym ', isym
           write(6,*) '       nbas  ', nbas(isym)
           write(6,*) '       ncsh  ', ncsh(isym)
           write(6,*) '       nosh  ', nosh(isym)
           write(6,*) '       nccup ', nccup(isym)
           write(6,*) '       n1    ', n1(isym)
           write(6,*) '       nbc   ', nbc(isym)
           do ibas = lobas, lobas+nbc(isym)-1
              write(6,*) '       function ', ibas
              do iprim = nstrt(ibas), nstrt(ibas+1)-1
                 write(6,1) iprim-nstrt(ibas)+1,zeta(iprim), cont(iprim)
 1               format('            ', i5,2f20.6)
              enddo
           enddo
           lobas = lobas + nbc(isym)
        enddo
        write(6,2) ajmn
 2      format(1x,4d16.8)
        write(6,*)
        write(6,*)
        end


