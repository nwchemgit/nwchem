C $Id$
      subroutine drdy_ratoms(string,istrt)
      implicit none ! double precision (a-h,o-z)
#include "errquit.fh"
#include "drdyP.fh"
#include "stdio.fh"
#include "nwc_drdyP.fh"
c
*::functions:
      Integer  drdy_icint
      External drdy_icint
      Character*2 drdy_atsym
      External drdy_atsym
      Double Precision  drdy_cfloat
      External drdy_cfloat
*::passed:
      Integer istrt
      Character*80 string
*::local:
      Integer i,iatm,ib,ie,ierr,itot,itry,j,label(natoms),nfind
      Integer isect,iend
      Double Precision ams(103), dummy
c
c     array of atomic masses
c
      data (ams(i),i=1,54)  /
     *   1.007825d+00,4.0026d+00,7.01600d+00,9.01218d+00,11.00931d+00,
     *   12.0d+00,14.00307d+00,15.99491d+00,18.99840d+00,19.99244d+00,
     *   22.9898d+00,23.98504d+00,26.98153d+00,27.97693d+00,
     *   30.97376d+00,31.97207d+00,34.96885d+00,39.948d+00,
     *   38.96371d+00,39.96259d+00,44.95592d+00,47.90d+00,50.9440d+00,
     *   51.9405d+00,54.9381d+00,55.9349d+00,58.9332d+00,57.9353d+00,
     *   62.9298d+00,63.9291d+00,68.9257d+00,73.9219d+00,74.9216d+00,
     *   79.9165d+00,78.9183d+00,83.9115d+00,
     *   84.9117d+00,87.9056d+00,89.9054d+00,89.9043d+00,92.9060d+00,
     *   97.9055d+00,97.0d+00,101.9037d+00,102.9048d+00,105.9032d+00,
     *   106.9041d+00,113.9036d+00,114.9041d+00,119.9022d+00,
     *   120.9038d+00,129.9067d+00,126.9044d+00,131.9042d+00/
      data (ams(i),i=55,103)  /
     *   132.9054d+00,137.9052d+00,138.9063d+00,139.9054d+00,
     *   140.9076d+00,141.9077d+00,144.9127d+00,151.9197d+00,
     *   152.9212d+00,157.9241d+00,158.9253d+00,163.9292d+00,
     *   164.9303d+00,165.9303d+00,168.9342d+00,173.9389d+00,
     *   174.9408d+00,179.9465d+00,180.9480d+00,183.9509d+00,
     *   186.9557d+00,191.9615d+00,192.9629d+00,194.9648d+00,
     *   196.9665d+00,201.9706d+00,
     *   204.9744d+00,207.9766d+00,208.9804d+00,208.9824d+00,
     *   209.9871d+00,222.0176d+00,
     *   223.0197d+00,226.0254d+00,
     *   227.0278d+00,232.0381d+00,231.0359d+00,238.0508d+00,
     *   237.0482d+00,244.0642d+00,243.0614d+00,247.0703d+00,
     *   247.0703d+00,251.0796d+00,252.0829d+00,257.0751d+00,
     *   258.0986d+00,259.1009d+00,260.1053d+00/
c
c     read in each line in this set until an 'end' is found.  For
c     each line check if the atomic name or number is given.  If
c     an associated mass is not give assign the atom the mass of
c     the most common isotope.
c
      itot = 1
      call drdy_rline(string,istrt,isect,iend)
      do while (string(istrt:istrt+2).ne.'end' .and. itot.lt.natoms)
C Check for input error
         if (isect.eq.1.or.iend.eq.1) then
            write(fu6,1000)
            call errquit('Needs an END',555, INPUT_ERR)
         end if
c
c first read in atom number
        ib = istrt
        call drdy_rword(string,istrt,ierr)
        ie = istrt - 1
        iatm = drdy_icint(string(ib:ie))
c See if the first non-blank character is numeric; if so
c   assume atomic number was given.
         j = ichar(string(istrt:istrt)) - ichar('0')
         if (j.ge.0.and.j.le.9) then
             label(iatm) = drdy_icint(string(istrt:80))
c
c              check that atomic number is less than or equal to 103
             if (label(iatm).le.0.or.label(iatm).gt.103) then
                write(fu6,1100)label(iatm),string
             end if
         else
c
c           for non-numeric assume it is the atomic symbol
c
             nfind = 0
             i = 1
             do while(nfind.eq.0.and.i.le.103)
               if (string(istrt:istrt+1).eq.drdy_atsym(i)) then
                   label(iatm) = i
                   nfind = 1
               end if
               i = i+1
             end do
c
             if (nfind.eq.0) then 
                write(fu6,1200)string(istrt:istrt+1)
                call errquit('ratoms: fatal error',911, INPUT_ERR)
             end if
         end if
c
c        find atomic masses -  search for mass after end of label
c
         call drdy_rword(string,istrt,ierr)
         if (ierr.eq.1) then
             xmass(iatm) = ams(label(iatm))
             nwcmass(iatm) = xmass(iatm)
         else
             dummy = drdy_cfloat(string(istrt:80))
             xmass(iatm) = dummy
             nwcmass(iatm) = xmass(iatm)
         end if
         chat(iatm) = drdy_atsym(label(iatm))
C Capitalize first letter in atomic label
         itry = ichar(chat(iatm)(1:1))
         chat(iatm)(1:1) = char(itry-32)
c
c store atomic charge (number) and symbol for nwchem
c
         atomic_charge(iatm) = dble(label(iatm))
         atomic_labels(iatm) = chat(iatm)
         atomic_atomct(iatm) = 0.0d0
c        atomic_labels(iatm) = asymb(label(iatm))
c
c        read in next atom
c
         itot = itot+1
         call drdy_rline(string,istrt,isect,iend)
      end do
c
c     set number of atoms variable
c
      natom = itot - 1
      if (natom.gt.natoms) then
        write(fu6,*)' number of atoms too large '
        write(fu6,*)' natom (number read) = ',natom
        write(fu6,*)' natoms (dimension)  = ',natoms
        call errquit('drdy_ratoms: fatal error',911, INPUT_ERR)
      endif
c
      return
1000  format(3x,'ERROR:  the atoms keyword must be terminated with',
     *  ' END')
1100  format(3x,'Incorrect atomic number:  ',i2,
     * ' determined from the following line in file 5:'
     * /1x,a80)
1200  format(3x,'ERROR: illegal atomic symbol: ',A2)
      end
