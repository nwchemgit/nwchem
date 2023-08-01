      subroutine bse_davidson_guess(pars,diagonal,g_trials,npoles,
     $                              ntrials)
      
      implicit none

#include "bse.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "cdft.fh"

      type(bse_params_t) :: pars
      integer npoles,ntrials,g_trials
      double precision diagonal(npoles)
      double precision amb(npoles)

      character*4 oname,vname
      logical found
      integer amo,imo,itrial,ipole,jpole,indeces(ntrials)
      integer isym,jsym,i,j,pol
      double precision minimum

      call nga_get(pars%g_apb,1,npoles,diagonal,npoles)
      call nga_get(pars%g_amb,1,npoles,amb,npoles)

      if (.not.pars%tda) then
        diagonal(:) = diagonal(:)*amb(:)
      endif
      indeces(:) = -1

      ! Get the indeces of the smallest diagonal elements
      do itrial=1,ntrials
        minimum = 1d6
        do ipole=1,npoles
          if (abs(diagonal(ipole)).lt.minimum) then
            found = .true.
            do jpole=1,itrial
              if (ipole.eq.indeces(jpole)) then
                found = .false.
                exit
              endif
            enddo
            if (found) then
              indeces(itrial) = ipole
              minimum = abs(diagonal(ipole))
            endif
          endif
        enddo
      enddo

      ! Get symmetry of ground state
      isym = 1
      if (pars%ipol.eq.2) then
        do i=1,pars%ipol
          do j=1,pars%nocc(i)
            isym = ieor(isym-1,int_mb(pars%k_irs(i)+j-1)-1)+1
          enddo
        enddo
      endif
      call sym_irrepname(geom, isym, oname)


      ! Start the guess vector
      write(LuOut,*)
      write(LuOut,9000) ntrials
      write(LuOut,9020)
      write(LuOut,9010)
      write(LuOut,9020)
      do itrial=1,ntrials
        ipole = indeces(itrial)
        call ga_put(g_trials,ipole,ipole,itrial,itrial,1d0,1)

        ! Print information
        if (ipole.le.pars%npoles(1)) then
          imo = (ipole-1)/pars%nvir(1) + 1
          amo = ipole - (imo-1)*pars%nvir(1) + pars%nocc(1)
          jsym = ieor(isym-1,int_mb(pars%k_irs(1)+imo-1)-1)+1
          jsym = ieor(jsym-1,int_mb(pars%k_irs(1)+amo-1)-1)+1
          pol = 1
        else
          imo = (ipole - pars%npoles(1) - 1)/pars%nvir(2) + 1
          amo = ipole - pars%npoles(1) - (imo-1)*pars%nvir(2) +
     $          pars%nocc(2)
          jsym = ieor(isym-1,int_mb(pars%k_irs(2)+imo-1)-1)+1
          jsym = ieor(jsym-1,int_mb(pars%k_irs(2)+amo-1)-1)+1
          pol = 2
        endif
        call sym_irrepname(geom, jsym, vname)
        write(luout,9030) itrial,pol,imo,amo,vname,
     $                    dsqrt(diagonal(ipole))*27.2114d0
      enddo

 9000 format(i5,' smallest diagonals (eV) ')
 9010 format('  No. Spin  Occ  Vir  Irrep     Diag')
 9020 format('------------------------------------')
 9030 format(4i5,3x,a4,1x,f8.3)

      end subroutine

