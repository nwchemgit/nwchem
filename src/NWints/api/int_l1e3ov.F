      subroutine int_l1e3ov(i_basis, ish, j_basis, jsh, k_basis, ksh,
     &       zerotol, lov3, OV3, ilab, jlab, klab, numov3, lscr, scr)
c
c routine to calculate labels for 3 center overlap integrals 
c
      implicit none 
#include "basP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c..  dcopy
c..  errquit
c::function
      logical int_chk_init
      external int_chk_init
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer i_basis            ! [input] basis set handle for ish
      integer ish                ! [input] i shell/contraction
      integer j_basis            ! [input] basis set handle for jsh
      integer jsh                ! [input] j shell/contraction
      integer k_basis            ! [input] basis set handle for ksh
      integer ksh                ! [input] k shell/contraction
      double precision zerotol   ! [input] tolerance threshold for keeping integrals
      integer lov3               ! [input] length of 3c overlap integrals array
      double precision OV3(lov3) ! [output] 3c overlap integrals
      integer ilab(lov3)         ! [output] i labels for 3c overlap integrals
      integer jlab(lov3)         ! [output] j labels for 3c overlap integrals
      integer klab(lov3)         ! [output] k labels for 3c overlap integrals
      integer numov3             ! [output] number of integrals and labels generated
      integer lscr               ! [input] length of scratch array
      double precision scr(lscr) ! [input] scratch array
c::local
      integer ibas, jbas, kbas  ! real basis set array index
      integer i, j, k, icount   ! and do loop variables counters
      logical shells_ok         ! all shells valid
      logical nonzero           ! for nonzero check
c
      logical ijbas, jkbas, ikbas, allbas  ! check of basis set handles for canonicalization
c
c
c check initialization and shells
c
      if (.not.int_chk_init('int_l1e3ov'))
     &       call errquit('int_l1e3ov: int_init was not called' ,0)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      shells_ok = shells_ok .and. int_chk_sh(k_basis,ksh)
      if (.not.shells_ok)
     &       call errquit('int_l1e3ov: invalid contraction/shell',0)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      kbas = k_basis + BASIS_HANDLE_OFFSET
      ijbas  = ibas.eq.jbas
      jkbas  = jbas.eq.kbas
      ikbas  = ibas.eq.kbas
      allbas = ijbas .and. jkbas
c
c.. copy vector to scratch array
c
      if (lov3.gt.lscr)
     &       call errquit('int_l1e3ov: scratch array too small',lov3)
      call dcopy(lov3,OV3,1,scr,1)
c
      numov3 = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,ibas), ibs_cn2bfr(2,ish,ibas)
        do 00200 j = ibs_cn2bfr(1,jsh,jbas), ibs_cn2bfr(2,jsh,jbas)
          do 00300 k = ibs_cn2bfr(1,ksh,kbas), ibs_cn2bfr(2,ksh,kbas)
            icount = icount + 1
            nonzero = abs(scr(icount)).ge.zerotol
            if (nonzero) then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .!  keep the integral
              if     (allbas) then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .!  i>=j>=k
                if(i.ge.j.and.j.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (ijbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! i>=j ; free k
                if (i.ge.j) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (jkbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! j>=k ; free i
                if (j.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (ikbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! i>=k ; free j  (least likely!!)
                if (i.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              else
                numov3 = numov3 + 1
                OV3(numov3)  = scr(icount)
                ilab(numov3) = i
                jlab(numov3) = j
                klab(numov3) = k
              endif
            endif
00300     continue
00200   continue
00100 continue
c      
      end
