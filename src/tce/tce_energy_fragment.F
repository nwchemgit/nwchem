      logical function tce_energy_fragment(rtdb)
c
c $Id: tce_energy_fragment.F,v 1.2 2004-03-15 20:43:31 sohirata Exp $
c
c Wrapper for tce_energy for fragment MO calculations.
c
c Written by So Hirata, Jan 2004. (c) Battelle, PNNL, 2004.
c
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
c
      integer rtdb             ! Run-time database
      integer fragment         ! -1 if not fragment MO
      integer nfrags, ifrag, jfrag
      double precision energy
      integer k_fdipole, l_fdipole
      integer k_fdipole_2, l_fdipole_2
      integer k_efrags, l_efrags
      integer k_efragpairs, l_efragpairs
      logical dipole
      logical scf_only
      integer iter
      integer maxiter
      parameter (maxiter = 10)
      double precision residual
      double precision thresh
      parameter (thresh = 1.0d-4)
      logical nodezero
      logical tce_energy
      external tce_energy
      logical scf
      external scf
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.rtdb_get(rtdb,'fragment:dipole',mt_log,1,dipole))
     1  dipole = .true.
      if (.not.rtdb_get(rtdb,'fragment:scf_only',mt_log,1,scf_only))
     1  scf_only = .false.
      if (.not.rtdb_get(rtdb,'tce:fragment',mt_int,1,fragment))
     1  call errquit('tce_energy_fragment',0,RTDB_ERR)
      if (fragment.eq.-1) then
        tce_energy_fragment = tce_energy(rtdb,.true.)
        return
      else
        tce_energy_fragment = .false.
        call geom_fragment(rtdb,nfrags)
        if (nfrags.eq.1) then
          tce_energy_fragment = tce_energy(rtdb,.true.)
          return
        else
c
c         SCF monomer fragments
c
          if (.not.ma_push_get(mt_dbl,3*nfrags,'fragment dipole',
     1      l_fdipole,k_fdipole)) call errquit
     2      ('tce_energy_fragment: failed to allocate fdipole',0,MA_ERR)
          if (.not.ma_push_get(mt_dbl,3*nfrags,'fragment dipole',
     1      l_fdipole_2,k_fdipole_2)) call errquit
     2      ('tce_energy_fragment: failed to allocate fdipole',0,MA_ERR)
          if (.not.ma_push_get(mt_dbl,nfrags,'fragment',
     1      l_efrags,k_efrags)) call errquit
     2      ('tce_energy_fragment: failed to allocate efrags',0,MA_ERR)
          do ifrag = 1, nfrags
            dbl_mb(k_fdipole+(ifrag-1)*3)   = 0.0d0
            dbl_mb(k_fdipole+(ifrag-1)*3+1) = 0.0d0
            dbl_mb(k_fdipole+(ifrag-1)*3+2) = 0.0d0
          enddo
          do iter = 1, maxiter
            do ifrag = 1, nfrags
              call select_fragment(rtdb,ifrag,dbl_mb(k_fdipole))
              if (.not.scf(rtdb)) call errquit
     1          ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
              if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,energy))
     1          call errquit('tce_energy_fragment',0,RTDB_ERR)
              dbl_mb(k_efrags+ifrag-1) = energy
              if (nodezero.and.scf_only) write(LuOut,9010) ifrag,
     1          dbl_mb(k_efrags+ifrag-1)
              if (dipole) then
                call dcopy(3*nfrags,dbl_mb(k_fdipole),  1,
     1                              dbl_mb(k_fdipole_2),1) 
                if (.not.rtdb_get(rtdb,'scf:dipole',mt_dbl,3,
     1            dbl_mb(k_fdipole+(ifrag-1)*3)))
     2            call errquit('tce_energy_fragment',0,RTDB_ERR)
              endif
            enddo
            residual = 0.0d0
            do ifrag = 1, nfrags
              if (nodezero) write(LuOut,9040) ifrag,
     1                                dbl_mb(k_fdipole+(ifrag-1)*3)
     2                               ,dbl_mb(k_fdipole+(ifrag-1)*3+1)
     3                               ,dbl_mb(k_fdipole+(ifrag-1)*3+2)
              residual=residual+(dbl_mb(k_fdipole  +(ifrag-1)*3)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3))**2
              residual=residual+(dbl_mb(k_fdipole  +(ifrag-1)*3+1)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3+1))**2
              residual=residual+(dbl_mb(k_fdipole  +(ifrag-1)*3+2)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3+2))**2
            enddo
            residual = dsqrt(residual)
            if (nodezero) write(LuOut,9050) residual
            if (residual.lt.thresh) goto 10
          enddo
          call errquit('tce_energy_fragment: failed to converge',
     1                 0,CALC_ERR)
c
c         TCE monomer fragments
c
   10     continue
          if (.not.scf_only) then
            do ifrag = 1, nfrags
              call select_fragment(rtdb,ifrag,dbl_mb(k_fdipole))
              if (fragment.eq.0) then
                if (.not.tce_energy(rtdb,.false.)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,energy))
     1            call errquit('tce_energy_fragment',0,RTDB_ERR)
                dbl_mb(k_efrags+ifrag-1) = energy
                if (nodezero) write(LuOut,9010) ifrag,
     1            dbl_mb(k_efrags+ifrag-1)
              else if (fragment.eq.ifrag) then
                if (.not.tce_energy(rtdb,.true.)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'tce:omega',mt_dbl,1,energy))
     1            call errquit('tce_energy_fragment',0,RTDB_ERR)
                dbl_mb(k_efrags+ifrag-1) = energy
                if (nodezero) write(LuOut,9110) ifrag,
     1            dbl_mb(k_efrags+ifrag-1)
              endif
            enddo
          endif
c
c         Dimer fragments
c
          if (.not.ma_push_get(mt_dbl,nfrags*nfrags,'fragment pair',
     1      l_efragpairs,k_efragpairs)) call errquit
     2      ('tce_energy_fragment: failed to allocate efrags',0,MA_ERR)
          do ifrag = 1, nfrags-1
            do jfrag = ifrag+1, nfrags
              call select_fragmentpair(rtdb,ifrag,jfrag,
     1          dbl_mb(k_fdipole))
              if (.not.scf_only) then
                if (fragment.eq.0) then
                  if (.not.tce_energy(rtdb,.false.)) call errquit
     1              ('tce_energy_fragment: fragment pair',0,UNKNOWN_ERR)
                  if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,energy))
     1              call errquit('tce_energy_fragment',0,RTDB_ERR)
                  dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1) = energy
                  if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1            dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                else if ((fragment.eq.ifrag).or.
     1                   (fragment.eq.jfrag)) then
                  if (.not.tce_energy(rtdb,.true.)) call errquit
     1              ('tce_energy_fragment: fragment pair',0,UNKNOWN_ERR)
                  if (.not.rtdb_get(rtdb,'tce:omega',mt_dbl,1,energy))
     1              call errquit('tce_energy_fragment',0,RTDB_ERR)
                  dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1) = energy
                  if (nodezero) write(LuOut,9130) ifrag,jfrag,
     1            dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                endif
              else
                if (.not.scf(rtdb)) call errquit
     1            ('tce_energy_fragment: fragment pair',0,UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,energy))
     1            call errquit('tce_energy_fragment',0,RTDB_ERR)
                dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1) = energy
                if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1          dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
              endif
            enddo
          enddo
          if (nodezero) write(LuOut,*)
          if (nodezero) write(LuOut,9000)
          if (nodezero) write(LuOut,*)
          if ((fragment.eq.0).or.(scf_only)) then
            energy = 0.0d0
            do ifrag = 1, nfrags-1
              do jfrag = ifrag+1, nfrags
                if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1          dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                energy = energy
     1                 + dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
              enddo
            enddo
            do ifrag = 1, nfrags
              if (nodezero) write(LuOut,9010) ifrag,
     1          dbl_mb(k_efrags+ifrag-1)
              energy = energy-dfloat(nfrags-2)*dbl_mb(k_efrags+ifrag-1)
            enddo
            if (nodezero) write(LuOut,9020) energy
          else
            energy = 0.0d0
            do ifrag = 1, nfrags-1
              do jfrag = ifrag+1, nfrags
                if ((fragment.eq.ifrag).or.(fragment.eq.jfrag)) then
                  if (nodezero) write(LuOut,9130) ifrag,jfrag,
     1            dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                  energy = energy
     1            +dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                endif
              enddo
            enddo
            do ifrag = 1, nfrags
              if (fragment.eq.ifrag) then
                if (nodezero) write(LuOut,9110) ifrag,
     1            dbl_mb(k_efrags+ifrag-1)
                energy = energy-dfloat(nfrags-2)
     1                 * dbl_mb(k_efrags+ifrag-1)
              endif
            enddo
            if (nodezero) write(LuOut,9120) energy
          endif
          if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,energy))
     1      call errquit('tce_energy_fragment: RTDB problem',0,RTDB_ERR)
          if (.not.ma_pop_stack(l_efragpairs))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_efrags))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_fdipole_2))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_fdipole))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          tce_energy_fragment = .true.
          return
        endif
      endif
 9000 format(1x,'@ Fragment energies')
 9010 format(1x,'@         Fragment',i3,' energy = ',F25.15,' hartree')
 9020 format(1x,'@               Total energy = ',F25.15,' hartree')
 9030 format(1x,'@ Fragment pair',i3,i3,' energy = ',F25.15,' hartree')
 9040 format(1x,'@ Fragment',i3,' dipole = ',3F10.6,' a.u.')
 9050 format(1x,'@ Dipole residual = ',3F20.15)
 9110 format(1x,'@         Fragment',i3,' excitation energy = ',
     1  F15.5,' eV')
 9120 format(1x,'@               Total excitation energy = ',
     1  F15.5,' eV')
 9130 format(1x,'@ Fragment pair',i3,i3,' excitation energy = ',
     1  F15.5,' eV')
 9140 format(1x,'@ Fragment',i3,' dipole = ',3F10.6,' a.u.')
      end
