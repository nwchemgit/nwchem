      logical function tce_energy_fragment(rtdb)
c
c $Id$
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
      integer nfrags, ifrag, jfrag, kfrag
      double precision energy
      integer k_fcount, l_fcount
      integer k_fdipole, l_fdipole
      integer k_fdipole_2, l_fdipole_2
      integer k_efrags, l_efrags
      integer k_efragpairs, l_efragpairs
      integer k_efragtriplets, l_efragtriplets
      logical dipole
      logical scf_only
      logical dft_only
      logical mp2_only
      logical binding
      integer iter
      double precision r
      integer maxiter
      parameter (maxiter = 40)
      double precision residual
      double precision max_residual
      double precision damping
      double precision thresh
      parameter (thresh = 1.0d-4)
      double precision radius
      logical nodezero
      logical tce_energy
      external tce_energy
      logical scf
      external scf
      logical dft_energy
      external dft_energy
      logical semi_dir_mp2
      external semi_dir_mp2

#ifdef MRCC_METHODS
ckbn mrcc-r-1 -3
      logical tce_mrcc_energy
      external tce_mrcc_energy
      integer imrcc
#endif

c
      nodezero = (ga_nodeid().eq.0)
      if (.not.rtdb_get(rtdb,'fragment:dipole',mt_log,1,dipole))
     1  dipole = .true.
      if (.not.rtdb_get(rtdb,'fragment:scf_only',mt_log,1,scf_only))
     1  scf_only = .false.
      if (.not.rtdb_get(rtdb,'fragment:dft_only',mt_log,1,dft_only))
     1  dft_only = .false.
      if (.not.rtdb_get(rtdb,'fragment:mp2_only',mt_log,1,mp2_only))
     1  mp2_only = .false.
      if (.not.rtdb_get(rtdb,'fragment:binding',mt_log,1,binding))
     1  binding = .false.
      if (.not.rtdb_get(rtdb,'fragment:radius',mt_dbl,1,radius))
     1  radius = 1.0d99
      if (.not.rtdb_get(rtdb,'fragment:damping',mt_dbl,1,damping))
     1  damping = 1.0d0

#ifdef MRCC_METHODS
ckbn mrcc-r-1 -3
      if (.not.rtdb_get(rtdb,'tce:mrcc',mt_int,1,imrcc))
     1  imrcc = 0
#endif

      if (.not.rtdb_get(rtdb,'tce:fragment',mt_int,1,fragment))
     1  call errquit('tce_energy_fragment',0,RTDB_ERR)
      if (fragment.eq.-1) then
#ifdef MRCC_METHODS
ckbn mrcc-r-1 -3
        if (imrcc.eq.0) then
#endif
          tce_energy_fragment = tce_energy(rtdb,.true.)
#ifdef MRCC_METHODS
ckbn mrcc-r-1 -3
        else
          tce_energy_fragment = tce_mrcc_energy(rtdb,.true.)
        endif
#endif
        return
        write(*,*) "I am here2"
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
          if (.not.ma_push_get(mt_int,nfrags,'fragment counts',
     1      l_fcount,k_fcount)) call errquit
     2      ('tce_energy_fragment: failed to allocate fcount',0,MA_ERR)
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
            int_mb(k_fcount+ifrag-1) = 0
            dbl_mb(k_fdipole+(ifrag-1)*3)   = 0.0d0
            dbl_mb(k_fdipole+(ifrag-1)*3+1) = 0.0d0
            dbl_mb(k_fdipole+(ifrag-1)*3+2) = 0.0d0
          enddo
          do iter = 1, maxiter
            do ifrag = 1, nfrags
              dbl_mb(k_fdipole_2+(ifrag-1)*3)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3)
              dbl_mb(k_fdipole_2+(ifrag-1)*3+1)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3+1)
              dbl_mb(k_fdipole_2+(ifrag-1)*3+2)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3+2)
            enddo
            do ifrag = 1, nfrags
              call select_fragment(rtdb,ifrag,dbl_mb(k_fdipole))
              if (.not.dft_only) then
              if (.not.scf(rtdb)) call errquit
     1          ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
              if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,energy))
     1          call errquit('tce_energy_fragment: scf energy',
     2          0,RTDB_ERR)
              else
              if (.not.dft_energy(rtdb)) call errquit
     1          ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
              if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,energy))
     1          call errquit('tce_energy_fragment: dft energy',
     2          0,RTDB_ERR)
              endif
              dbl_mb(k_efrags+ifrag-1) = energy
              if (nodezero.and.(scf_only.or.dft_only)) write(LuOut,9010)
     1          ifrag,dbl_mb(k_efrags+ifrag-1)
              if (dipole) then
                if (.not.dft_only) then
                if (.not.rtdb_get(rtdb,'scf:dipole',mt_dbl,3,
     1            dbl_mb(k_fdipole+(ifrag-1)*3)))
     2            call errquit('tce_energy_fragment: scf dipole',
     3            0,RTDB_ERR)
                else
                if (.not.rtdb_get(rtdb,'dft:dipole',mt_dbl,3,
     1            dbl_mb(k_fdipole+(ifrag-1)*3)))
     2            call errquit('tce_energy_fragment: dft dipole',
     3            0,RTDB_ERR)
                endif
              endif
            enddo
            max_residual = 0.0d0
            do ifrag = 1, nfrags
              residual=         (dbl_mb(k_fdipole  +(ifrag-1)*3)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3))**2
              residual=residual+(dbl_mb(k_fdipole  +(ifrag-1)*3+1)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3+1))**2
              residual=residual+(dbl_mb(k_fdipole  +(ifrag-1)*3+2)
     1                          -dbl_mb(k_fdipole_2+(ifrag-1)*3+2))**2
              residual = dsqrt(residual)
              if (nodezero) write(LuOut,9040) ifrag,
     1                                dbl_mb(k_fdipole+(ifrag-1)*3)
     2                               ,dbl_mb(k_fdipole+(ifrag-1)*3+1)
     3                               ,dbl_mb(k_fdipole+(ifrag-1)*3+2)
     4                               ,residual
              if (residual.gt.max_residual) max_residual = residual
            enddo
            if (nodezero) write(LuOut,9050) max_residual
            if (max_residual.lt.thresh) goto 10
            do ifrag = 1, nfrags
              dbl_mb(k_fdipole  +(ifrag-1)*3)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3)*damping
     2        + dbl_mb(k_fdipole_2+(ifrag-1)*3)*(1.0d0-damping)
              dbl_mb(k_fdipole  +(ifrag-1)*3+1)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3+1)*damping
     2        + dbl_mb(k_fdipole_2+(ifrag-1)*3+1)*(1.0d0-damping)
              dbl_mb(k_fdipole  +(ifrag-1)*3+2)
     1        = dbl_mb(k_fdipole  +(ifrag-1)*3+2)*damping
     2        + dbl_mb(k_fdipole_2+(ifrag-1)*3+2)*(1.0d0-damping)
            enddo
          enddo
          call errquit('tce_energy_fragment: failed to converge',
     1                 0,CALC_ERR)
c
c         TCE monomer fragments
c
   10     continue
          if (.not.(scf_only.or.dft_only)) then
            do ifrag = 1, nfrags
              call select_fragment(rtdb,ifrag,dbl_mb(k_fdipole))
              if (fragment.eq.0) then
                if (.not.mp2_only) then
                if (.not.tce_energy(rtdb,.false.)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,energy))
     1            call errquit('tce_energy_fragment: tce energy',
     2            0,RTDB_ERR)
                else
                if(.not.rtdb_cput(rtdb,'sdmp2:task',1,'energy'))
     1            call errquit('tce_energy_fragment',0, RTDB_ERR)
                if (.not.semi_dir_mp2(rtdb)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'mp2:energy',mt_dbl,1,energy))
     1            call errquit('tce_energy_fragment: mp2 energy',
     2            0,RTDB_ERR)
                endif
                dbl_mb(k_efrags+ifrag-1) = energy
                if (nodezero) write(LuOut,9010) ifrag,
     1            dbl_mb(k_efrags+ifrag-1)
              else if (fragment.eq.ifrag) then
                if (.not.mp2_only) then
                if (.not.tce_energy(rtdb,.true.)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                else
                if(.not.rtdb_cput(rtdb,'sdmp2:task',1,'energy'))
     1            call errquit('tce_energy_fragment',0, RTDB_ERR)
                if (.not.semi_dir_mp2(rtdb)) call errquit
     1            ('tce_energy_fragment: fragment',0,UNKNOWN_ERR)
                endif
                if (.not.binding) then
                  if (mp2_only) call errquit
     1              ('tce_energy_fragment: mp2_only',1,UNKNOWN_ERR)
                  if (.not.rtdb_get(rtdb,'tce:omega',mt_dbl,1,energy))
     1              call errquit('tce_energy_fragment: tce omega',
     2              1,RTDB_ERR)
                  dbl_mb(k_efrags+ifrag-1) = energy
                  if (nodezero) write(LuOut,9110) ifrag,
     1              dbl_mb(k_efrags+ifrag-1)
                else
                  if (.not.mp2_only) then
                  if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,energy))
     1              call errquit('tce_energy_fragment: tce energy',
     2              1,RTDB_ERR)
                  else
                  if (.not.rtdb_get(rtdb,'mp2:energy',mt_dbl,1,energy))
     1              call errquit('tce_energy_fragment: mp2 energy',
     2              1,RTDB_ERR)
                  endif
                  dbl_mb(k_efrags+ifrag-1) = energy
                  if (nodezero) write(LuOut,9010) ifrag,
     1              dbl_mb(k_efrags+ifrag-1)
                endif
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
     1          dbl_mb(k_fdipole),r)
              if (r.le.radius) then
                if (fragment.eq.0) then
                  int_mb(k_fcount+ifrag-1) = int_mb(k_fcount+ifrag-1)+1
                  int_mb(k_fcount+jfrag-1) = int_mb(k_fcount+jfrag-1)+1
                  if (.not.(scf_only.or.dft_only)) then
                    if (.not.mp2_only) then
                    if (.not.tce_energy(rtdb,.false.)) call errquit
     1                ('tce_energy_fragment: fragment pair',0,
     2                UNKNOWN_ERR)
                    if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,
     1                energy))
     2                call errquit('tce_energy_fragment: tce energy',
     3                3,RTDB_ERR)
                    else
                    if(.not.rtdb_cput(rtdb,'sdmp2:task',1,'energy'))
     1                call errquit('tce_energy_fragment',0, RTDB_ERR)
                    if (.not.semi_dir_mp2(rtdb)) call errquit
     1                ('tce_energy_fragment: fragment pair',0,
     2                UNKNOWN_ERR)
                    if (.not.rtdb_get(rtdb,'mp2:energy',mt_dbl,1,
     1                energy))
     2                call errquit('tce_energy_fragment: mp2 energy',
     3                3,RTDB_ERR)
                    endif
                  else
                    if (.not.dft_only) then
                    if (.not.scf(rtdb)) call errquit
     1                ('tce_energy_fragment: scf',0,
     2                UNKNOWN_ERR)
                    if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,
     1                energy))
     2                call errquit('tce_energy_fragment: scf energy',
     3                0,RTDB_ERR)
                    else
                    if (.not.dft_energy(rtdb)) call errquit
     1                ('tce_energy_fragment: dft',0,
     2                UNKNOWN_ERR)
                    if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,
     1                energy))
     2                call errquit('tce_energy_fragment: dft energy',
     3                0,RTDB_ERR)
                    endif
                  endif
                  dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1) = energy
                  if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1            dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                else if ((fragment.eq.ifrag).or.
     1                   (fragment.eq.jfrag)) then
                  int_mb(k_fcount+ifrag-1) = int_mb(k_fcount+ifrag-1)+1
                  int_mb(k_fcount+jfrag-1) = int_mb(k_fcount+jfrag-1)+1
                  if (mp2_only) call errquit
     1              ('tce_energy_fragment: mp2_only',1,UNKNOWN_ERR)
                  if (.not.binding) then
                    if (.not.tce_energy(rtdb,.true.)) call errquit
     1                ('tce_energy_fragment: tce',0,
     2                UNKNOWN_ERR)
                    if (.not.rtdb_get(rtdb,'tce:omega',mt_dbl,1,energy))
     1                call errquit('tce_energy_fragment: tce omega',
     2                0,RTDB_ERR)
                    dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
     1                = energy
                    if (nodezero) write(LuOut,9130) ifrag,jfrag,
     1              dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                  else
                    if (.not.(scf_only.or.dft_only)) then
                      if (.not.mp2_only) then
                      if (.not.tce_energy(rtdb,.false.)) call errquit
     1                  ('tce_energy_fragment: tce',0,
     2                  UNKNOWN_ERR)
                      if (.not.rtdb_get(rtdb,'tce:energy',mt_dbl,1,
     1                  energy))
     2                  call errquit('tce_energy_fragment: tce energy',
     3                  1,RTDB_ERR)
                      else
                      if(.not.rtdb_cput(rtdb,'sdmp2:task',1,'energy'))
     1                  call errquit('tce_energy_fragment',0, RTDB_ERR)
                      if (.not.semi_dir_mp2(rtdb,.false.)) call errquit
     1                  ('tce_energy_fragment: mp2',0,
     2                  UNKNOWN_ERR)
                      if (.not.rtdb_get(rtdb,'mp2:energy',mt_dbl,1,
     1                  energy))
     2                  call errquit('tce_energy_fragment: mp2 energy',
     3                  1,RTDB_ERR)
                      endif
                    else
                      if (.not.dft_only) then
                      if (.not.scf(rtdb)) call errquit
     1                  ('tce_energy_fragment: scf',0,
     2                  UNKNOWN_ERR)
                      if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,
     1                  energy))
     2                  call errquit('tce_energy_fragment: scf energy',
     3                  0,RTDB_ERR)
                      else
                      if (.not.dft_energy(rtdb)) call errquit
     1                  ('tce_energy_fragment: dft',0,
     2                  UNKNOWN_ERR)
                      if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,
     1                  energy))
     2                  call errquit('tce_energy_fragment: dft energy',
     3                  0,RTDB_ERR)
                      endif
                    endif
                    dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
     1                = energy
                    if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1              dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                  endif
                endif
              else
                dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1) = 0.0d0
                if ((fragment.eq.0).and.nodezero)
     1            write(LuOut,9035) ifrag,jfrag
              endif
            enddo
          enddo
c
c         Trimer fragments
c
          if (.not.ma_push_get(mt_dbl,nfrags**3,'fragment triplet',
     1      l_efragtriplets,k_efragtriplets)) call errquit
     2      ('tce_energy_fragment: failed to allocate efrags',0,MA_ERR)
          do ifrag = 1, nfrags-2
            do jfrag = ifrag+1, nfrags-1
              do kfrag = jfrag+1, nfrags
                call select_fragmenttriplet(rtdb,ifrag,jfrag,kfrag,
     1            dbl_mb(k_fdipole))
                if (.not.dft_energy(rtdb)) call errquit
     1            ('tce_energy_fragment: dft',0,
     2            UNKNOWN_ERR)
                if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,
     1            energy))
     2            call errquit('tce_energy_fragment: dft energy',
     3            0,RTDB_ERR)
                dbl_mb(k_efragtriplets+(kfrag-1)*nfrags*nfrags+
     1            (ifrag-1)*nfrags+jfrag-1) = energy
                if (nodezero) write(LuOut,9330) ifrag,jfrag,kfrag,
     1            dbl_mb(k_efragtriplets+(kfrag-1)*nfrags*nfrags+
     2            (ifrag-1)*nfrags+jfrag-1)
              enddo
            enddo
          enddo
          if (nodezero) write(LuOut,*)
          if (nodezero) write(LuOut,9000)
          if (nodezero) write(LuOut,*)
          if (fragment.eq.0) then
            energy = 0.0d0
            do ifrag = 1, nfrags-2
              do jfrag = ifrag+1, nfrags-1
                do kfrag = jfrag+1, nfrags
                if (nodezero) write(LuOut,9330) ifrag,jfrag,kfrag,
     1          dbl_mb(k_efragtriplets+(kfrag-1)*nfrags*nfrags+
     2            (ifrag-1)*nfrags+jfrag-1)
                energy = energy
     1          + dbl_mb(k_efragtriplets+(kfrag-1)*nfrags*nfrags+
     2            (ifrag-1)*nfrags+jfrag-1)
     3          - dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
     4          - dbl_mb(k_efragpairs+(ifrag-1)*nfrags+kfrag-1)
     5          - dbl_mb(k_efragpairs+(jfrag-1)*nfrags+kfrag-1)
     6          + dbl_mb(k_efrags+ifrag-1)
     7          + dbl_mb(k_efrags+jfrag-1)
     8          + dbl_mb(k_efrags+kfrag-1)
                enddo
              enddo
            enddo
            do ifrag = 1, nfrags-1
              do jfrag = ifrag+1, nfrags
                if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1          dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                energy = energy
     1                 + dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
     2          - dbl_mb(k_efrags+ifrag-1)
     3          - dbl_mb(k_efrags+jfrag-1)
              enddo
            enddo
            do ifrag = 1, nfrags
              if (nodezero) write(LuOut,9010) ifrag,
     1          dbl_mb(k_efrags+ifrag-1),int_mb(k_fcount+ifrag-1)
              energy = energy
     1               + dbl_mb(k_efrags+ifrag-1)
            enddo
            if (nodezero) write(LuOut,9020) energy
          else
            energy = 0.0d0
            do ifrag = 1, nfrags-1
              do jfrag = ifrag+1, nfrags
                if ((fragment.eq.ifrag).or.(fragment.eq.jfrag)) then
                  if (.not.binding) then
                    if (nodezero) write(LuOut,9130) ifrag,jfrag,
     1              dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                  else
                    if (nodezero) write(LuOut,9030) ifrag,jfrag,
     1              dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                  endif
                  energy = energy
     1            +dbl_mb(k_efragpairs+(ifrag-1)*nfrags+jfrag-1)
                endif
              enddo
            enddo
            do ifrag = 1, nfrags
              if (.not.binding) then
                if (fragment.eq.ifrag) then
                  if (nodezero) write(LuOut,9110) ifrag,
     1              dbl_mb(k_efrags+ifrag-1),int_mb(k_fcount+ifrag-1)
                  energy = energy
     1                   - dble(max(0,int_mb(k_fcount+ifrag-1)-1))
     2                   * dbl_mb(k_efrags+ifrag-1)
                endif
              else
                if (nodezero) write(LuOut,9010) ifrag,
     1            dbl_mb(k_efrags+ifrag-1),int_mb(k_fcount+ifrag-1)
                energy = energy
     1               - dble(int_mb(k_fcount+ifrag-1))
     2               * dbl_mb(k_efrags+ifrag-1)
              endif
            enddo
            if (.not.binding) then
              if (nodezero) write(LuOut,9120) energy
            else
              if (nodezero) write(LuOut,9025) energy
            endif
          endif
          if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,energy))
     1      call errquit('tce_energy_fragment: RTDB problem',0,RTDB_ERR)
          if (.not.ma_pop_stack(l_efragtriplets))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_efragpairs))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_efrags))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_fdipole_2))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_fdipole))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          if (.not.ma_pop_stack(l_fcount))
     1      call errquit("tce_energy_fragment: MA problem",1,MA_ERR)
          tce_energy_fragment = .true.
          return
        endif
      endif
 9000 format(1x,'@ Fragment energies')
 9010 format(1x,'@         Fragment',i3,' energy = ',F25.15,' hartree '
     1         ,i4,' times')
 9020 format(1x,'@               Total energy = ',F25.15,
     1          ' hartree')
 9025 format(1x,'@             Binding energy = ',F25.15,
     1          ' hartree')
 9030 format(1x,'@ Fragment pair',i3,i3,' energy = ',F25.15,' hartree')
 9035 format(1x,'@ Fragment pair',i3,i3,' skipped')
 9040 format(1x,'@ Fragment',i3,' dipole = ',3F10.6,
     1          ' a.u. residual = ',F20.15)
 9330 format(1x,'@ Fragment triplet',3i3,' energy = ',F25.15,' hartree')
 9050 format(1x,'@ Dipole residual = ',3F20.15)
 9110 format(1x,'@         Fragment',i3,' excitation energy = ',
     1  F15.5,' eV ',i4,' times')
 9120 format(1x,'@    Total excitation energy = ',F15.5,
     1          ' eV')
 9130 format(1x,'@ Fragment pair',i3,i3,' excitation energy = ',
     1  F15.5,' eV')
 9140 format(1x,'@ Fragment',i3,' dipole = ',3F10.6,' a.u.')
      end
