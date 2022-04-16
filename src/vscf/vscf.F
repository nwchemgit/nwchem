c $Id$
      subroutine vscfm(NGRID,DMDR,NCOUP,VCFCT,IEXC,
     *                 rtdb,nat,vec,eigen)
      implicit double precision(a-h,o-z)
#include "errquit.fh"
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      character*255 dipole_file
      logical dmdr, dipole_file_exists
      logical linear, restart
      integer geom, rtdb, nat, geom_ref, nnm, geom_new
      dimension vec(3*nat,3*nat), eigen(3*nat)
      integer NSTART(7)
c
c     ----- allocate memory and carry out vscf -----
c
      nc1 = 3*nat
      nc2 = (nc1*nc1+nc1)/2
      nc3 = nc1*nc1
c
c     nnm2 is number of pairs of modes (not triangular storage)
c
      if (.not.rtdb_get(rtdb,'vib:linear',mt_log,1,linear))
     $     linear = .false.
      nnm = 3*nat - 6
      if(linear) nnm = 3*nat-5
      nnm2 = (nnm*nnm-nnm)/2
      nnm3 = nnm*nnm
c
c     when we have only one frequency we can only do ncoup=1
c
      IF (NNM.eq.1) NCOUP=1
c
      ngrid2 = (ngrid*ngrid+ngrid)/2
      ngrid3 = ngrid*ngrid
c
c     number of vibr. states for which wavefunctions are stored
c
      NST=NGRID/2
      IF(NST.LT.4) NST=4
c
c     maximum excitation level used in vibrational MP2
c
      NMAX=NST-1
c
c     number of virtual states in MP2
c
      NVIRST=NMAX*NMAX*NNM2 + NMAX*NNM + 1
C
C     number of triples of normal modes
C     no. of virtual states in MP2 includes triple excitations
c
      if (ncoup.gt.2) then
         NTR=NNM*(NNM-1)*(NNM-2)/6
         NVIRST=NVIRST+NMAX*NMAX*NMAX*NTR
      end if
c
c     memory for vscf
c     the first two of these hold harmonic normal modes and frequencies
c
      if (.not.
     &  ma_push_get(mt_dbl,NNM*NNM*NNM*NGRID*NGRID*NGRID,
     &  'XTRIPV',l_XTRIPV,k_XTRIPV))
     &  call errquit('vscfm: could not allocate l_XTRIPV',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm,'vscf freqs',l_fvscf,k_fvscf))
     &  call errquit('vscfm: could not allocate l_fvscf',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm,'mp2 freqs',l_fmp2,k_fmp2))
     &  call errquit('vscfm: could not allocate l_fmp2',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm,'delta Qs',l_dq,k_dq))
     &  call errquit('vscfm: could not allocate l_dq',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm*ngrid,'grid Qs',l_rq,k_rq))
     &  call errquit('vscfm: could not allocate l_rq',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm*(nnm+1)*ngrid,'awave',l_awave,k_awave))
     &  call errquit('vscfm: could not allocate l_awave',555, MA_ERR)
c
c     Get dipole information from previous frequency calculation
c
      call util_file_name('fd_ddipole',.false., .false.,dipole_file)
      dipole_file_exists = .false.
      inquire(file=dipole_file,exist=dipole_file_exists)
      if (DMDR.and.iexc.eq.1.and.dipole_file_exists) then
         ndipmo=0
         if (.not. 
     &     ma_push_get(mt_dbl,3*nc1,'dipole deriv',l_ddm,k_ddm))
     &     call errquit('vscfm: could not allocate l_ddm',555, MA_ERR)
         open(unit=70,file=dipole_file,form='formatted',status='old',
     &       err=89900,access='sequential')
         do iii = 0,((3*nc1)-1)
           read(70,*,err=89901,end=89902) dbl_tmp
           dbl_mb(k_ddm) = dbl_tmp
         enddo
         close(unit=70,status='keep')
      else 
         ndipmo=1
         if (.not. 
     &     ma_push_get(mt_dbl,nnm*ngrid,'diag dipole x',l_dmx,k_dmx))
     &     call errquit('vscfm: could not allocate l_dmx',555, MA_ERR)
         if (.not. 
     &     ma_push_get(mt_dbl,nnm*ngrid,'diag dipole y',l_dmy,k_dmy))
     &     call errquit('vscfm: could not allocate l_dmy',555, MA_ERR)
         if (.not. 
     &     ma_push_get(mt_dbl,nnm*ngrid,'diag dipole z',l_dmz,k_dmz))
     &     call errquit('vscfm: could not allocate l_dmz',555, MA_ERR)
c
c     Calculation of dipole should be at center of mass
c
         if (.not.rtdb_put(rtdb,'prop:center',mt_int,1,2))
     &       call errquit('vscfm: rtdb_put of dipole:center failed',555,
     &       RTDB_ERR)
      endif 
c
c     Allocate additional memory 
c
      if (.not. 
     &  ma_push_get(mt_dbl,nnm,'scaled freqs',l_sfreq,k_sfreq))
     &  call errquit('vscfm: could not allocate l_sfreq',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm*ngrid,'diagonal V',l_diagv,k_diagv))
     &  call errquit('vscfm: could not allocate l_diagv',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm3*ngrid3,'coupled V',l_coupv,k_coupv))
     &  call errquit('vscfm: could not allocate l_coupv',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm3*ngrid3,'coup dipole x',l_dm2x,k_dm2x))
     &  call errquit('vscfm: could not allocate l_dm2x',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm3*ngrid3,'coup dipole y',l_dm2y,k_dm2y))
     &  call errquit('vscfm: could not allocate l_dm2y',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nnm3*ngrid3,'coup dipole z',l_dm2z,k_dm2z))
     &  call errquit('vscfm: could not allocate l_dm2z',555, MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl,nc1,'delta Xs',l_dx,k_dx))
     &  call errquit('vscfm: could not allocate l_dx',555, MA_ERR)
c
c     Check if we can restart:
c     - check where where we are on diag  when ncoup = 1, 2, 3
c     - check where where we are on coupl when ncoup = 1, 2
c     - check where where we are on tripl when ncoup = 3
c     - get reference geometry from rtdb
c     - if any of the diag, coupl, or tripl is incomplete, deterimine
c       where we were and start from that point, minimizing the overhead
c
c     nstart will be used to house restart data:
c     nstart(1)   = level that was completed, i.e. none=0, diag=1, coupl=2, tripl=3
c     nstart(2-7) = restart point in level, needed are mode and grid point
c                   i.e., two values for diag, four for coupl and six for tripl
c
c     if restart, get reference geometry
c     if not a restart, get the geometry and make a copy to save
c     while we do potential scans
c
      call vscf_restart(rtdb,restart,ncoup,nstart,ldiag,lcoup,ltrip,
     &                     nnm,ngrid,dbl_mb(k_rq),dbl_mb(k_dq),
     &                     dbl_mb(k_diagv),dbl_mb(k_dmx),dbl_mb(k_dmy),
     &                     dbl_mb(k_dmz),dbl_mb(k_dm2x),dbl_mb(k_dm2y),
     &                     dbl_mb(k_dm2z),dbl_mb(k_coupv),
     &                     dbl_mb(k_XTRIPV),ndipmo)
      if (restart) then
        if (.not.geom_create(geom_ref,'vscf_reference')) call errquit
     *    ('vscfm: vscf_reference geom_create failed',551, GEOM_ERR)
        if (.not.geom_rtdb_load(rtdb,geom_ref,'vscf_reference')) 
     *     call errquit('vscfm: geom_rtdb_load failed',551, RTDB_ERR)
        if (.not.geom_rtdb_store(rtdb,geom_ref,'geometry'))
     *     call errquit('vscfm: geom_rtdb_store failed',551, RTDB_ERR)
      else
        if (.not.geom_create(geom_ref,'geometry')) call errquit
     *    ('vscfm: geometry geom_create failed',552, GEOM_ERR)
        if (.not.geom_rtdb_load(rtdb,geom_ref,'geometry')) 
     *    call errquit('vscfm: geometry geom_rtdb_load failed',552,
     &       RTDB_ERR)
        if (.not.geom_rtdb_store(rtdb,geom_ref,'vscf_reference'))
     *    call errquit('vscfm: geom_rtdb_store failed',552, RTDB_ERR)
        if (.not. rtdb_put(rtdb,'vscf:restart',mt_log,1,restart))
     $    call errquit('vscfm: rtdb_put restart',552, RTDB_ERR)
      end if
c
c     pick up a new geometry that is used for the potential scan
c
      if (.not.geom_create(geom_new,'geometry'))
     *  call errquit('vscfm: geom_new geom_create failed',555, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom_new,'geometry'))
     *  call errquit('vscfm: geom_rtdb_load failed',555, RTDB_ERR)
      if (.not.geom_strip_sym(geom_new)) call errquit
     *  ('vscfm: geom_strip_sym failed',555, GEOM_ERR)
C
C     Calculate diagonal and coupling potentials on grids
C
      CALL VGRID(vec,eigen,dbl_mb(k_sfreq),dbl_mb(k_rq),dbl_mb(k_dq),
     *           dbl_mb(k_dx),coords(1,1,geom_ref),dbl_mb(k_diagv),
     *           dbl_mb(k_coupv),dbl_mb(k_XTRIPV),dbl_mb(k_dmx),
     *           dbl_mb(k_dmy),
     *           dbl_mb(k_dmz),dbl_mb(k_dm2x),dbl_mb(k_dm2y),
     *           dbl_mb(k_dm2z),NC1,NAT,NNM,NGRID,NDIPMO,NSTART,
     *           NCOUP,rtdb,geom_new,ldiag,lcoup,ltrip)
      IF (ga_nodeid().eq.0) WRITE(LuOut,9060)
c
c     get rid of some memory that we no longer need
c
      if (.not. ma_pop_stack(l_dx))
     &  call errquit('vscfm:ma_pop of l_dx failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_dm2z))
     &  call errquit('vscfm:ma_pop of l_dm2z failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_dm2y))
     &  call errquit('vscfm:ma_pop of l_dm2y failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_dm2x))
     &  call errquit('vscfm:ma_pop of l_dm2x failed',555, MA_ERR)
c
c     get memory for the next section
c
      if (.not.
     &  ma_push_get(mt_int,nnm,'state indices',l_state,k_state))
     &  call errquit('vscfm: could not allocate l_state',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'temp V',l_tv,k_tv))
     &  call errquit('vscfm: could not allocate l_tv',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm*ngrid,'vscf',l_vscf,k_vscf))
     &  call errquit('vscfm: could not allocate l_vscf',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'vhf',l_vhf,k_vhf))
     &  call errquit('vscfm: could not allocate l_vhf',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm,'enrgy',l_enrgy,k_enrgy))
     &  call errquit('vscfm: could not allocate l_enrgy',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm+1,'ediag',l_ediag,k_ediag))
     &  call errquit('vscfm: could not allocate l_ediag',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm+1,'escf',l_escf,k_escf))
     &  call errquit('vscfm: could not allocate l_escf',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm+1,'emppt',l_emppt,k_emppt))
     &  call errquit('vscfm: could not allocate l_emppt',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm*ngrid,'wave',l_wave,k_wave))
     &  call errquit('vscfm: could not allocate l_wave',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm*ngrid*nst,'vwave',l_vwave,k_vwave))
     &  call errquit('vscfm: could not allocate l_vwave',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm*nst,'virte',l_virte,k_virte))
     &  call errquit('vscfm: could not allocate l_virte',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_int,nnm*nvirst,'vst',l_vst,k_vst))
     &  call errquit('vscfm: could not allocate l_vst',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_int,ngrid,'pvt',l_pvt,k_pvt))
     &  call errquit('vscfm: could not allocate l_pvt',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'xx',l_xx,k_xx))
     &  call errquit('vscfm: could not allocate l_xx',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'a',l_a,k_a))
     &  call errquit('vscfm: could not allocate l_a',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'phi',l_phi,k_phi))
     &  call errquit('vscfm: could not allocate l_phi',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'r',l_r,k_r))
     &  call errquit('vscfm: could not allocate l_r',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'rr',l_rr,k_rr))
     &  call errquit('vscfm: could not allocate l_rr',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'g',l_g,k_g))
     &  call errquit('vscfm: could not allocate l_g',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'v',l_v,k_v))
     &  call errquit('vscfm: could not allocate l_v',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid2,'h',l_h,k_h))
     &  call errquit('vscfm: could not allocate l_h',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'ec',l_ec,k_ec))
     &  call errquit('vscfm: could not allocate l_ec',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid*ngrid,'vecc',l_vecc,k_vecc))
     &  call errquit('vscfm: could not allocate l_vecc',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid*8,'scr1',l_scr1,k_scr1))
     &  call errquit('vscfm: could not allocate l_scr1',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_int,ngrid,'ia1',l_ia1,k_ia1))
     &  call errquit('vscfm: could not allocate l_ia1',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid3,'gr',l_gr,k_gr))
     &  call errquit('vscfm: could not allocate l_gr',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'tmp',l_tmp,k_tmp))
     &  call errquit('vscfm: could not allocate l_tmp',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,ngrid,'twave',l_twave,k_twave))
     &  call errquit('vscfm: could not allocate l_twave',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nvirst,'emp0',l_emp0,k_emp0))
     &  call errquit('vscfm: could not allocate l_emp0',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nvirst,'vmp',l_vmp,k_vmp))
     &  call errquit('vscfm: could not allocate l_vmp',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_dbl,nnm,'ovrlp',l_ovrlp,k_ovrlp))
     &  call errquit('vscfm: could not allocate l_ovrlp',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_int,nvirst*nnm,'virt',l_virt,k_virt))
     &  call errquit('vscfm: could not allocate l_virt',555, MA_ERR)
      if (.not.
     &  ma_push_get(mt_int,nnm,'ref',l_ref,k_ref))
     &  call errquit('vscfm: could not allocate l_ref',555, MA_ERR)
C
C     Perform VSCF and MP2-VSCF
C
      IF (ga_nodeid().eq.0) THEN
         CALL VSCFMP(dbl_mb(k_sfreq),dbl_mb(k_rq),dbl_mb(k_dq),
     *               dbl_mb(k_diagv),dbl_mb(k_coupv),int_mb(k_state),
     *               dbl_mb(k_tv),dbl_mb(k_vscf),dbl_mb(k_vhf),
     *               dbl_mb(k_enrgy),dbl_mb(k_ediag),dbl_mb(k_escf),
     *               dbl_mb(k_emppt),dbl_mb(k_wave),dbl_mb(k_awave),
     *               dbl_mb(k_vwave),dbl_mb(k_virte),int_mb(k_vst),
     *               int_mb(k_pvt),dbl_mb(k_xx),dbl_mb(k_a),
     *               dbl_mb(k_phi),dbl_mb(k_r),dbl_mb(k_rr),
     *               dbl_mb(k_g),dbl_mb(k_v),dbl_mb(k_h),
     *               dbl_mb(k_ec),dbl_mb(k_vecc),dbl_mb(k_scr1),
     *               int_mb(k_ia1),dbl_mb(k_gr),dbl_mb(k_tmp),
     *               dbl_mb(k_twave),dbl_mb(k_emp0),dbl_mb(k_vmp),
     *               dbl_mb(k_ovrlp),int_mb(k_virt),int_mb(k_ref),
     *               dbl_mb(k_fvscf),dbl_mb(k_fmp2),dbl_mb(k_XTRIPV),
     *               vcfct,NNM,NGRID,NGRID2,NST,NVIRST,ncoup,nmax,iexc)
         WRITE(LuOut,9070)
      END IF
c
c     get rid of some more memory that isn't needed
c
      if (.not. ma_pop_stack(l_ref))
     &  call errquit('vscfm:ma_pop of l_ref failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_virt))
     &  call errquit('vscfm:ma_pop of l_virt failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_ovrlp))
     &  call errquit('vscfm:ma_pop of l_ovrlp failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vmp))
     &  call errquit('vscfm:ma_pop of l_vmp failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_emp0))
     &  call errquit('vscfm:ma_pop of l_emp0 failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_twave))
     &  call errquit('vscfm:ma_pop of l_twave failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_tmp))
     &  call errquit('vscfm:ma_pop of l_tmp failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_gr))
     &  call errquit('vscfm:ma_pop of l_gr failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_ia1))
     &  call errquit('vscfm:ma_pop of l_ia1 failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_scr1))
     &  call errquit('vscfm:ma_pop of l_scr1 failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vecc))
     &  call errquit('vscfm:ma_pop of l_vecc failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_ec))
     &  call errquit('vscfm:ma_pop of l_ec failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_h))
     &  call errquit('vscfm:ma_pop of l_h failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_v))
     &  call errquit('vscfm:ma_pop of l_v failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_g))
     &  call errquit('vscfm:ma_pop of l_g failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_rr))
     &  call errquit('vscfm:ma_pop of l_rr failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_r))
     &  call errquit('vscfm:ma_pop of l_r failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_phi))
     &  call errquit('vscfm:ma_pop of l_phi failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_a))
     &  call errquit('vscfm:ma_pop of l_a failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_xx))
     &  call errquit('vscfm:ma_pop of l_xx failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_pvt))
     &  call errquit('vscfm:ma_pop of l_pvt failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vst))
     &  call errquit('vscfm:ma_pop of l_vst failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_virte))
     &  call errquit('vscfm:ma_pop of l_virte failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vwave))
     &  call errquit('vscfm:ma_pop of l_vwave failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_wave))
     &  call errquit('vscfm:ma_pop of l_wave failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_emppt))
     &  call errquit('vscfm:ma_pop of l_emppt failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_escf))
     &  call errquit('vscfm:ma_pop of l_escf failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_ediag))
     &  call errquit('vscfm:ma_pop of l_ediag failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_enrgy))
     &  call errquit('vscfm:ma_pop of l_enrgy failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vhf))
     &  call errquit('vscfm:ma_pop of l_vhf failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_vscf))
     &  call errquit('vscfm:ma_pop of l_vscf failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_tv))
     &  call errquit('vscfm:ma_pop of l_tv failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_state))
     &  call errquit('vscfm:ma_pop of l_state failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_coupv))
     &  call errquit('vscfm:ma_pop of l_coupv failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_diagv))
     &  call errquit('vscfm:ma_pop of l_diagv failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_sfreq))
     &  call errquit('vscfm:ma_pop of l_sfreq failed',555, MA_ERR)
c
c     get memory for the next section
c
      if (.not.
     &  ma_push_get(mt_dbl,nnm,'intensity',l_int,k_int))
     &  call errquit('vscfm: could not allocate l_int',555, MA_ERR)
C
C     Calculate IR intensities
C
      if (ga_nodeid().eq.0) then
      IF (ndipmo.eq.0) THEN    
         if (.not.
     &     ma_push_get(mt_dbl,nnm,'dder',l_dder,k_dder))
     &     call errquit('vscfm: could not allocate l_dder',555, MA_ERR)
         CALL DINTENS(dbl_mb(k_int),dbl_mb(k_ddm),dbl_mb(k_dder),
     *                dbl_mb(k_fvscf),dbl_mb(k_fmp2),vec,dbl_mb(k_rq),
     *                dbl_mb(k_dq),dbl_mb(k_awave),NNM,NGRID,NC1)
         if (.not. ma_pop_stack(l_dder))
     &     call errquit('vscfm:ma_pop of l_dder failed',555, MA_ERR)
      ELSE
         CALL INTENS(dbl_mb(k_int),dbl_mb(k_dmx),dbl_mb(k_dmy),
     *               dbl_mb(k_dmz),dbl_mb(k_fvscf),dbl_mb(k_fmp2),
     *               dbl_mb(k_dq),dbl_mb(k_awave),NNM,NGRID)
      END IF
      endif
c
c clean up memory
c
      if (.not. ma_pop_stack(l_int))
     &  call errquit('vscfm:ma_pop of l_int failed',555, MA_ERR)
      if (ndipmo.eq.0) then
         if (.not. ma_pop_stack(l_ddm))
     &     call errquit('vscfm:ma_pop of l_dmz failed',555, MA_ERR)
      else
         if (.not. ma_pop_stack(l_dmz))
     &     call errquit('vscfm:ma_pop of l_dmz failed',555, MA_ERR)
         if (.not. ma_pop_stack(l_dmy))
     &     call errquit('vscfm:ma_pop of l_dmy failed',555, MA_ERR)
         if (.not. ma_pop_stack(l_dmx))
     &     call errquit('vscfm:ma_pop of l_dmx failed',555, MA_ERR)
      endif
      if (.not. ma_pop_stack(l_awave))
     &  call errquit('vscfm:ma_pop of l_awave failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_rq))
     &  call errquit('vscfm:ma_pop of l_rq failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_dq))
     &  call errquit('vscfm:ma_pop of l_dq failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_fmp2))
     &  call errquit('vscfm:ma_pop of l_fmp2 failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_fvscf))
     &  call errquit('vscfm:ma_pop of l_fvscf failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_XTRIPV))
     &  call errquit('vscfm:ma_pop of l_XTRIPV failed',555, MA_ERR)
c
c     Reset geometry to default (reference) geometry
c
      if (.not.geom_rtdb_load(rtdb,geom_ref,'vscf_reference')) 
     *   call errquit('vscfm: geom_rtdb_load failed',221, RTDB_ERR)
      if (.not.geom_rtdb_store(rtdb,geom_ref,'geometry'))
     *   call errquit('vscfm: geom_rtdb_store failed',221, RTDB_ERR)
      if (.not.geom_rtdb_delete(rtdb,'vscf_reference'))
     *   call errquit('vscfm: geom_rtdb_store failed',221, RTDB_ERR)
      if (.not.geom_destroy(geom_new))
     *   call errquit('vscfm: geom_destroy',221, RTDB_ERR)
      if (.not.geom_destroy(geom_ref))
     *   call errquit('vscfm: geom_destroy',221, RTDB_ERR)
c
      IF (ga_nodeid().eq.0) WRITE(LuOut,9080)
      return
c
89900 call errquit('vscfm: cannot open fd_ddipole file',221, DISK_ERR)
89901 call errquit('vscfm: error reading fd_ddipole file',221, DISK_ERR)
89902 call errquit('vscfm: incomplete fd_ddipole file',221, DISK_ERR)
 9060 FORMAT(1X,'......Done with potentials on grids......')
 9070 FORMAT(1X,'......Finished vibrational SCF......')
 9080 FORMAT(1X,'......Finished calculating IR intensities......')
      end
c*module vscf    *deck vgrid
      SUBROUTINE VGRID(VEC,EIG,FREQ,RQ,DQ,DX,C0,DIAGV,COUPV,TRIPV,
     *                 DM1X,DM1Y,DM1Z,DM2X,DM2Y,DM2Z,
     *                 NC1,NAT1,NNM,NGRID,NDIPMO,
     *                 NFIN,NCOUP,rtdb,geom,
     *                 ldiag,lcoup,ltrip)
C
      implicit none
#include "errquit.fh"
c
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
C
      integer NC1,NAT1,NNM,NGRID,NDIPMO,NCOUP
      LOGICAL psi0,status
      integer ldiag, lcoup, ltrip ! these are unit numbers
      character*32 theory
      character*255 vectors_in
c
      logical task_energy, property
      external task_energy, property
C
      integer rtdb, geom
c
      double precision EIG(NC1),VEC(NC1,NC1),DX(NC1),DIAGV(NNM,NGRID)
      double precision FREQ(NNM),C0(3,NAT1),COUPV(NNM,NNM,NGRID,NGRID)
      double precision DQ(NNM),RQ(NNM,NGRID),TRIPV(NNM,NNM,NNM,
     *                                        NGRID,NGRID,NGRID)
      double precision DM1X(NNM,NGRID),DM1Y(NNM,NGRID),DM1Z(NNM,NGRID) 
      double precision DM2X(NNM,NNM,NGRID,NGRID),
     *                 DM2Y(NNM,NNM,NGRID,NGRID),
     *                 DM2Z(NNM,NNM,NGRID,NGRID)
      double precision dmtemp(3), dmx, dmy, dmz, dmx0, dmy0, dmz0
      double precision qmin, qrange, wave, fact, e, e0
      integer          NFIN(7),nstart, n, m, j, i, k
      integer          mn, jm, km, jl, kl, im, il, jk, ii
C
      double precision zero, one, four, two, amu
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00, FOUR=4.0D+00)
      PARAMETER (TWO=2.0D+00, AMU=1.8229D+03)
C
      FACT = one/sqrt(AMU)
      NSTART=NC1-NNM+1
c
c        mp2 runs must execute gradient code to get density matrix
c        unfortunately, parallel runs don't have fast cutout yet.
c
c     if(mplevl.gt.0) then
c        mpprop=1
c        if(nproc.gt.1) then
c           if(ga_nodeid().eq.0) 
c    &        write(LuOut,*) 'need parallel property cutout'
c           call errquit('vgrid: need mp2 nos for mp2 properties',555)
c        end if
c     end if
C      
c     do we need to get starting orbitals at the input geometry?
c
      psi0 = .true.
c
        if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +      call errquit('vgrid: no input for theory?',0, RTDB_ERR)
      if(psi0) then
         if(ga_nodeid().eq.0) write(LuOut,9005) 
         if (.not.task_energy(rtdb))  ! equilibrium geometry energy
     &     call errquit('vgrid: first energy did not converge!',555,
     &       CALC_ERR)
         if (theory.eq.'dft') then
            if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,e))
     &     call errquit('vgrid: failed to get energy from rtdb 1',555,
     &       RTDB_ERR)
         else
            if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,e))
     &     call errquit('vgrid: failed to get energy from rtdb 1',555,
     &       RTDB_ERR)
         endif
c
c     Point prop:input vectors to the scf or dft vectors so it does
c     not have to redo the scf/dft everytime
c
         if (theory.eq.'dft') then
            if (.not. rtdb_cget(rtdb,'dft:output vectors',1,vectors_in))
     &         call errquit('vgrid: rtdb_cget dft:output failed', 100,
     &                       RTDB_ERR)
         else
            if (.not. rtdb_cget(rtdb,'scf:output vectors',1,vectors_in))
     &          call errquit('vgrid: rtdb_cget scf:output failed', 100,
     &                        RTDB_ERR)
         end if
         if (.not. rtdb_cput(rtdb,'prop:vectors',1,vectors_in))
     &       call errquit('vgrid: rtdb_cput prop:input failed', 100,
     &                     RTDB_ERR)
c
         if (.not.rtdb_put(rtdb,'prop:dipole',mt_int,1,0))
     &     call errquit('vgrid: rtdb_put of dipole failed',555,
     &       RTDB_ERR)
         if (.not.property(rtdb))  ! equilibrium dipole moment
     &     call errquit('vgrid: first dipole did not work!',555,
     &       RTDB_ERR)
         if (.not.rtdb_get(rtdb,'prop:dipval',mt_dbl,3,dmtemp))
     &     call errquit('vgrid: rtdb_get of dipval failed',555,
     &                  RTDB_ERR)
         dmx = dmtemp(1)
         dmy = dmtemp(2)
         dmz = dmtemp(3)
      else
         e   = zero
         dmx = zero
         dmy = zero
         dmz = zero
      end if
C
C     save potential energy, dipole moment at equilibrium 
C
      E0=E
      dmx0=dmx
      dmy0=dmy
      dmz0=dmz
C
C     Frequencies in atomic units
C
      WAVE = 1.6021892D00*1.6021892D00*1.0D2
      WAVE = WAVE/(8.0d0*ACOS(0.0D00)*8.85418782D00*0.5291771D00)
      WAVE = WAVE*6.022045D00/(0.5291771D00*0.5291771D00)
      WAVE = SQRT(WAVE)
      WAVE = WAVE/(4.0d0*ACOS(0.0D00)*2.9979245D00)
      WAVE = WAVE*1.0D4
c
c     Bring in frequencies and convert back from cm-1, and then convert to atomic units
c
      do i=NC1,NSTART,-1   
         ii=NC1+1-i
         freq(ii)=(eig(i)/WAVE)/sqrt(AMU)
      end do
c
c     If all diag values were done, goto coupl
c
      if (nfin(1).gt.0) goto 100
C
C     Calculate diagonal V on a grid
C
      if (ga_nodeid().eq.0) write(LuOut,9040) ngrid, nc1-nstart+1 
C
      do i=NC1,NSTART,-1
         im=nc1-i+1
         Qrange=FOUR
         Qrange=Qrange/sqrt(freq(im))
         Qmin=-Qrange
         dQ(im)=(TWO*Qrange)/(ngrid-1)
         do il=1, ngrid
            RQ(im,il)=Qmin+(il-1)*dQ(im)
            jk= 0
            do j=1, nat1
               do k=1, 3
                  jk=jk+1
                  dx(jk) = FACT*VEC(jk,i)*RQ(im,il)
                  coords(k,j,geom) = C0(k,j) + dx(jk)
               end do
            end do
c
c     Skip pieces on restart until restart point is reached
c     At that point restart point reset to zero to indicate we
c     just want to continue to calculate points
c
c     If nfin(2).gt.zero, then we have not reached the restart
c     point and we need to skip one more
c
                  if (im.eq.nfin(2).and.il.eq.nfin(3)) then
                     nfin(2)=0
                     goto  99
                  else if (nfin(2).gt.0) then
                     goto  99
                  endif
c
            if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     *        call errquit('vgrid: geom_rtdb_store failed',555,
     &       RTDB_ERR)
            if(ga_nodeid().eq.0) write(LuOut,9045) il,i
            status = task_energy(rtdb)
            if (theory.eq.'dft') then
            if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,e))
     &      call errquit('vgrid: failed to get energy from rtdb 2',555,
     &       RTDB_ERR)
            else
            if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,e))
     &      call errquit('vgrid: failed to get energy from rtdb 2',555,
     &       RTDB_ERR)
            endif
            IF(.not.status .AND. ga_nodeid().eq.0) WRITE(LuOut,9000)
c
            diagV(im,il) = e - e0
            if (ndipmo.eq.1) then
              if (status) then
              if (.not.property(rtdb))  
     &          call errquit('vgrid: first dipole did not work!',555,
     &       RTDB_ERR)
              if (.not.rtdb_get(rtdb,'prop:dipval',mt_dbl,3,dmtemp))
     &          call errquit('vgrid: rtdb_get of dipval failed',555,
     &       RTDB_ERR)
              dmx = dmtemp(1)
              dmy = dmtemp(2)
              dmz = dmtemp(3)
              else
                dmx = dmx0
                dmy = dmy0
                dmz = dmz0
              end if
              dm1x(im,il) = dmx - dmx0 
              dm1y(im,il) = dmy - dmy0 
              dm1z(im,il) = dmz - dmz0 
              if (ga_nodeid().eq.0) write(ldiag,9015) im, RQ(im,il), 
     *             diagV(im,il),dm1x(im,il),dm1y(im,il),dm1z(im,il)
              if(ga_nodeid().eq.0) call util_flush(ldiag)
            else
              if (ga_nodeid().eq.0) 
     &           write(ldiag,9010) im, RQ(im,il), diagV(im,il)
              if(ga_nodeid().eq.0) call util_flush(ldiag)
            end if
  99        continue
         end do
      end do
      if (ga_nodeid().eq.0) write(LuOut,9050)
C
C     Off-diagonal pair coupling potential
C
  100 CONTINUE
      if ((ncoup.le.1) .or. (nfin(1).gt.1)) goto 200
C
      if (ga_nodeid().eq.0) write(LuOut,9060) ngrid,ngrid,
     *                           (nc1-nstart+1)*(nc1-nstart)/2
C
      do i=NC1, NSTART+1,-1
         im=nc1-i+1
         do j=i-1, NSTART,-1
            jm=nc1-j+1
            do il=1, ngrid
               do jl=1, ngrid
                  mn=0
                  do m=1, nat1
                     do n=1, 3
                        mn=mn+1
                        dx(mn)= FACT*VEC(mn,i)*RQ(im,il) +
     *                          FACT*VEC(mn,j)*RQ(jm,jl)
                        coords(n,m,geom) = C0(n,m) + dx(mn)
                     end do
                  end do
c
c     Skip pieces on restart until restart point is reached
c     At that point restart point reset to zero to indicate we
c     just want to continue to calculate points
c
c     If nfin(2).gt.zero, then we have not reached the restart
c     point and we need to skip one more
c
                  if (nfin(1) .le. 1) then ! Do them all
                  else if (im.eq.nfin(2).and.jm.eq.nfin(4).and.
     &                il.eq.nfin(3).and.jl.eq.nfin(5)) then
                     nfin(2)=0
                     goto 199
                  else if (nfin(2).gt.0) then
                     goto 199
                  endif
c
                  if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     *              call errquit('vgrid: geom_rtdb_store failed',555,
     &       RTDB_ERR)
                  if(ga_nodeid().eq.0) write(LuOut,9065) il,jl,i,j
                  status = task_energy(rtdb)
                  if (theory.eq.'dft') then
                  if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,e))
     *              call errquit
     *                ('vgrid: failed to get energy from rtdb 3',555,
     &       RTDB_ERR)
                  else
                  if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,e))
     *              call errquit
     *                ('vgrid: failed to get energy from rtdb 3',555,
     &       RTDB_ERR)
                  endif
                  IF(.not.status .AND. ga_nodeid().eq.0) 
     *               WRITE(LuOut,9000)
c
                  coupV(im,jm,il,jl)=e-diagV(im,il)-diagV(jm,jl)-e0
                  if (ndipmo.eq.1) then
                   if (status) then
                    if (.not.property(rtdb))  
     &               call errquit('vgrid: dipole did not work!',555,
     &       RTDB_ERR)
c
c  Need to not errquit in above line and write an error message
c
                    if (.not.rtdb_get(rtdb,'prop:dipval',
     &               mt_dbl,3,dmtemp))
     &                call errquit
     &                ('vgrid: rtdb_get of dipval failed',555, RTDB_ERR)
                    dmx = dmtemp(1)
                    dmy = dmtemp(2)
                    dmz = dmtemp(3)
                   else
                    dmx = dmx0
                    dmy = dmy0
                    dmz = dmz0
                   end if
                   dm2x(im,jm,il,jl) = dmx - dmx0
                   dm2y(im,jm,il,jl) = dmy - dmy0
                   dm2z(im,jm,il,jl) = dmz - dmz0
                   if (ga_nodeid().eq.0) then
#if defined(FUJITSU) || defined(SOLARIS) || defined(GCC46)
                      backspace 82
#endif
                    write(lcoup,9022) im, jm, RQ(im,il), RQ(jm,jl),
     *                  coupV(im,jm,il,jl),dm2x(im,jm,il,jl),
     *                  dm2y(im,jm,il,jl),dm2z(im,jm,il,jl)
                      call util_flush(lcoup)
                   endif
                  else
                   if(ga_nodeid().eq.0) 
     &              write(lcoup,9020) im,jm,RQ(im,il),RQ(jm,jl),
     *                             coupV(im,jm,il,jl)
                   if(ga_nodeid().eq.0) call util_flush(lcoup)
                  end if
  199             continue
               end do
            end do
         end do
      end do
c
      e=e0
      if (ga_nodeid().eq.0) write(LuOut,9070) 
C
C     3-body coupling potential
C
  200 CONTINUE
      if ((ncoup.le.2) .or. (nfin(1).eq.3)) goto 300
C
      if (ga_nodeid().eq.0) then
         write(LuOut,*)
         write(LuOut,*) ' 3-body coupling potential'
         write(LuOut,*)
      end if
C
      if (ga_nodeid().eq.0) write(LuOut,9080) 
 
      do i=NC1, NSTART+2,-1
        im=nc1-i+1
        do j=i-1, NSTART+1,-1
          jm=nc1-j+1
          do k=j-1, NSTART,-1
            km=nc1-k+1
            do il=1, ngrid
              do jl=1, ngrid
                do kl=1, ngrid
                  mn=0
                  do m=1, nat1
                    do n=1, 3
                      mn=mn+1
                      dx(mn)= FACT*VEC(mn,i)*RQ(im,il) +
     *                        FACT*VEC(mn,j)*RQ(jm,jl) +
     *                        FACT*VEC(mn,k)*RQ(km,kl)
                      coords(n,m,geom) = C0(n,m) + dx(mn)
                    end do
                  end do
c
c     Skip pieces on restart until restart point is reached
c     At that point restart point reset to zero to indicate we
c     just want to continue to calculate points
c
c     If nfin(2).gt.zero, then we have not reached the restart
c     point and we need to skip one more
c
                  if (nfin(1) .le. 2) then ! Do them all
                  else if (im.eq.nfin(2).and.jm.eq.nfin(4).and.
     &                il.eq.nfin(3).and.jl.eq.nfin(5).and.
     &                km.eq.nfin(6).and.kl.eq.nfin(7)) then
                     nfin(2)=0
                     goto 299
                  else if (nfin(2).gt.0) then
                     goto 299
                  endif
c
                  if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     *              call errquit('vgrid: geom_rtdb_store failed',555,
     *                           RTDB_ERR)
                  if(ga_nodeid().eq.0) 
     *              write(LuOut,9085) il,jl,kl,i,j,k
                  status = task_energy(rtdb)
                  if (theory.eq.'dft') then
                  if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,e))
     *              call errquit
     *                ('vgrid: failed to get energy from rtdb',555,
     *                 RTDB_ERR)
                  else
                  if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,e))
     *              call errquit
     *                ('vgrid: failed to get energy from rtdb',555,
     *                 RTDB_ERR)
                  endif
                  IF(.not.status .AND. ga_nodeid().eq.0) 
     *              WRITE(LuOut,9000)
 
                  tripV(im,jm,km,il,jl,kl)=e-e0-
     *                  diagV(im,il)-diagV(jm,jl)-diagV(km,kl)-
     *                  coupV(im,jm,il,jl)-coupV(im,km,il,kl)-
     *                  coupV(jm,km,jl,kl)
                  if(ga_nodeid().eq.0) write(LuOut,9025) im,jm,km,
     *                            RQ(im,il),RQ(jm,jl),RQ(km,kl),
     *                            tripV(im,jm,km,il,jl,kl)
                  if(ga_nodeid().eq.0) call util_flush(LuOut)
  299             continue
                end do
              end do
            end do
          end do
        end do
      end do
 
      e=e0
      if (ga_nodeid().eq.0) write(LuOut,9090) 
c
  300 CONTINUE
      if (.not. rtdb_delete(rtdb,'prop:vectors'))
     &    call errquit('vgrid: rtdb_delete failed',100,RTDB_ERR)
c
      RETURN
C
 9000 FORMAT(1x,4(2h*-),'*'/1X,'warning !'/1x,4(2h*-),'*'/
     *       1x,'SCF HAS NOT CONVERGED at this vscf point!')
 9005 format(//1x,'evaluating wavefunction at original geometry...'//)
 9010 format(1x,I2,2x,f14.8,2x,e14.7)
 9015 format(1x,I2,2x,f14.8,2x,e14.7,3(2x,e14.7))
 9020 format(1x,2(I2,2x),2(f14.8,2x),e14.7) 
 9022 format(1x,2(I2,2x),2(f14.8,2x),e14.7,3(e14.7,2x))
 9025 format(1x,3(I2,2x),3(f14.8,2x),e14.7) 
 9030 format(1x,f14.9)
 9040 FORMAT(//1X,'Starting diagonal potential on a grid of ',I3,
     *       ' points for ',I3,' normal modes') 
 9045 format(//1x,'VSCF: Energy and dipole for grid point',i3,
     *          ' along mode',i3)
 9050 FORMAT(/1X,'Done with diagonal potential')
 9060 FORMAT(//1X,'Starting pair coupling potential on a square grid'/
     *  ' of',I3,' by',I3,' points for',I6,' pairs of normal modes') 
 9065 format(//1x,'VSCF: Energy and dipole for pair grid points',2i3,
     *          ' for mode pair',2i3)
 9070 FORMAT(/1X,'Done with pair coupling potential')
 9080 FORMAT(//1X,'Starting 3-body coupling potential')
 9085 format(//1x,'VSCF: Energy and dipole for 3-body grid points',3i3,
     *          ' for mode triplet',3i3)
 9090 FORMAT(/1X,'Done with 3-body coupling potential')
      END
C
      SUBROUTINE VSCFMP(FREQ,RQ,DQ,DIAGV,COUPV,ISTATE,TV,VSCF,VHF,
     *            E,EDIAG,ESCF,EMPPT,WAVE,ALLWAVE,VIRTWAVE,VIRTE,
     *            IVST,IPVT,XX,A,PHI,R,RR,G,V,H,EC,VEC,SCR,IA,GR,
     *            TEMP,TWAVE,EMP0,VMP,OVRLP,JVIRT,JREF,frscf,frmp2,
     *            tripv,vcfct,NNM,NGRID,NGRID2,NST,NVST,ncoup,nmax,
     *            iexc)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
c
#include "global.fh"
#include "stdio.fh"
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION DIAGV(NNM,NGRID),COUPV(NNM,NNM,NGRID,NGRID)
      DIMENSION TRIPV(NNM,NNM,NNM,NGRID,NGRID,NGRID)
      DIMENSION FREQ(NNM),DQ(NNM),RQ(NNM,NGRID) 
      DIMENSION tV(NGRID),vscf(NNM,NGRID),VHF(NGRID),E(NNM) 
      dimension ediag(NNM+1), escf(NNM+1), emppt(NNM+1) 
      dimension wave(NNM,NGRID),allwave(NNM+1,NNM,NGRID)
      dimension virtwave(NNM,NST,NGRID),virtE(NNM,NST)
      dimension ISTATE(NNM),IVST(NVST,NNM),ipvt(NGRID)
      dimension xx(NGRID),A(NGRID),phi(NGRID,NGRID)
      dimension R(NGRID,NGRID),RR(NGRID,NGRID),G(NGRID,NGRID)
      dimension V(NGRID,NGRID),H(NGRID2),EC(NGRID),VEC(NGRID,NGRID)
      dimension GR(NGRID,NGRID),temp(NGRID),SCR(NGRID,8),IA(NGRID)
      dimension twave(NGRID),emp0(NVST),vmp(NVST)
      dimension jvirt(NVST,NNM),jref(NNM),ovrlp(NNM)
      dimension frscf(nnm),frmp2(nnm) 
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (EPS=1.0D-06, cm=2.19474D+05)
      PARAMETER (mxiter=100)
      logical print_wv
c
c     Print vibrational wave functions
c
      print_wv=.false.
C
C     SCALE PAIR COUPLING POTENTIAL IF NEEDED 
C
      vcfct1=one
  500 fct=vcfct/vcfct1 
      if (ncoup.gt.1) then
      do i=1, NNM-1
         do j=i+1, NNM
            do l=1, NGRID
               do m=1, NGRID
C
                  if ( (abs(coupV(i,j,l,m)).ge.diagV(i,l)).or.
     *                 (abs(coupV(i,j,l,m)).ge.diagV(j,m)) ) 
     *               coupV(i,j,l,m) = fct*coupV(i,j,l,m)
C
                  coupV(j,i,m,l) = coupV(i,j,l,m)
C
               end do
            end do
         end do
      end do
      end if
C
C     SCALE 3-BODY COUPLING POTENTIAL IF NEEDED
C
      if (ncoup.gt.2) then
      do i=1, NNM-2
         do j=i+1, NNM-1
            do k=j+1, NNM
               do l1=1, NGRID
               do l2=1, NGRID
               do l3=1, NGRID
C
                 if ( (abs(tripV(i,j,k,l1,l2,l3)).ge.diagV(i,l1)).or.
     *                (abs(tripV(i,j,k,l1,l2,l3)).ge.diagV(j,l2)).or.
     *                (abs(tripV(i,j,k,l1,l2,l3)).ge.diagV(k,l3)) )
     *            tripV(i,j,k,l1,l2,l3) = fct*tripV(i,j,k,l1,l2,l3)
C
                  tripV(i,k,j,l1,l3,l2) = tripV(i,j,k,l1,l2,l3)
                  tripV(j,i,k,l2,l1,l3) = tripV(i,j,k,l1,l2,l3)
                  tripV(j,k,i,l2,l3,l1) = tripV(i,j,k,l1,l2,l3)
                  tripV(k,j,i,l3,l2,l1) = tripV(i,j,k,l1,l2,l3)
                  tripV(k,i,j,l3,l1,l2) = tripV(i,j,k,l1,l2,l3)
C
               end do
               end do
               end do
            end do
         end do
      end do
      end if
C
C     Start loop over states
C
      do index=0, NNM
c
         if (index.eq.0) then
            write(LuOut,9019)
         else
            write(LuOut,9029) nnm
         endif
         write(LuOut,*) ''
C      
C        array of state indices
C
         do i=1, NNM
            if (index.eq.0) then
               istate(i)=0
            else if (index.eq.i) then
               istate(i)=IEXC
            else
               istate(i)=0
            end if
         end do
C
c        INITIALIZE THE WAVEFUNCTIONS
C
         CALL dfill(nnm*ngrid,zero,wave,1)
         CALL dfill(nnm*ngrid,zero,vscf,1)
C
c        CALCULATE THE EIGENVALUES AND EIGENFUNCTIONS OF
c        DIAGONAL POTENTIAL 
C
         sumE=zero
         do mode=1, NNM
            do l=1, NGRID
               tV(l)=DIAGV(mode,l)
            end do
C
            call collocat(mode,tV,tE,wave,allwave,virtwave,virtE,
     *           ipvt,rq,dq,xx,a,phi,r,rr,g,v,h,ec,vec,scr,ia,
     *           gr,temp,twave,nnm,ngrid,ngrid2,nst,istate,index)
C
            E(mode)=tE
            sumE=sumE+tE
         end do
         ediag(index+1)=sumE
         Etot=sumE
         emp1=zero
         if (ncoup.le.1) goto 50
C
C        START SCF ITERATIONS
C
         write(LuOut,9050)
c
         iter=0
 1000    continue 
         iter=iter+1
         Eprev=Etot
         if (iter.eq.1) Eprev=zero
         Etot=zero
         sumE=zero
         do mode=1, NNM
C
C           effective potential
C
            call Veffect (mode,vhf,wave,dq,COUPV,TRIPV,NNM,NGRID,ncoup)
C
            do l=1, NGRID
               tV(l)=vhf(l)+diagV(mode,l)
               Vscf(mode,l)=vhf(l)
            end do
C
            call collocat(mode,tV,tE,wave,allwave,virtwave,virtE,
     *           ipvt,rq,dq,xx,a,phi,r,rr,g,v,h,ec,vec,scr,ia,gr,
     *           temp,twave,nnm,ngrid,ngrid2,nst,istate,index)
C
            if (tE.lt.zero) then 
               vcfct1=fct
               vcfct=vcfct-1.0D-02
               if (vcfct.lt.2.0D-01) then
                  IF(ga_nodeid().eq.0) 
     *            WRITE(LuOut,*) 'SCALING FACTOR IS LESS THAN 0.2'
                  CALL errquit('vscfmp: problem with potential',555,
     &       INPUT_ERR)
               end if
               if (ga_nodeid().eq.0) write(LuOut,9010) vcfct
               goto 500 
            endif
            E(mode)=tE
            sumE=sumE + tE
         end do
C
C        calculate the SCF correction
C
         call scfcorr (wave,emp1,dq,coupV,TRIPV,NNM,NGRID,ncoup)
C
         Etot=sumE-emp1
C
         write(LuOut,9051) iter,sumE*cm,emp1*cm,Etot*cm
c
         if (ABS(Etot-Eprev).gt.eps .and. iter.le.mxiter) goto 1000 
C
C        exit SCF iterations
C  
         if (ga_nodeid().eq.0) then
            if (iter.le.mxiter) then 
               write(LuOut,*) "VSCF converged in",iter,"iterations"
            else
               write(LuOut,*) "*** VSCF DID NOT CONVERGE ***"
            end if
         end if
   50    continue
c
c        punch VSCF wavefunctions 
c
         if (ga_nodeid().eq.0.and.print_wv) then 
            if (index.eq.0) then
               write(LuOut,9020)
               write(LuOut,9040) vcfct
            else
               write(LuOut,9030) index
               write(LuOut,9040) vcfct
            end if
            write(LuOut,*)
         end if
         sumE=ZERO
         do mode=1, NNM
            if (ga_nodeid().eq.0.and.print_wv) 
     &         write(LuOut,*) mode,istate(mode),E(mode)*cm
            do l=1, NGRID
               if (ga_nodeid().eq.0.and.print_wv) 
     &           write(LuOut,*) mode,rq(mode,l),wave(mode,l)
            end do
            sumE=sumE+E(mode)
            if (ga_nodeid().eq.0.and.print_wv) write(LuOut,*) 
         end do
c
         Etot=sumE-emp1
         escf(index+1)=Etot
         if (ncoup.le.1) goto 100
C
C        MP2 correction 
C
         call virtstate(index,nnm,nmax,nvst,ivst,ncoup,iexc)
C
         call mppt(emp2,emp0,virtwave,vmp,virtE,coupV,tripV,Vscf,dq,
     *           ovrlp,nnm,ngrid,nst,nvst,ivst,istate,jref,jvirt,ncoup)
C
         Etot=Etot+emp2
         emppt(index+1)=Etot
         if (ga_nodeid().eq.0) then
            write(LuOut,*) 
            write(LuOut,*) "E(DIAG) :",ediag(index+1)*cm,
     &                     " cm-1 (without mode coupling)"
            write(LuOut,*) "E(VSCF) :", (Etot-emp2)*cm, 'cm-1'
            if (ncoup.gt.1) write(LuOut,*) "E(VMP2) :", Etot*cm, 'cm-1'
            write(LuOut,*) 
         endif
  100    continue
C
      end do
C
c     calculate vibrational frequencies in cm-1 and print them out
c
      Nstate=NNM+1
      call endiff(freq,frscf,frmp2,ediag,escf,emppt,nnm,Nstate,
     *                                              ncoup,iexc)
C 
      RETURN
C
 9010 FORMAT(/1X,'Scaling coupling potential by ',f4.2) 
 9019 FORMAT(/1X,'Solving VSCF for the ground vibrational state') 
 9020 FORMAT(/1X,'Wavefunctions for the ground vibrational state') 
 9029 FORMAT(/1X,'Solving VSCF for the excited state of mode',I3) 
 9030 FORMAT(/1X,'Wavefunctions for the excited state of mode',I3) 
 9040 FORMAT(1X,'(Scaling factor for coupling potential ',f4.2,')') 
 9050 FORMAT(/1X,'Iteration    E(SCF)      E(MP1)    E(TOTAL)')
 9051 FORMAT(4x,i3,5x,f8.2,4x,f8.2,4x,f8.2)
      end
C
c*module vscf    *deck collocat
      SUBROUTINE collocat(mode,tV,tE,wave,allwave,virtwave,virtE,
     *                    ipvt,x,dx,xx,a,phi,r,rr,g,v,h,e,vec,scr,ia,
     *                    gr,temp,twave,nm,n,n2,nst,istate,index)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
c
#include "global.fh"
#include "stdio.fh"
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION ISTATE(NM)
      DIMENSION x(NM,N),dx(NM)
      dimension wave(NM,N),allwave(NM+1,NM,N)
      dimension tV(N),virtwave(NM,NST,N),virtE(NM,NST)
      dimension ipvt(N)
      dimension xx(N),A(N),phi(N,N)
      dimension R(N,N),RR(N,N),G(N,N),V(N,N)
      dimension H(N2),E(N),VEC(N,N)
      DIMENSION SCR(N,8),IA(N)
      dimension GR(N,N),temp(N),det(2)
      dimension twave(N)
C
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00, FOUR=4.0D+00)
      PARAMETER (TWO=2.0D+00, HALF=0.5D+00, QTR=0.25D+00)
      PARAMETER (c=0.7D+00)
C
c     THIS ROUTINE USES THE COLLOCATIONS METHOD TO GENERATE
C     THE EIGENVALUES AND WAVEFUNCTIONS:  E, wave
c     REF. CPL V153, 1988 pg.98   Yang & Peet
c
      pi=four*atan(ONE)
c
c     NORMAL COORDINATE ON A GRID
c
      do l=1, N
         xx(l)=x(mode,l)
      end do
c
c     GENERATE THE PARAMETERS A'S
c
      A(1)=(c**2)/((xx(2)-xx(1))**2)
      do i=2, N-1
         A(i)=(4*(c**2))/((xx(i+1)-xx(i-1))**2)
      end do
      A(N)=(c**2)/((xx(N)-xx(N-1))**2)
c
c     GENERATE N GAUSSIAN WAVEFUNCTIONS R(i,j)
c
      do i=1, N
         fac=(TWO*A(i)/pi)**QTR
         do j=1, N
            R(i,j)=fac*exp(-A(i)*(xx(j)-xx(i))**2)
            if (abs(R(i,j)).lt.1.0d-99) R(i,j)=zero
            RR(i,j)=R(i,j)
         end do
      end do
c
c     GENERATE THE POTENTIAL MATRIX V(i,j)
c
      do i=1, N
         do j=1, N
            V(i,j)=zero
            if (j.eq.i) V(i,j)=tV(i)
         end do
      end do
c
c     GENERATE THE 2ND ORDER DERIVATIVE OF WAVEFUNCTIONS
c
      do i=1, N
         do j=1, N
            G(i,j)=zero
            dR2=(4*(A(i)**2)*((xx(j)-xx(i))**2))-2*A(i)
            G(i,j)=(-HALF)*dR2*R(i,j)
            if (abs(G(i,j)).lt.1.0d-99) G(i,j)=zero
         end do
      end do
c
c     GENERATE THE INVERSE MATRIX OF R(I,J):  R-1
c
      INFO=0
      CALL DGEFA(R,N,N,IPVT,INFO)
C
      IF(INFO.NE.0) THEN
         IF(ga_nodeid().eq.0) WRITE(LuOut,*) 'MATRIX R IS SINGULAR'
         CALL errquit('collocat: problem generating inverse',555,
     &       CALC_ERR)
      END IF
C
c     options for dgedi routine
c       job = 11 both determinant and inverse
c       job = 01 inverse only
c       job = 10 determinant only
c
      job=01
      call dgedi(R,N,N,ipvt,det,temp,job)
c
c     Multiply kinetic energy matrix G and inverse R
c 
      CALL MRARBR(G,N,N,N,R,N,N,GR,N)

c     GENERATE THE HAMILTONIAN TO BE DIAGONALIZED : GR-1 + V
c     (average the off-diagonal terms and put in triangular
c      form for diagonalization)
c
      ij=0
      do i=1, N
        do j=1, i
           ij=ij+1
           if(i.eq.j) then
              H(ij)=GR(i,j)+V(i,j)
              VEC(i,j)=H(ij)
           else
              H(ij)=(GR(i,j)+GR(j,i))/two
              VEC(i,j)=H(ij)
              VEC(j,i)=H(ij)
           endif
        end do
      end do
c
c     DIAGONALIZE THE PSEDOSPECTRAL MATRIX FOR EIGENVALUES AND EIGENVECTORS
c
      call util_jacobi(N,VEC,N,E)
c
c     rearrange into descending order
c
      call eigsort(E,VEC,N,N)
c
c     find wavefunctions
c
      do i=1, N
         do j=1, N
            phi(i,j)=VEC(j,i)*RR(j,j)
         end do
      end do
c
      do i=1, N
         if (i.ge.N-NST+1) then
            ist=N-i+1
c
            area=zero
            do j=1, N/2
               area=area+phi(i,j)
            end do
c
            if ((mod(ist,2).eq.1).and.(area.lt.zero)) then
               do j=1, N
                  temp(j)=-phi(i,j)
               end do
            else if ((mod(ist,2).eq.0).and.(area.gt.zero)) then
               do j=1, N
                  temp(j)=-phi(i,j)
               end do
            else
               do j=1, N
                  temp(j)=phi(i,j)
               end do
            end if
c
            call norm (mode,temp,dx,nm,n)
c
            do j=1, N
               virtwave(mode,ist,j)=temp(j)
            end do
            virtE(mode,ist)=E(i)
c
            if (i.eq.(N-istate(mode))) then
               do j=1, N
                 twave(j)=temp(j)
               end do
               tE=E(i)
            end if
         end if
      end do
c
      do l=1, N
         wave(mode,l)=twave(l)
         allwave(index+1,mode,l)=twave(l)
      end do
c
      return
      end
c*module vscf    *deck norm
      subroutine norm (mode,twave,dx,nm,n)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c     this routine normalizes the wavefunction obtained from the 
c     collocation grid points for (i) th mode.
c
      dimension dx(NM)
      dimension twave(N)
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
C     CALCULATE THE NORMALIZATION CONSTANTS USING THE GRID WAVEFUNCTIONS
C
      wnorm=zero
      do m=1, N
         temp=dx(mode)*(twave(m)**2)
         wnorm=wnorm+temp
      end do
c
c     NORMALIZATION OF THE WAVEFUNCTIONS ON THE GRID POINTS
c
      do m=1, N
         twave(m)=(one/sqrt(wnorm))*twave(m)
      end do
C
C     CHECK FOR NORMALIZATION:
C
      area=zero
      do l=1, N
         area=area+dx(mode)*(twave(l)**2)
      end do
c
      return
      end
c*module vscf    *deck veffect
      subroutine Veffect (mode,vscf,wave,dx,Vc,Vtr,nm,n,ncoup)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      dimension vscf(N)
      dimension wave(NM,N)
      dimension dx(NM)
      dimension Vc(NM,NM,N,N)
      dimension Vtr(NM,NM,NM,N,N,N)
c
      PARAMETER (ZERO=0.0D+00)
c
c     this routine calculates the scf averaged potential 
c     on the GH-quadrature grid points for (i) th mode.
c
      do l=1, N
         vscf(l)=zero
         do j=1, NM
            if (j.ne.mode) then
               call SCFAVG(sum,wave,dx,Vc,mode,j,l,nm,n)
               vscf(l)=vscf(l)+sum
            end if
         end do
      end do
      if (ncoup.le.2) goto 100
c
      do l=1, N
         do j=1, NM-1
            do k=j+1, NM
              if (j.ne.mode .and. k.ne.mode) then
                 call SCFAVGT(sum,wave,dx,Vtr,mode,j,k,l,nm,n)
                 vscf(l)=vscf(l)+sum
              end if
            end do
         end do
      end do
c
  100 continue
      return
      end
c*module vscf    *deck scfavg
      subroutine SCFAVG(sum,wave,dx,Vc,i,j,l,nm,n)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     calculate the Vscf(i,l) for mode i and pt. l
c
      dimension dx(NM)
      dimension Vc(NM,NM,N,N)
      dimension wave(NM,N)
C
      PARAMETER (ZERO=0.0D+00)
C
      sum=zero
      do m=1, N
         temp1=dx(j)*Vc(i,j,l,m)*(wave(j,m)**2)
         sum=sum+temp1 
      end do
C
      return
      end
c*module vscf    *deck scfavgt
      subroutine SCFAVGT(sum,wave,dx,Vtr,i,j,k,l,nm,n)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     calculate the Vscf(i,l) for mode i and pt. l
c
      dimension dx(NM)
      dimension Vtr(NM,NM,NM,N,N,N)
      dimension wave(NM,N)
C
      PARAMETER (ZERO=0.0D+00)
C
      sum=zero
      do m1=1, N
         do m2=1, N
            temp1=dx(j)*dx(k)*Vtr(i,j,k,l,m1,m2)*
     *                        (wave(j,m1)**2)*(wave(k,m2)**2)
            sum=sum+temp1 
         end do
      end do
C
      return
      end
c*module vscf    *deck scfcorr
      subroutine scfcorr (wave,emp1,dx,Vc,Vtr,nm,n,ncoup)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     this routine calculates the scf correction 
c     on the GH-quadrature grid points. 
c
      dimension dx(NM)
      dimension Vc(NM,NM,N,N)
      dimension Vtr(NM,NM,NM,N,N,N)
      dimension wave(NM,N)
c
      PARAMETER (ZERO=0.0D+00)
c
      emp1=zero
      do i=1, NM-1
        do j=i+1, NM
          do l=1, N
            do m=1, N
              sum=dx(i)*dx(j)*Vc(i,j,l,m)*(wave(i,l)**2)*(wave(j,m)**2)
              emp1=emp1+sum
            end do
         end do
        end do
      end do
      if (ncoup.le.2) goto 100
c
      emp1t=zero
      do i=1, NM-2
         do j=i+1, NM-1
            do k=j+1, NM
               do l1=1, N
               do l2=1, N
               do l3=1, N
                  sum=dx(i)*dx(j)*dx(k)*Vtr(i,j,k,l1,l2,l3)*
     *                (wave(i,l1)**2)*(wave(j,l2)**2)*(wave(k,l3)**2)
                  emp1t=emp1t+sum
               end do
               end do
               end do
            end do
         end do
      end do
      emp1 = emp1 + emp1t*2.0D+00
c
  100 continue
      return
      end
c*module vscf    *deck virtstate
      subroutine virtstate (kstate,nm,nmax,nvst,ivst,ncoup,iexc) 
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      dimension IVST(NVST,NM)
c
      nost=0
c
c     generate single excitations
C
      do istate=0, NMAX
        do i=1, NM
          nost=nost+1
          do mode=1, NM
            if (mode.eq.i) then
              ivst(nost,mode)=istate
            else
              ivst(nost,mode)=0
            end if
          end do 
          if (nost.eq.1) go to 1000
        end do
1000  end do
c
c     generate double excitations
C
      do jstate=1, NMAX
         do istate=1, NMAX
c           isum=istate+jstate
            do i=1, NM-1
               do j=i+1, NM
                  nost=nost+1
                  do mode=1, NM
                     if (mode.eq.i) then
                        ivst(nost,mode)=istate
                     else if (mode.eq.j) then
                        ivst(nost,mode)=jstate
                     else if (mode.eq.kstate) then
                        ivst(nost,mode)=IEXC
                     else
                        ivst(nost,mode)=0
                     end if
                  end do
               end do
            end do
         end do
      end do
      if (ncoup.le.2) goto 2000
c
c     triple excitations if 3body coupling is included
C
      do istate=1, NMAX
         do jstate=1, NMAX
            do lstate=1, NMAX
c              isum=istate+jstate+lstate
               do i=1, NM-2
                  do j=i+1, NM-1
                     do l=j+1, NM
                        nost=nost+1
                        do mode=1, NM
                           if (mode.eq.i) then
                              ivst(nost,mode)=istate
                           else if (mode.eq.j) then
                              ivst(nost,mode)=jstate
                           else if (mode.eq.l) then
                              ivst(nost,mode)=lstate
                           else if (mode.eq.kstate) then
                              ivst(nost,mode)=IEXC
                           else
                              ivst(nost,mode)=0
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
c
2000  continue 
      return
      end
c*module vscf    *deck mppt
      subroutine mppt(emp2,emp0,wave,V,E,Vc,Vtr,Vscf,dx,ovrlp,
     *                nm,n,nst,nvst,ivst,iref,jref,jvirt,ncoup)
c
c     This program calculates the 2nd order MPPT energy correction.
c
      implicit double precision (a-h,o-z)
c
#include "global.fh"
#include "stdio.fh"
c
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      dimension wave(NM,NST,N)
      dimension Vc(NM,NM,N,N)
      dimension Vtr(NM,NM,NM,N,N,N)
      dimension Vscf(NM,N)
      dimension E(NM,NST)
      dimension dx(NM)
      dimension iref(NM)
      dimension IVST(NVST,NM)
      dimension jvirt(NVST,NM),jref(NM)
      dimension emp0(NVST)
      dimension v(NVST)
      dimension ovrlp(NM)
c
      PARAMETER (ZERO=0.0D+00, small=1.0D-05, cm=2.19474D+05)
c
c     reference state
c
      do i=1, NM
         iref(i)=iref(i)+1
      end do
c
c     compare virtual states with the reference state
c
       index=0
       do j=1, NVST  
          isum=0
          emp0(j)=zero
          index=index+1
          do i=1, NM
             jref(i)=ivst(index,i)
             jvirt(j,i)=jref(i)+1
             emp0(j)=emp0(j)+E(i,jvirt(j,i))
             if (iref(i).eq.jvirt(j,i)) isum=isum+1 
          end do
          if (isum.eq.NM) istate=j
       end do
c
c      generate the overlap between two states i and j
c
       do j=1, NVST
         v(j)=zero
         do jmode=1, NM
           jref(jmode)=jvirt(j,jmode)
         end do
         do k=1, NM-1
            do l=k+1, NM
               call AVG(sum,iref,jref,k,l,dx,wave,Vc,ovrlp,nm,n,nst) 
               v(j)=v(j)+sum
            end do
         end do 
         if (ncoup.gt.2) then
            do l1=1, NM-2
              do l2=l1+1, NM-1
                do l3=l2+1, NM
                  call AVG3(sum3,iref,jref,l1,l2,l3,dx,wave,Vtr,
     *                                               ovrlp,nm,n,nst) 
                  v(j)=v(j)+sum3
                end do
              end do
            end do 
         end if
         do k=1, NM
            call AVGHF(sum2,iref,jref,k,dx,wave,Vscf,ovrlp,nm,n,nst)
            v(j)=v(j)-sum2
         end do
       end do
c
c      1st order energy correction:
c
       emp1=v(istate)
c
c      2nd order energy correction:
c
       nkeep=0
       emp2=zero
       do j=1, NVST    !  state index
         if (j.ne.istate) then
            sumV=v(j)*v(j)
            sumE=emp0(istate)-emp0(j)
            if (abs(sumE).lt.small) then
               sumV=zero
            end if
            emp2=emp2+(sumV/sumE)
c
c      to cut down on the computation of coefficients for
c      2nd order wavefunction correction, select only the states
c      which contribute most to the energy correction and use these 
c      states as the basis.
c
            if (abs(sumV/sumE).gt.1.0D-08) nkeep=nkeep+1
         end if
       end do
       if (ga_nodeid().eq.0) then
          write(LuOut,9030) NVST
          write(LuOut,9040) nkeep
       endif
c
       return
c
 9030 format(1x,'Number of virtual states',I6)
 9040 format(1x,'Number of selected states',I6)
       end
c
c*module vscf    *deck avg
      subroutine AVG(sum,iref,jref,k,l,dx,wave,Vc,ovrlp,nm,n,nst)
c
c     this routine calculates the
c     avg = < Psi | Vc(k,l) | Psi >
c
      implicit double precision (a-h,o-z)
c
      dimension iref(NM),jref(NM)
      dimension ovrlp(NM)
      dimension wave(NM,NST,N)
      dimension Vc(NM,NM,N,N)
      dimension dx(NM)
c
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
c
      sum=zero
      do kl=1, N
        Si=dx(k)*wave(k,iref(k),kl)*wave(k,jref(k),kl)
        do ll=1, N
           Sj=dx(l)*wave(l,iref(l),ll)*wave(l,jref(l),ll) 
           sum=sum+Vc(k,l,kl,ll)*Si*Sj
        end do
      end do
c
      Sij=one
      do m=1, NM
        if ((m.eq.k).or.(m.eq.l)) then
           ovrlp(m)=one
        else
           ovrlp(m)=zero
           if (iref(m).eq.jref(m)) then
             do ll=1, N
               ovrlp(m)=ovrlp(m)+
     *                  dx(m)*wave(m,iref(m),ll)*wave(m,jref(m),ll)
             end do
           else
             Sij=zero
             goto 1000
           end if
        end if
        Sij=Sij*ovrlp(m)
      end do
1000  sum=sum*Sij
      return
      end
c*module vscf    *deck avg3
      subroutine AVG3(sum,iref,jref,l1,l2,l3,dx,wave,Vtr,
     *                                          ovrlp,nm,n,nst)
c
c     this routine calculates the
c     avg = < Psi | Vtrip(l1,l2,l3) | Psi >
c
      implicit double precision (a-h,o-z)
c
      dimension iref(NM),jref(NM)
      dimension ovrlp(NM)
      dimension wave(NM,NST,N)
      dimension Vtr(NM,NM,NM,N,N,N)
      dimension dx(NM)
c
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
c
      sum=zero
      do ll1=1, N
        S1=dx(l1)*wave(l1,iref(l1),ll1)*wave(l1,jref(l1),ll1)
        do ll2=1, N
          S2=dx(l2)*wave(l2,iref(l2),ll2)*wave(l2,jref(l2),ll2)
          do ll3=1, N
            S3=dx(l3)*wave(l3,iref(l3),ll3)*wave(l3,jref(l3),ll3) 
            sum=sum+Vtr(l1,l2,l3,ll1,ll2,ll3)*S1*S2*S3
          end do
        end do
      end do
c
      Sij=one
      do m=1, NM
        if ((m.eq.l1).or.(m.eq.l2).or.(m.eq.l3)) then
           ovrlp(m)=one
        else
           ovrlp(m)=zero
           if (iref(m).eq.jref(m)) then
             do ll=1, N
               ovrlp(m)=ovrlp(m)+
     *                  dx(m)*wave(m,iref(m),ll)*wave(m,jref(m),ll)
             end do
           else
             Sij=zero
             goto 1000
           end if
        end if
        Sij=Sij*ovrlp(m)
      end do
1000  sum=sum*Sij
      return
      end
c*module vscf    *deck avghf
      subroutine AVGHF(sum,iref,jref,k,dx,wave,Vscf,ovrlp,nm,n,nst)
c
c     this routine calculates the
c     avg = < Psi | Vhf(k) | Psi >
c
      implicit double precision (a-h,o-z)
c
      dimension iref(NM),jref(NM)
      dimension wave(NM,NST,N)
      dimension Vscf(NM,N)
      dimension dx(NM)
      dimension ovrlp(NM)
c
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
c
      sum=zero
      Sij=one
      do m=1, NM
         if (m.eq.k) then
           ovrlp(m)=one
           do ll=1, N
             temp=dx(m)*wave(m,jref(m),ll)*wave(m,iref(m),ll)
             sum=sum+Vscf(m,ll)*temp
           end do
         else
           ovrlp(m)=zero
           if (iref(m).eq.jref(m)) then
             do ll=1, N
               ovrlp(m)=ovrlp(m)+
     *              dx(m)*wave(m,jref(m),ll)*wave(m,iref(m),ll)
             end do
           else
             Sij=zero
             goto 1000
           end if
         end if
         Sij=Sij*ovrlp(m)
      end do
1000  sum=sum*Sij
      return
      end
c*module vscf    *deck eigsort
      SUBROUTINE EIGSORT(D,V,N,NP)
c
c     this routine resorts the eigenvalues from ascending 
c     to descending order, and rearranges the columns of V 
c     correspondingly. 
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION  D(NP),V(NP,NP)
C
      DO 130 I=1,N/2
         K=N-I+1
         P=D(K)
         D(K)=D(I)
         D(I)=P
         DO 120 J=1,N
            P=V(J,K)
            V(J,K)=V(J,I)
            V(J,I)=P
  120    CONTINUE
  130 CONTINUE
C
      RETURN
      END
c*module vscf    *deck endiff
      subroutine endiff(freq,frscf,frmp2,diag,emp1,emp2,nm,nstate,
     *                                                  ncoup,iexc)
C
      implicit double precision (a-h,o-z)
c
#include "global.fh"
#include "stdio.fh"
C
      dimension freq(nm), emp1(Nstate), emp2(Nstate), diag(Nstate)
      dimension frscf(nm), frmp2(nm) 
C
      PARAMETER (cm=2.19474D+05)
C
      if (ga_nodeid().eq.0) write(LuOut,*) 
      if (ga_nodeid().eq.0) 
     *    write(LuOut,*) " RESULTS OF VIBRATIONAL SCF CALCULATION:" 
     *                       ," Frequencies, cm-1" 
      if (ncoup.le.1) then
         if (ga_nodeid().eq.0) 
     *     write(LuOut,*)  " Mode      Harmonic     Diagonal" 
         do imode = 1, NM
            hfr = freq(imode)*cm
            if (iexc.gt.1) hfr=hfr*iexc
            ddiag = (diag(imode+1) - diag(1))*cm
            if (ga_nodeid().eq.0) write(LuOut,9010) imode, hfr, ddiag
            frscf(imode)=ddiag/cm
            frmp2(imode)=frscf(imode)
         end do
      else
         if (ga_nodeid().eq.0) write(LuOut,*) 
     *             " Mode      Harmonic     Diagonal      VSCF",
     *             "       PT2-VSCF"
         do imode = 1, NM
            hfr = freq(imode)*cm
            if (iexc.gt.1) hfr=hfr*iexc
            ddiag = (diag(imode+1) - diag(1))*cm
            dvscf = (emp1(imode+1) - emp1(1))*cm
            demp2 = (emp2(imode+1) - emp2(1))*cm
            if (ga_nodeid().eq.0) 
     *        write(LuOut,9020) imode,hfr,ddiag,dvscf,demp2
            frscf(imode)=dvscf/cm
            frmp2(imode)=demp2/cm
         end do
      end if
C
      return
C
 9010 format (2x,i2,5x,2(f10.2,3x))
 9020 format (2x,i2,5x,4(f10.2,3x))
      end
c*module vscf    *deck rddiag
      SUBROUTINE VSCF_RDDIAG(RQ,DQ,diagV,DMX,DMY,DMZ,NM,N,ndip,
     *                  nstart,ldiag)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "errquit.fh"
c
#include "nwc_const.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
      DIMENSION DQ(NM),RQ(NM,N),DIAGV(NM,N)
      DIMENSION DMX(NM,N),DMY(NM,N),DMZ(NM,N)
      dimension nstart(7)
C
      IF (ga_nodeid().eq.0) THEN
        if (ndip.eq.1) then
          do i=1, NM
            do l=1, N
              READ (ldiag,9012,end=9001,err=9001) RQ(i,l), diagV(i,l),
     *           dmx(i,l),dmy(i,l),dmz(i,l)
              nstart(2)=i
              nstart(3)=l
            enddo
            dq(i)=rq(i,2)-rq(i,1)
          enddo
        else
         do i=1, NM
            do l=1, N
              READ (ldiag,9010,end=9001,err=9001) RQ(i,l), diagV(i,l)
              nstart(2)=i
              nstart(3)=l
            enddo
            dq(i)=rq(i,2)-rq(i,1)
          enddo
        endif
 9001 if ((nstart(2).eq.nm).and.(nstart(3).eq.n)) nstart(1)=1
      ENDIF
c
c     Now broadcast information to all processors
c
      master = 0
      CALL ga_brdcst(msg_ns1,nstart,7*ma_sizeof(mt_int,1,mt_byte),
     *                  MASTER)
      CALL ga_brdcst(msg_rq,RQ,NM*N*ma_sizeof(mt_dbl,1,mt_byte),
     *                  MASTER)
      CALL ga_brdcst(msg_dq,DQ,NM*ma_sizeof(mt_dbl,1,mt_byte),
     *                  MASTER)
      CALL ga_brdcst(msg_diagv,DIAGV,
     *        NM*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
      if (ndip.eq.1) then
         CALL ga_brdcst(msg_dmx,DMX,
     *           NM*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
         CALL ga_brdcst(msg_dmy,DMY,
     *           NM*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
         CALL ga_brdcst(msg_dmz,DMZ,
     *           NM*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
      end if
C
      RETURN
 9010 format(5x,f14.8,2x,e14.7)
 9012 format(5x,f14.8,2x,e14.7,3(2x,e14.7))
      END
c*module vscf    *deck rdcoup
      SUBROUTINE RDCOUP(coupV,NM,N,NSTART,dm2x,dm2y,dm2z,lcoup,ndip)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
C
      DIMENSION COUPV(NM,NM,N,N)
      DIMENSION dm2x(NM,NM,N,N)
      DIMENSION dm2y(NM,NM,N,N)
      DIMENSION dm2z(NM,NM,N,N)
      dimension nstart(7)
C
      IF (ga_nodeid().eq.0) THEN
        if (ndip.eq.1) then
         do i=1, NM
            do j=i+1, NM
               do l=1, N
                  do m=1, N
                     read(lcoup,9022,end=9019,err=9019) coupV(i,j,l,m),
     *                       dm2x(i,j,l,m),dm2y(i,j,l,m),dm2z(i,j,l,m)
                     nstart(2)=i
                     nstart(3)=l
                     nstart(4)=j
                     nstart(5)=m
                  enddo
               enddo
            enddo
         enddo
        else
         do i=1, NM
            do j=i+1, NM
               do l=1, N
                  do m=1, N
                     read(lcoup,9020,end=9019,err=9019) coupV(i,j,l,m)
                     nstart(2)=i
                     nstart(3)=l
                     nstart(4)=j
                     nstart(5)=m
                  enddo
               enddo
            enddo
         enddo
        end if
 9019 if ((nstart(2).eq.nm).and.(nstart(3).eq.n).and.(nstart(4).eq.nm)
     &    .and.(nstart(5).eq.n)) nstart(1)=2
      END IF
c
c     Now broadcast information to all processors
c
      master = 0
      CALL ga_brdcst(msg_ns2,nstart,7*ma_sizeof(mt_int,1,mt_byte),
     &               MASTER)
      CALL ga_brdcst(msg_coupv,COUPV,
     *      NM*NM*N*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
      if (ndip.eq.1) then
        CALL ga_brdcst(msg_dm2x,dm2x,
     *      NM*NM*N*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
        CALL ga_brdcst(msg_dm2y,dm2y,
     *      NM*NM*N*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
        CALL ga_brdcst(msg_dm2z,dm2z,
     *      NM*NM*N*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
      end if
     
C
      RETURN
 9020 format(41x,e14.7) 
 9022 format(41x,e14.7,3(e14.7,2x))
      END
c*module vscf    *deck rdtrip
      SUBROUTINE RDTRIP(tripV,NM,N,NSTART,ltrip)
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
#include "global.fh"
 
      integer ltrip
      DIMENSION TRIPV(NM,NM,NM,N,N,N)
      dimension nstart(7)
 
      IF (ga_nodeid().eq.0) THEN
         do i=1, NM 
            do j=i+1, NM
               do k=j+1, NM
                  do l1=1, N
                     do l2=1, N
                        do l3=1, N
                           read(ltrip,9025,end=9024,err=9024) 
     &                           tripV(i,j,k,l1,l2,l3)
                           nstart(2)=i
                           nstart(3)=l1
                           nstart(4)=j
                           nstart(5)=l2
                           nstart(6)=k
                           nstart(7)=l3
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
 9024 if ((nstart(2).eq.nm).and.(nstart(3).eq.n).and.(nstart(4).eq.nm)
     &    .and.(nstart(5).eq.n).and.(nstart(6).eq.nm)
     &    .and.(nstart(7).eq.n)) nstart(1)=3
      END IF
c
c     Now broadcast information to all processors
c
      master = 0
      CALL ga_brdcst(msg_ns3,nstart,7*ma_sizeof(mt_int,1,mt_byte),
     &               MASTER)
      CALL ga_brdcst(msg_tripv,TRIPV,
     *      NM*NM*NM*N*N*N*ma_sizeof(mt_dbl,1,mt_byte),MASTER)
      RETURN
 9025 format(61x,e14.7)
      END
c*module vscf    *deck intens
      SUBROUTINE INTENS(SINTSCF,UX,UY,UZ,FRSCF,FRMP2,DQ,WAVE,NM,N)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
#include "global.fh"
#include "stdio.fh"
C
c     THIS ROUTINE CALCULATES IR intensity for i-th MODE
c     USING THE following equation:
c
c     Intensity(i) = |<ground psi(i)| u(i) |excited psi(i)>|^2
c
c             where   psi(i) = wavefunction corresponding to mode i
c                     u(i)   = dipole moment along normal coordinate i
c
      dimension dQ(NM)
      dimension wave(NM+1,NM,N)            ! VSCF wavefunctions
      dimension frscf(NM)                  ! VSCF frequencies
      dimension frmp2(NM)                  ! VMP2 frequencies
      dimension ux(NM,N),uy(NM,N),uz(NM,N) ! dipole moment components
      dimension sIntscf(NM)                ! VSCF Intensity
c
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (const=2.5048D+00, cm=2.19474D+05) 
c
c     const = 8*pi^3*Navogadro*(E-18)^2*E-5/(3*h*c) in CGS
c
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) " IR INTENSITIES CALCULATED USING DIPOLE MOMENT:"
        write(LuOut,*) " mode   frequency, cm-1  intensity, km/mol" 
      endif
c
c     Note: dipole assumed to be in debeye
c
      do j=1, NM
         sIntx = zero 
         sInty = zero 
         sIntz = zero 
         do l=1, N
            sIntx = sIntx +
     *                 ux(j,l)*wave(1,j,l)*wave(j+1,j,l)*dQ(j)
            sInty = sInty +
     *                 uy(j,l)*wave(1,j,l)*wave(j+1,j,l)*dQ(j)
            sIntz = sIntz +
     *                 uz(j,l)*wave(1,j,l)*wave(j+1,j,l)*dQ(j)
         end do
         sIntscf(j) = sIntx*sIntx+sInty*sInty+sIntz*sIntz
         sIntscf(j) = frscf(j)*cm*const*sIntscf(j)
         if (ga_nodeid().eq.0) 
     *     write(LuOut,9000) j, frmp2(j)*cm, sIntscf(j)
      end do
c
      if (ga_nodeid().eq.0) then
         write(LuOut,*) 
         write(LuOut,*) 
     *   " (Note: SCF dipole moments used to calculate intensities)"
      endif
      return
c
 9000 format (2x,i2,5x,f10.2,4x,f10.2)
      end
c*module vscf    *deck dintens
      SUBROUTINE DINTENS(SINTSCF,DDM,DDER,FRSCF,FRMP2,
     *                   VEC,Q,DQ,WAVE,NM,N,NC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
#include "global.fh"
#include "stdio.fh"
C
c     THIS ROUTINE CALCULATES IR intensity for i-th MODE
c     USING THE following equation:
c
c     Intensity(i) = ( du/dQ(i) |<gr. psi(i)| Q(i) |ex. psi(i)>| )^2
c
c             where   du/dQ(i) = dipole derivative at equilibrium
c
      dimension dQ(NM)
      dimension Q(NM,N)
      dimension VEC(NC,NC)
      dimension wave(NM+1,NM,N)        ! VSCF wavefunctions
      dimension frscf(NM)              ! VSCF frequencies
      dimension frmp2(NM)              ! VMP2 frequencies
      dimension sIntscf(NM)            ! VSCF Intensity
      dimension ddm(NC*3)              ! dipole derivative matrix in cartezians
      dimension dder(NM)               ! square of dipole moment derivative
c
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (AMU=1.8229D+03, cm=2.19474D+05)
      PARAMETER (const=2.5048D+00, bohr=5.2918D-01)
c
c     const = 8*pi^3*Navogadro*(E-18)^2*E-5/(3*h*c) in CGS
c
      FACT = one/AMU
      NSTART=NC-NM+1
C
      if (ga_nodeid().eq.0) write(LuOut,*) 
      if (ga_nodeid().eq.0) write(LuOut,*) 
     * " IR INTENSITIES CALCULATED USING DIPOLE DERIVATIVE:" 
c
      do i = NSTART,NC
         im=nc-i+1
         DDX = DDOT(NC,VEC(1,i),1,DDM(1),3)
         DDY = DDOT(NC,VEC(1,i),1,DDM(2),3)
         DDZ = DDOT(NC,VEC(1,i),1,DDM(3),3)
         dder(im) = DDX*DDX + DDY*DDY + DDZ*DDZ
      enddo
c
      if (ga_nodeid().eq.0) write(LuOut,*) 
     *   " mode   frequency, cm-1  intensity, km/mol" 
c
      do j=1, NM
         sIntscf(j) = zero 
         do l=1, N
            sIntscf(j) = sIntscf(j) +
     *                   Q(j,l)*wave(1,j,l)*wave(j+1,j,l)*dQ(j)
         end do
c
         sIntscf(j) = frscf(j)*cm*sIntscf(j)**2
         sIntscf(j) = sIntscf(j)*dder(j)*const*fact*bohr*bohr
         if (ga_nodeid().eq.0) 
     *     write(LuOut,9000) j, frmp2(j)*cm, sIntscf(j)
c
      end do
c
      return
c
 9000 format (2x,i2,5x,f10.2,4x,f10.2)
      end
C*MODULE MTHLIB  *DECK DGEFA
      SUBROUTINE DGEFA_off(A,LDA,N,IPVT,INFO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,*),IPVT(*)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D+00) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D+00/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D+00) INFO = N
      RETURN
      END
C*MODULE MTHLIB  *DECK DGEDI
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,*),DET(2),WORK(*),IPVT(*)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN ABS,MOD
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D+00
         DET(2) = 0.0D+00
         TEN = 10.0D+00
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D+00) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0D+00) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D+00
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D+00
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D+00/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D+00
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D+00
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
C*MODULE MTHLIB  *DECK MRARBR
      SUBROUTINE MRARBR(A,LDA,NA,MA,B,LDB,MB,C,LDC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "errquit.fh"
C
#include "global.fh"
c
      DIMENSION A(LDA,MA),B(LDB,MB),C(LDC,MB)
C
      PARAMETER (ZERO=0.0D+00)
*BL3  PARAMETER (ONE=1.0D+00)
C
C     ----- SEQUENTIAL MATRIX MULTIPLY: C = A * B -----
C     ----- SEE ALSO -PRARBR- ROUTINE -----
C        A    - THE INPUT REAL NA X MA RECTANGULAR MATRIX
C        B    - THE INPUT REAL MA X MB RECTANGULAR MATRIX
C        C    - THE OUTPUT PRODUCT NA X MB MATRIX
C        LDA  - THE LEADING DIMENSION OF MATRIX A
C        NA   - THE EFFECTIVE ROW DIMENSION OF MATRICES A AND C
C        MA   - THE COLUMN DIMENSION OF MATRIX A,
C               AND EFFECTIVE ROW DIMENSION OF MATRIX B
C        LDB  - THE LEADING DIMENSION OF MATRIX B
C        MB   - THE COLUMN DIMENSION OF MATRICES B AND C
C        LDC  - THE LEADING DIMENSION OF MATRIX C
C     AUTHOR = STEVE ELBERT, 31 OCT 1979
C
      IF(LDA.LT.NA .OR. LDB.LT.MA .OR. LDC.LT.NA) GO TO 800
C
*BL3  CALL DGEMM('N','N',NA,MB,MA,ONE,A,LDA,B,LDB,ZERO,C,LDC)
*BL3  IF(ONE.NE.ZERO) RETURN
C
      M=MA
      IF(MOD(M,2).NE.0) M=M-1
C
      DO 300 I=1,NA
         DO 200 J=1,MB
            CIJ=ZERO
            IF(M.NE.MA) CIJ=A(I,MA)*B(MA,J)
            IF(MA.GT.1) THEN
               DO 100 K=1,M,2
                  CIJ=CIJ + A(I,K)*B(K,J) + A(I,K+1)*B(K+1,J)
  100          CONTINUE
            END IF
            C(I,J)=CIJ
  200    CONTINUE
  300 CONTINUE
      RETURN
C
  800 CONTINUE
      IF(ga_nodeid().eq.0) WRITE(6,900) LDA,NA,MA,LDB,MB,LDC
      CALL errquit("mrarbr: bad value of a leading dimension",555,
     &       INPUT_ERR)
      STOP
C
  900 FORMAT(/1X,'ERROR IN CALL TO -MRARBR-'/
     *        1X,'LDA,NA,MA,LDB,MB,LDC=',6I10)
      END
c
      subroutine vscf_restart(rtdb,restart,ncoup,nstart,ldiag,lcoup,
     &                            ltrip,nnm,ngrid,rq,dq,diagv,dm1x,dm1y,
     &                            dm1z,dm2x,dm2y,dm2z,coupv,tripv,ndip)
      implicit double precision(a-h,o-z)
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      dimension rq(nnm,ngrid),dq(nnm),diagv(nnm,ngrid)
      dimension dm1x(nnm,ngrid),dm1y(nnm,ngrid),dm1z(nnm,ngrid)
      dimension dm2x(nnm,nnm,ngrid,ngrid),dm2y(nnm,nnm,ngrid,ngrid)
      dimension dm2z(nnm,nnm,ngrid,ngrid),coupv(nnm,nnm,ngrid,ngrid)
      dimension tripv(nnm,nnm,nnm,ngrid,ngrid,ngrid)
c
      character*255 namediag, namecoup, nametrip
      dimension nstart(7)
      logical restart,diag_exists,coup_exists,trip_exists
      integer rtdb
c
      do ii=1,7
        nstart(ii) = -1
      enddo
c
c     Open files for diagonal, coupling, and triple potentials 
c
      ldiag = 81
      lcoup = 82
      ltrip = 83
c
      if (ga_nodeid().eq.0) then
        diag_exists=.false.
        coup_exists=.false.
        trip_exists=.false.
        if (ncoup.ge.1) then
          call util_file_name('diag',.false.,.false.,namediag)
          inquire(file=namediag,exist=diag_exists)
          OPEN (UNIT=ldiag, FORM='FORMATTED', FILE=namediag,
     *          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
          REWIND (UNIT=ldiag)
        endif
        if (ncoup.ge.2) then
          call util_file_name('coup',.false.,.false.,namecoup)
          inquire(file=namecoup,exist=coup_exists)
          OPEN (UNIT=lcoup, FORM='FORMATTED', FILE=namecoup,
     *          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
          REWIND (UNIT=lcoup)
        endif
        if (ncoup.ge.3) then
          call util_file_name('trip',.false.,.false.,nametrip)
          inquire(file=nametrip,exist=trip_exists)
          OPEN (UNIT=ltrip, FORM='FORMATTED', FILE=nametrip,
     *          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
          REWIND (UNIT=ltrip)
        endif
      end if
c
c     If there is restart data, make sure all the processors
c     are aware of it and get it from the files.
c
      master=0
      CALL ga_brdcst(msg_diagv,diag_exists,
     *        ma_sizeof(mt_log,1,mt_byte),MASTER)
      CALL ga_brdcst(msg_diagv,coup_exists,
     *        ma_sizeof(mt_log,1,mt_byte),MASTER)
      CALL ga_brdcst(msg_diagv,trip_exists,
     *        ma_sizeof(mt_log,1,mt_byte),MASTER)
      if (diag_exists) then
         if(ga_nodeid().eq.0) write(LuOut,9030)
         CALL VSCF_RDDIAG(rq,dq,diagv,dm1x,dm1y,dm1z,nnm,ngrid,ndip, 
     *               nstart,ldiag)
      endif
      if (coup_exists.and.ncoup.ge.2) then
         if(ga_nodeid().eq.0) write(LuOut,9040)
         CALL RDCOUP(coupv,nnm,ngrid,nstart,dm2x,dm2y,dm2z,lcoup,ndip)
      endif
      if (trip_exists.and.ncoup.eq.3) then
         if(ga_nodeid().eq.0) write(LuOut,9050)
         CALL RDTRIP(TRIPV,NNM,NGRID,nstart,ltrip)
      endif
c
c     Restart if user specified restart keyword and we are reusing 
c     the old rtdb
c
      if (.not. rtdb_get(rtdb,'vscf:restart',mt_log,1,restart))
     &   restart=.false.
c
      return
 9030 format(/1x,'Reading diagonal energy data for restart...')
 9040 format(/1x,'Reading coupling energy data for restart...')
 9050 format(/1x,'Reading 3-body energy data for restart...')
      end
