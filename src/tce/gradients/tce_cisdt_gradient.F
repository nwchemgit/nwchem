      subroutine tce_cisdt_gradient(d_f1,d_v2,d_t1,d_t2,d_t3,
     1     k_f1_offset,k_v2_offset,
     1     k_t1_offset,k_t2_offset,k_t3_offset,
     1     size_t1,size_t2,size_t3,rtdb)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "rtdb.fh"
c     
      integer d_f1,k_f1_offset
      integer d_v2,k_v2_offset
      integer d_t1,k_t1_offset,size_t1
      integer d_t2,k_t2_offset,size_t2
      integer d_t3,k_t3_offset,size_t3
c
      integer d_t0,k_t0_offset,l_t0_offset,size_t0
      integer d_y0,k_y0_offset,l_y0_offset,size_y0
      integer d_y1,k_y1_offset,l_y1_offset,size_y1
      integer d_y2,k_y2_offset,l_y2_offset,size_y2
      integer d_y3,k_y3_offset,l_y3_offset,size_y3
      integer d_d0,k_d0_offset,l_d0_offset,size_d0
c
c     integer k_force
c     
      integer rtdb
      integer ii
      double precision corr,denominator
      double precision corr_a,corr_b,corr_c
      double precision corr_d,corr_e,corr_f
      double precision corr_g,corr_h,corr_i
      character*255 filename
c
c     variables for 1PDM, 2PDM, and 1PWDM
c     ====
c     1PDM
c     ====
      integer d_hh,l_hh_offset,k_hh_offset,size_hh
      integer d_pp,l_pp_offset,k_pp_offset,size_pp
      integer d_ph,l_ph_offset,k_ph_offset,size_ph
      integer d_1pdm,l_1pdm_offset,k_1pdm_offset,size_1pdm
c     ====
c     2PDM
c     ====
      integer d_hhhh,l_hhhh_offset,k_hhhh_offset,size_hhhh
      integer d_hhhp,l_hhhp_offset,k_hhhp_offset,size_hhhp
      integer d_hhpp,l_hhpp_offset,k_hhpp_offset,size_hhpp
      integer d_hphh,l_hphh_offset,k_hphh_offset,size_hphh
      integer d_hphp,l_hphp_offset,k_hphp_offset,size_hphp
      integer d_hppp,l_hppp_offset,k_hppp_offset,size_hppp
      integer d_pphh,l_pphh_offset,k_pphh_offset,size_pphh
      integer d_pphp,l_pphp_offset,k_pphp_offset,size_pphp
      integer d_pppp,l_pppp_offset,k_pppp_offset,size_pppp
      integer d_2pdm,l_2pdm_offset,k_2pdm_offset,size_2pdm
c     =====
c     1PWDM
c     ===== 
      integer d_whh,l_whh_offset,k_whh_offset,size_whh
      integer d_whp,l_whp_offset,k_whp_offset,size_whp
      integer d_wpp,l_wpp_offset,k_wpp_offset,size_wpp
      integer d_1pwdm,l_1pwdm_offset,k_1pwdm_offset,size_1pwdm
c     ===========================
c     add references contribution
c     ===========================
      integer l_a,k_a,l_b,k_b,dima
c     ===========================
c     atomic orbital tiles
c     ===========================
c      integer atpart2
c      integer nalength2(200)
c      integer a2length2(200)
c     ======
c     Z-VEC
c     ======
      integer d_x1,l_x1_offset,k_x1_offset,size_x1
      integer d_r1,l_r1_offset,k_r1_offset,size_r1
c
      double precision cpu,wall,r1,residual
      integer dummy, iter
c     
c     memory allocation
c     
      integer nh,np
c     
c     mos
c     
      integer dim_mo
      integer l_mo,k_mo
c     
c     1pdm_ao, 2pdm_ao, and 1pwdm_ao
c     
      integer dim_1pdm_ao,l_1pdm_ao,k_1pdm_ao
      integer dim_2pdm_ao,l_2pdm_ao,k_2pdm_ao
      integer dim_1pwdm_ao,l_1pwdm_ao,k_1pwdm_ao
c     ===========================
c     GA for 1PDM_AO and 1PWDM_AO
c     ===========================
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer  g_1pdm, g_1pwdm
c     
      logical tce_gradients_doit
      external tce_gradients_doit
      logical status
c     
      integer i
      logical nodezero
c     
      integer nbf_verify        ! nbf read from movecs_read_header
      integer ipol_verify       ! ipol read from movecs_read_header
      character*255 movecs_out  ! MO vector file name
      character*255 title       ! Title
      character*255 basisname   ! Basis set name
      character*255 scftype     ! Type of SCF
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
c      integer ga_create_atom_blocked
c      external ga_create_atom_blocked
      logical rohf
      integer j
c     
      integer l_test,k_test,l_r2,k_r2,size,l_r1,k_r1
      integer g1,g2,g3,g4,g1b,g2b,g3b,g4b
      integer ig1,ig2,ig3,ig4,igg1,igg2,igg3,igg4
      integer k,l,ind1,ind2
      integer dim_2pdm, k_2pdm, l_2pdm
      integer dim_1pdm, k_1pdm, l_1pdm
      integer dim_1pwdm,k_1pwdm,l_1pwdm
      integer k_2eint,l_2eint,k_1eint,l_1eint
c
      double precision res
      integer l_zvec_x,k_zvec_x
      integer m,e,a,b,q,r,s,n,f
c
c     debug
      logical zvec_debug, wdm_debug
c
c 
      nodezero=(ga_nodeid().eq.0)
c
      irrep_y=irrep_t

      call tce_filename('t0',filename)
      call tce_e_offset(l_t0_offset,k_t0_offset,size_t0)
      call createfile(filename,d_t0,size_t0)
      call put_block(d_t0,1.d0,1,0)
      call reconcilefile(d_t0,size_t0)
c
      call tce_filename('y0',filename)
      call tce_e_offset(l_y0_offset,k_y0_offset,size_y0)
      call createfile(filename,d_y0,size_y0)
      call put_block(d_y0,1.d0,1,0)
      call reconcilefile(d_y0,size_y0)
c
      call tce_filename('y1',filename)
      call tce_y1_offset(l_y1_offset,k_y1_offset,size_y1)
      call createfile(filename,d_y1,size_y1)
      call trans_t1_to_t1dagger(d_t1,k_t1_offset,d_y1,k_y1_offset)
      call reconcilefile(d_y1,size_y1)
c
      call tce_filename('y2',filename)
      call tce_y2_offset(l_y2_offset,k_y2_offset,size_y2)
      call createfile(filename,d_y2,size_y2)
      call trans_t2_to_t2dagger(d_t2,k_t2_offset,d_y2,k_y2_offset)
      call reconcilefile(d_y2,size_y2)
c
      call tce_filename('y3',filename)
      call tce_y3_offset(l_y3_offset,k_y3_offset,size_y3)
      call createfile(filename,d_y3,size_y3)
      call trans_t3_to_t3dagger(d_t3,k_t3_offset,d_y3,k_y3_offset)
      call reconcilefile(d_y3,size_y3)
c
      call tce_filename('d0',filename)
      call tce_e_offset(l_d0_offset,k_d0_offset,size_d0)
      call createfile(filename,d_d0,size_d0)
      call tce_cisdt_denominator(d_d0,d_t1,d_t2,d_t3,d_y1,d_y2,d_y3,
     1     k_d0_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset,k_y3_offset)
      call reconcilefile(d_d0,size_d0)
      call get_block(d_d0,denominator,1,0)
      denominator=denominator+1.d0
c
      call dscalfile(1.d0/dsqrt(denominator),d_t0,size_t0)
      call dscalfile(1.d0/dsqrt(denominator),d_t1,size_t1)
      call dscalfile(1.d0/dsqrt(denominator),d_t2,size_t2)
      call dscalfile(1.d0/dsqrt(denominator),d_t3,size_t3)
      call dscalfile(1.d0/dsqrt(denominator),d_y0,size_y0)
      call dscalfile(1.d0/dsqrt(denominator),d_y1,size_y1)
      call dscalfile(1.d0/dsqrt(denominator),d_y2,size_y2)
      call dscalfile(1.d0/dsqrt(denominator),d_y3,size_y3)
c    
c     ====================================
c     1PDM, 2PDM, and 1PWDM with MO basis
c     ====================================
c
      print *,'1PDM'
c
c     
c     ====================
c     1PDM Hole-Hole Block
c     ====================
      call tce_filename('hh',filename)
      
      call tce_dens_hh_offset(l_hh_offset,k_hh_offset,size_hh)
      call createfile(filename,d_hh,size_hh)
      call cisdt_1pdm_hh_mo(d_hh,d_t1,d_t2,d_t3,d_y1,d_y2,d_y3,
     1     k_hh_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset,k_y3_offset)
      call reconcilefile(d_hh,size_hh)
      print *,'1PDM Hole-Hole Block'
c     ============================
c     1PDM Particle-Particle Block
c     ============================
      call tce_filename('pp',filename)
      call tce_dens_pp_offset(l_pp_offset,k_pp_offset,size_pp)
      call createfile(filename,d_pp,size_pp)
      call cisdt_1pdm_pp_mo(d_pp,d_t1,d_t2,d_t3,d_y1,d_y2,d_y3,
     1     k_pp_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset,k_y3_offset)
      call reconcilefile(d_pp,size_pp)
c
      print *,'2PDM'
c     ==============================
c     2PDM Hole-Hole-Hole-Hole Block
c     ==============================
      call tce_filename('hhhh',filename)
      call tce_dens_hhhh_offset(l_hhhh_offset,k_hhhh_offset,size_hhhh)
      call createfile(filename,d_hhhh,size_hhhh)
      call cisdt_2pdm_hhhh_mo(d_hhhh,d_t2,d_t3,d_y2,d_y3,
     1     k_hhhh_offset,k_t2_offset,k_t3_offset,
     1     k_y2_offset,k_y3_offset)
      call reconcilefile(d_hhhh,size_hhhh)
c     ==================================
c     2PDM Hole-Hole-Hole-Particle Block
c     ==================================
      call tce_filename('hhhp',filename)
      call tce_dens_hhhp_offset(l_hhhp_offset,k_hhhp_offset,size_hhhp)
      call createfile(filename,d_hhhp,size_hhhp)
      call cisdt_2pdm_hhhp_mo(d_hhhp,d_t1,d_t2,d_y2,d_y3,
     1     k_hhhp_offset,k_t1_offset,k_t2_offset,
     1     k_y2_offset,k_y3_offset)
      call reconcilefile(d_hhhp,size_hhhp)
c     ======================================
c     2PDM Hole-Hole-Particle-Particle Block
c     ======================================
      call tce_filename('hhpp',filename)
      call tce_dens_hhpp_offset(l_hhpp_offset,k_hhpp_offset,size_hhpp)
      call createfile(filename,d_hhpp,size_hhpp)
      call cisdt_2pdm_hhpp_mo(d_hhpp,d_t0,d_t1,d_y2,d_y3,
     1     k_hhpp_offset,k_t0_offset,k_t1_offset,
     1     k_y2_offset,k_y3_offset)
      call reconcilefile(d_hhpp,size_hhpp)
c     ==================================
c     2PDM Hole-Particle-Hole-Hole Block
c     ==================================
      call tce_filename('hphh',filename)
      call tce_dens_hphh_offset(l_hphh_offset,k_hphh_offset,size_hphh)
      call createfile(filename,d_hphh,size_hphh)
      call cisdt_2pdm_hphh_mo(d_hphh,d_t2,d_t3,d_y1,d_y2,
     1     k_hphh_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset)
      call reconcilefile(d_hphh,size_hphh)
c     ======================================
c     2PDM Hole-Particle-Hole-Particle Block
c     ======================================
      call tce_filename('hphp',filename)
      call tce_dens_hphp_offset(l_hphp_offset,k_hphp_offset,size_hphp)
      call createfile(filename,d_hphp,size_hphp)
      call cisdt_2pdm_hphp_mo(d_hphp,d_t1,d_t2,d_t3,d_y1,d_y2,d_y3,
     1     k_hphp_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset,k_y3_offset)
      call reconcilefile(d_hphp,size_hphp)
c     ==========================================
c     2PDM Hole-Particle-Particle-Particle Block
c     ==========================================
      call tce_filename('hppp',filename)
      call tce_dens_hppp_offset(l_hppp_offset,k_hppp_offset,size_hppp)
      call createfile(filename,d_hppp,size_hppp)
      call cisdt_2pdm_hppp_mo(d_hppp,d_t1,d_t2,d_y2,d_y3,
     1     k_hppp_offset,k_t1_offset,k_t2_offset,
     1     k_y2_offset,k_y3_offset)
      call reconcilefile(d_hppp,size_hppp)
c     ======================================
c     2PDM Particle-Particle-Hole-Hole Block
c     ======================================
      call tce_filename('pphh',filename)
      call tce_dens_pphh_offset(l_pphh_offset,k_pphh_offset,size_pphh)
      call createfile(filename,d_pphh,size_pphh)
      call cisdt_2pdm_pphh_mo(d_pphh,d_t2,d_t3,d_y0,d_y1,
     1     k_pphh_offset,k_t2_offset,k_t3_offset,
     1     k_y0_offset,k_y1_offset)
      call reconcilefile(d_pphh,size_pphh)
c     ==========================================
c     2PDM Particle-Particle-Hole-Particle Block
c     ==========================================
      call tce_filename('pphp',filename)
      call tce_dens_pphp_offset(l_pphp_offset,k_pphp_offset,size_pphp)
      call createfile(filename,d_pphp,size_pphp)
      call cisdt_2pdm_pphp_mo(d_pphp,d_t2,d_t3,d_y1,d_y2,
     1     k_pphp_offset,k_t2_offset,k_t3_offset,
     1     k_y1_offset,k_y2_offset)
      call reconcilefile(d_pphp,size_pphp)
c     ==============================================
c     2PDM Particle-Particle-Particle-Particle Block
c     ==============================================
      call tce_filename('pppp',filename)
      call tce_dens_pppp_offset(l_pppp_offset,k_pppp_offset,size_pppp)
      call createfile(filename,d_pppp,size_pppp)
      call cisdt_2pdm_pppp_mo(d_pppp,d_t2,d_t3,d_y2,d_y3,
     1     k_pppp_offset,k_t2_offset,k_t3_offset,
     1     k_y2_offset,k_y3_offset)
      call reconcilefile(d_pppp,size_pppp)
c
c
      call tce_filename('x1',filename)
      call tce_zvec_x1_offset(l_x1_offset,k_x1_offset,size_x1)
      call createfile(filename,d_x1,size_x1)
      call tce_zero(d_x1,size_x1)
      call zvec_x_a_a(d_hhhh,d_x1,d_v2,k_hhhh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_b_a(d_hhhp,d_x1,d_v2,k_hhhp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_b_b(d_hphh,d_x1,d_v2,k_hphh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_c_a(d_hphh,d_x1,d_v2,k_hphh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_c_b(d_hhhp,d_x1,d_v2,k_hhhp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_d(d_hphp,d_x1,d_v2,k_hphp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_e_a(d_pphh,d_x1,d_v2,k_pphh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_e_b(d_hhpp,d_x1,d_v2,k_hhpp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_f_a(d_pphp,d_x1,d_v2,k_pphp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_a_f_b(d_hppp,d_x1,d_v2,k_hppp_offset,
     1     k_x1_offset,k_v2_offset)
c
      call zvec_x_b_a_a(d_hphh,d_x1,d_v2,k_hphh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_a_b(d_hhhp,d_x1,d_v2,k_hhhp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_b_a(d_pphh,d_x1,d_v2,k_pphh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_b_b(d_hhpp,d_x1,d_v2,k_hhpp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_c(d_hphp,d_x1,d_v2,k_hphp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_d_a(d_pphp,d_x1,d_v2,k_pphp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_d_b(d_hppp,d_x1,d_v2,k_hppp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_e_a(d_hppp,d_x1,d_v2,k_hppp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_e_b(d_pphp,d_x1,d_v2,k_pphp_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_b_f(d_pppp,d_x1,d_v2,k_pppp_offset,
     1     k_x1_offset,k_v2_offset)
c
      call zvec_x_c_a(d_hh,d_x1,d_v2,k_hh_offset,
     1     k_x1_offset,k_v2_offset)
      call zvec_x_c_b(d_pp,d_x1,d_v2,k_pp_offset,
     1     k_x1_offset,k_v2_offset)
c     ==============================================================================
c     1PDM Particle-Hole Block (Hole-Particle Block is equal to Particle-Hole Block)
c     ==============================================================================
      call tce_filename('ph',filename)
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,size_ph)
      call createfile(filename,d_ph,size_ph)
      call tce_diis_init()
      call tce_dens_ph_offset(l_r1_offset,k_r1_offset,size_r1)
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        call tce_filename('r1',filename)
        call createfile(filename,d_r1,size_r1)
        call ccsd_zvec(d_f1,d_r1,d_v2,d_x1,d_ph,k_f1_offset,
     1         k_r1_offset,k_v2_offset,k_x1_offset,k_ph_offset)
        call reconcilefile(d_r1,size_r1)
        call tce_residual_t1(d_r1,k_r1_offset,r1)
        residual = r1
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(6,*) iter,residual
        if (residual .lt. thresh) then
          call deletefile(d_r1)
          call tce_diis_tidy()
          goto 320
        endif
        call tce_diis(.false.,iter,.true.,.false.,.false.,.false.,
     1                  d_r1,d_ph,k_ph_offset,size_ph,
     2                  dummy,dummy,dummy,dummy,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
        call deletefile(d_r1)
      enddo
      call errquit('tce_cisdt_gradient: maxiter exceeded',iter,CALC_ERR)
c
320   continue
      if(.not.ma_pop_stack(l_r1_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',120,ma_err)
c     =====================
c     1PWDM Hole-Hole Block
c     =====================
      call tce_filename('whh',filename)
      call tce_dens_hh_offset(l_whh_offset,k_whh_offset,size_whh)
      call createfile(filename,d_whh,size_whh)
      call tce_zero(d_whh,size_whh)
      call wdm_hh_mo_a(d_hh,d_f1,d_whh,k_hh_offset,k_f1_offset,
     1     k_whh_offset)
      call wdm_hh_mo_b_a_a(d_hh,d_whh,d_v2,k_hh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_b_a_b(d_hh,d_whh,d_v2,k_hh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_b_b_a(d_pp,d_whh,d_v2,k_pp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_b_b_b(d_pp,d_whh,d_v2,k_pp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_a(d_hhhh,d_whh,d_v2,k_hhhh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_b_a(d_hhhp,d_whh,d_v2,k_hhhp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_b_b(d_hphh,d_whh,d_v2,k_hphh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_c_a(d_hphh,d_whh,d_v2,k_hphh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_c_b(d_hhhp,d_whh,d_v2,k_hhhp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_d(d_hphp,d_whh,d_v2,k_hphp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_e_a(d_pphh,d_whh,d_v2,k_pphh_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_e_b(d_hhpp,d_whh,d_v2,k_hhpp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_f_a(d_pphp,d_whh,d_v2,k_pphp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_c_f_b(d_hppp,d_whh,d_v2,k_hppp_offset,k_whh_offset,
     1     k_v2_offset)
      call wdm_hh_mo_d(d_whh,d_v2,d_ph,k_whh_offset,k_v2_offset,
     1     k_ph_offset)
      call reconcilefile(d_whh,size_whh)
c     =========================
c     1PWDM Hole-Particle Block
c     =========================
      call tce_filename('whp',filename)
      call tce_dens_hp_offset(l_whp_offset,k_whp_offset,size_whp)
      call createfile(filename,d_whp,size_whp)
      call tce_zero(d_whp,size_whp)
      call wdm_hp_mo_a(d_pp,d_f1,d_whp,k_pp_offset,k_f1_offset,
     1     k_whp_offset)
      call wdm_hp_mo_b_a_a(d_hphh,d_whp,d_v2,k_hphh_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_a_b(d_hhhp,d_whp,d_v2,k_hhhp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_b_a(d_pphh,d_whp,d_v2,k_pphh_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_b_b(d_hhpp,d_whp,d_v2,k_hhpp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_c(d_hphp,d_whp,d_v2,k_hphp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_d_a(d_pphp,d_whp,d_v2,k_pphp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_d_b(d_hppp,d_whp,d_v2,k_hppp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_e_a(d_hppp,d_whp,d_v2,k_hppp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_e_b(d_pphp,d_whp,d_v2,k_pphp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_b_f(d_pppp,d_whp,d_v2,k_pppp_offset,k_whp_offset,
     1     k_v2_offset)
      call wdm_hp_mo_c(d_f1,d_whp,d_ph,k_f1_offset,k_whp_offset,
     1     k_ph_offset)
      call reconcilefile(d_whp,size_whp)
c     =============================
c     1PWDM Particle-Particle Block
c     =============================
      call tce_filename('wpp',filename)
      call tce_dens_pp_offset(l_wpp_offset,k_wpp_offset,size_wpp)
      call createfile(filename,d_wpp,size_wpp)
      call tce_zero(d_wpp,size_wpp)
      call wdm_pp_mo_a(d_pp,d_f1,d_wpp,k_pp_offset,k_f1_offset,
     1     k_wpp_offset)
      call wdm_pp_mo_b_a_a(d_hphh,d_wpp,d_v2,k_hphh_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_a_b(d_hhhp,d_wpp,d_v2,k_hhhp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_b_a(d_pphh,d_wpp,d_v2,k_pphh_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_b_b(d_hhpp,d_wpp,d_v2,k_hhpp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_c(d_hphp,d_wpp,d_v2,k_hphp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_d_a(d_pphp,d_wpp,d_v2,k_pphp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_d_b(d_hppp,d_wpp,d_v2,k_hppp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_e_a(d_hppp,d_wpp,d_v2,k_hppp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_e_b(d_pphp,d_wpp,d_v2,k_pphp_offset,k_wpp_offset,
     1     k_v2_offset)
      call wdm_pp_mo_b_f(d_pppp,d_wpp,d_v2,k_pppp_offset,k_wpp_offset,
     1     k_v2_offset)
      call reconcilefile(d_wpp,size_wpp)
c
c     ===========================================
c     put 2pdm 9 pieces into one file
c     ===========================================
      call tce_filename('2pdm',filename)
      call tce_2pdm_offset(l_2pdm_offset,k_2pdm_offset,size_2pdm)
      call createfile(filename,d_2pdm,size_2pdm)
c
      call put_2pdm_hhhh(d_2pdm,k_2pdm_offset,d_hhhh,k_hhhh_offset)
      call put_2pdm_hhhp(d_2pdm,k_2pdm_offset,d_hhhp,k_hhhp_offset)
      call put_2pdm_hhpp(d_2pdm,k_2pdm_offset,d_hhpp,k_hhpp_offset)
      call put_2pdm_hphh(d_2pdm,k_2pdm_offset,d_hphh,k_hphh_offset)
      call put_2pdm_hphp(d_2pdm,k_2pdm_offset,d_hphp,k_hphp_offset)
      call put_2pdm_hppp(d_2pdm,k_2pdm_offset,d_hppp,k_hppp_offset)
      call put_2pdm_pphh(d_2pdm,k_2pdm_offset,d_pphh,k_pphh_offset)
      call put_2pdm_pphp(d_2pdm,k_2pdm_offset,d_pphp,k_pphp_offset)
      call put_2pdm_pppp(d_2pdm,k_2pdm_offset,d_pppp,k_pppp_offset)
c
      call reconcilefile(d_2pdm,size_2pdm)
c     ===============================
c     put 1pdm 3 pieces into one file
c     ===============================
      call tce_filename('1pdm',filename)
      call tce_1pdm_offset(l_1pdm_offset,k_1pdm_offset,size_1pdm)
      call createfile(filename,d_1pdm,size_1pdm)
c
      call put_1pdm_hh(d_1pdm,k_1pdm_offset,d_hh,k_hh_offset)
      call put_1pdm_ph(d_1pdm,k_1pdm_offset,d_ph,k_ph_offset)
      call put_1pdm_pp(d_1pdm,k_1pdm_offset,d_pp,k_pp_offset)
c
      call reconcilefile(d_1pdm,size_1pdm)
c     ======================================
c     put 1pwdm 3 pieces into one file
c     ======================================
      call tce_filename('1pwdm',filename)
      call tce_1pdm_offset(l_1pwdm_offset,k_1pwdm_offset,size_1pwdm)
      call createfile(filename,d_1pwdm,size_1pwdm)
c
      call put_1pdm_hh(d_1pwdm,k_1pwdm_offset,d_whh,k_whh_offset)
      call put_1pwdm_hp(d_1pwdm,k_1pwdm_offset,d_whp,k_whp_offset)
      call put_1pdm_pp(d_1pwdm,k_1pwdm_offset,d_wpp,k_wpp_offset)
c
      call reconcilefile(d_1pwdm,size_1pwdm)
c
c     =================================
c     HF REFERENCE CONTRIBUTION TO 2PDM
c     =================================
      do g1b=1,noab
      do g2b=1,noab
         if (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),
     1       ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1)))) 
     1       .eq. irrep_e) THEN
         if ((.not.restricted).or.(int_mb(k_spin+g1b-1)+
     1       int_mb(k_spin+g2b-1).ne.4)) then
            dima=int_mb(k_range+g1b-1)*int_mb(k_range+g1b-1)*
     1           int_mb(k_range+g2b-1)*int_mb(k_range+g2b-1)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('tce_cisdt_gradient: ma problem',3,ma_err)
c
            call get_hash_block(d_2pdm,dbl_mb(k_a),dima,
     1           int_mb(k_2pdm_offset),(g2b-1)+(g1b-1)*(noab+nvab)+
     1           (g2b-1)*(noab+nvab)*(noab+nvab)+(g1b-1)*(noab+nvab)*
     1           (noab+nvab)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
            do g2=1,int_mb(k_range+g2b-1)
               dbl_mb(k_a+(g2-1)+(g1-1)*int_mb(k_range+g2b-1)+
     1         (g2-1)*int_mb(k_range+g2b-1)*int_mb(k_range+g1b-1)+
     1         (g1-1)*int_mb(k_range+g2b-1)*int_mb(k_range+g1b-1)*
     1                int_mb(k_range+g2b-1)) = 
     1         dbl_mb(k_a+(g2-1)+(g1-1)*int_mb(k_range+g2b-1)+
     1         (g2-1)*int_mb(k_range+g2b-1)*int_mb(k_range+g1b-1)+
     1         (g1-1)*int_mb(k_range+g2b-1)*int_mb(k_range+g1b-1)*
     1                int_mb(k_range+g2b-1)) + 1.d0
            enddo
            enddo
            call put_hash_block(d_2pdm,dbl_mb(k_a),dima,
     1           int_mb(k_2pdm_offset),(g2b-1)+(g1b-1)*(noab+nvab)+
     1           (g2b-1)*(noab+nvab)*(noab+nvab)+(g1b-1)*(noab+nvab)*
     1           (noab+nvab)*(noab+nvab))
c
            call get_hash_block(d_2pdm,dbl_mb(k_a),dima,
     1           int_mb(k_2pdm_offset),(g1b-1)+(g2b-1)*(noab+nvab)+
     1           (g2b-1)*(noab+nvab)*(noab+nvab)+(g1b-1)*(noab+nvab)*
     1           (noab+nvab)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
            do g2=1,int_mb(k_range+g2b-1)
               dbl_mb(k_a+(g1-1)+(g2-1)*int_mb(k_range+g1b-1)+
     1         (g2-1)*int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)+
     1         (g1-1)*int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1                int_mb(k_range+g2b-1)) =
     1         dbl_mb(k_a+(g1-1)+(g2-1)*int_mb(k_range+g1b-1)+
     1         (g2-1)*int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)+
     1         (g1-1)*int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1                int_mb(k_range+g2b-1)) - 1.d0
            enddo
            enddo
            call put_hash_block(d_2pdm,dbl_mb(k_a),dima,
     1           int_mb(k_2pdm_offset),(g1b-1)+(g2b-1)*(noab+nvab)+
     1           (g2b-1)*(noab+nvab)*(noab+nvab)+(g1b-1)*(noab+nvab)*
     1           (noab+nvab)*(noab+nvab))
            if (.not.ma_pop_stack(l_a))
     1      call errquit('tce_cisdt_gradient: ma problem',4,ma_err)
c
         endif
         endif
      enddo
      enddo
c    =================================
c    HF REFERENCE CONTRIBUTION TO 1PDM
c    =================================
      do g1b=1,noab
         if ((.not.restricted).or.(int_mb(k_spin+g1b-1).ne.2)) then
            dima=int_mb(k_range+g1b-1)*int_mb(k_range+g1b-1)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('tce_cisdt_gradient: ma problem',2,ma_err)
            call get_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))+1.d0
            enddo
            call put_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            if (.not.ma_pop_stack(l_a))
     1      call errquit('tce_cisdt_gradient: ma problem',3,ma_err)
         endif
      enddo
c     ==================================
c     HF REFERENCE CONTRIBUTION TO 1PWDM
c     ==================================
      do g1b=1,noab
         if ((.not.restricted).or.(int_mb(k_spin+g1b-1).ne.2)) then
            dima = int_mb(k_range+g1b-1)*int_mb(k_range+g1b-1)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('tce_cisdt_gradient: ma problem',4,ma_err)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_b,k_b))
     1      call errquit('tce_cisdt_gradient: ma problem',5,ma_err)
c
            call get_hash_block(d_1pwdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pwdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            call get_hash_block(d_f1,dbl_mb(k_b),dima,
     1           int_mb(k_f1_offset),(g1b-1)+(g1b-1)*(noab+nvab))
c
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1)) -
     1         dbl_mb(k_b+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))
            enddo
            call put_hash_block(d_1pwdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pwdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            if (.not.ma_pop_stack(l_b))
     1      call errquit('tce_cisdt_gradient: ma problem',6,ma_err)
            if (.not.ma_pop_stack(l_a))
     1      call errquit('tce_cisdt_gradient: ma problem',3,ma_err)
         endif
      enddo
c     ==============================
c
      call int_terminate()
      call schwarz_tidy()
c     =================
c     Do the gradients
c     =================
      status=tce_gradients_doit(rtdb, model,
     1       d_1pdm,  k_1pdm_offset,  size_1pdm,
     1       d_1pwdm, k_1pwdm_offset, size_1pwdm,
     1       d_2pdm,  k_2pdm_offset,  size_2pdm,
     1       noab,nvab,k_range,k_offset,k_spin,
     1       k_movecs_sorted,restricted)
c     =============================
c     clean up the files and arrays
c     =============================
      call deletefile(d_1pwdm)
      if(.not.ma_pop_stack(l_1pwdm_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',997,ma_err)
      call deletefile(d_1pdm)
      if(.not.ma_pop_stack(l_1pdm_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',996,ma_err)
      call deletefile(d_2pdm)
      if(.not.ma_pop_stack(l_2pdm_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',995,ma_err)     
c
      call deletefile(d_wpp)
      if(.not.ma_pop_stack(l_wpp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',1,ma_err)
c
      call deletefile(d_whp)
      if(.not.ma_pop_stack(l_whp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',3,ma_err)
c
      call deletefile(d_whh)
      if(.not.ma_pop_stack(l_whh_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',4,ma_err)
c
      call deletefile(d_ph)
      if(.not.ma_pop_stack(l_ph_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',6,ma_err)
c
      call deletefile(d_x1)
      if(.not.ma_pop_stack(l_x1_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',7,ma_err)
c
      call deletefile(d_pppp)
      if(.not.ma_pop_stack(l_pppp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',8,ma_err)
c
      call deletefile(d_pphp)
      if(.not.ma_pop_stack(l_pphp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',9,ma_err)
c
      call deletefile(d_pphh)
      if(.not.ma_pop_stack(l_pphh_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',10,ma_err)
c
      call deletefile(d_hppp)
      if(.not.ma_pop_stack(l_hppp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',11,ma_err)
c
      call deletefile(d_hphp)
      if(.not.ma_pop_stack(l_hphp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',12,ma_err)
c
      call deletefile(d_hphh)
      if(.not.ma_pop_stack(l_hphh_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',13,ma_err)
c
      call deletefile(d_hhpp)
      if(.not.ma_pop_stack(l_hhpp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',14,ma_err)
c
      call deletefile(d_hhhp)
      if(.not.ma_pop_stack(l_hhhp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',15,ma_err)
c
      call deletefile(d_hhhh)
      if(.not.ma_pop_stack(l_hhhh_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',16,ma_err)
c
      call deletefile(d_pp)
      if(.not.ma_pop_stack(l_pp_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',17,ma_err)
c
      call deletefile(d_hh)
      if(.not.ma_pop_stack(l_hh_offset))
     1  call errquit('tce_cisdt_gradient: ma problem',18,ma_err)
c
      call deletefile(d_d0)
      if (.not.ma_pop_stack(l_d0_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',19,ma_err)
c
      call deletefile(d_y3)
      if (.not.ma_pop_stack(l_y3_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',20,ma_err)
c
      call deletefile(d_y2)
      if (.not.ma_pop_stack(l_y2_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',21,ma_err)
c
      call deletefile(d_y1)
      if (.not.ma_pop_stack(l_y1_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',22,ma_err)
c
      call deletefile(d_y0)
      if (.not.ma_pop_stack(l_y0_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',23,ma_err)
c
      call deletefile(d_t0)
      if (.not. ma_pop_stack(l_t0_offset))
     1   call errquit('tce_cisdt_gradient: ma problem',24,ma_err)
c
      end
c $Id$
