      logical function pre_mktop(lfnout,ffield,
     + lfnseq,filseq,lfntop,filtop,lfnsgm,lfnpar,filpar,
     + lfnmod,filmod,lfnmat,dir_s,dir_x,dir_u,dir_t,ignore,
     + slvnam)
c
c $Id: pre_mktop.F,v 1.26 2001-07-06 03:45:13 d3j191 Exp $
c
c     in : integer lfnout    = logical file number output file
c          char*80 ffield    = force field from [amber]
c          integer lfnseq    = logical file number for seq file
c          char*80 filseq    = file name of seq file
c          integer lfntop    = logical file number for top file
c          char*80 filtop    = file name of top file
c          integer lfnsgm    = logical file number for sgm file
c          char*80 dir_s     = directory standard ffield files
c          char*80 dir_x     = directory extensions ffield files
c          char*80 dir_u     = directory user preferences ffield files
c          char*80 dir_t     = directory temporary ffield files
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
c
      logical pre_seqsiz,pre_rdseq,pre_term,pre_dimens,pre_mklist
      external pre_seqsiz,pre_rdseq,pre_term,pre_dimens,pre_mklist
      logical pre_params,pre_natyps,pre_wrttop,pre_third,pre_excl
      external pre_params,pre_natyps,pre_wrttop,pre_third,pre_excl
      logical pre_nonbon,pre_modify
      external pre_nonbon,pre_modify
c
      integer lfnout,lfnseq,lfntop,lfnsgm,lfnpar,lfnmod,ignore,lfnmat
      character*80 ffield
      character*255 filseq,filtop,filpar,filmod
      character*255 dir_s,dir_x,dir_u,dir_t
      character*3 slvnam
c
      integer mseq,nseq
      integer l_lseq,i_lseq,l_cseq,i_cseq
c
      integer mlnk,nlnk
      integer l_llnk,i_llnk,l_clnk,i_clnk
c
      integer matm,natm
      integer l_latm,i_latm,l_catm,i_catm,l_qatm,i_qatm
c
      integer mbnd,nbnd
      integer l_lbnd,i_lbnd,l_rbnd,i_rbnd
c
      integer mang,nang
      integer l_lang,i_lang,l_rang,i_rang
c
      integer mdih,ndih
      integer l_ldih,i_ldih,l_kdih,i_kdih,l_rdih,i_rdih
c
      integer mimp,nimp
      integer l_limp,i_limp,l_rimp,i_rimp
c
      integer matt,natt,nats
      integer l_latt,i_latt,l_catt,i_catt,l_patt,i_patt,l_ratt,i_ratt
      integer l_lats,i_lats
c
      integer nval
      integer l_ival,i_ival,l_rval,i_rval,l_ndx,i_ndx
c
      integer m3rd,n3rd
      integer i_l3rd,l_l3rd
c
      integer mexc,nexc
      integer i_lexc,l_lexc
c
      integer mnon,nnon
      integer i_lnon,l_lnon
c
      integer natmt,nbndt,nangt,ndiht,nimpt,n3rdt,mqu
c
      real*8 releps,q14fac,wcorr(3)
c
      real*8 timer_wall_total
      external timer_wall_total
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1000)
 1000 format('TOPOLOGY GENERATION')
      endif
c
      call timer_init()
      if(.not.pre_seqsiz(lfnout,lfnseq,filseq,nseq,nlnk))
     + call errquit('pre_seqsize failed',9999)
c
      mseq=nseq+2
      mlnk=nlnk+1
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1001) nseq
 1001 format('sequence length is',i10)
      endif
c
c     allocate memory for sequence
c     ----------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : link type
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
      if(.not.ma_push_get(mt_int,4*mseq,'lseq',l_lseq,i_lseq))
     + call errquit('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
      if(.not.ma_push_get(mt_byte,10*mseq,'cseq',l_cseq,i_cseq))
     + call errquit('Memory allocation failed for cseq',9999)
c
c     allocate memory for link list
c     -----------------------------
c
c     integer llnk(1,mlnk) : link segment 1
c                  2       : link segment 2
c                  3       : link type : 0: normal
c                                        1: forced
c     char*4  clnk(1,mlnk) : link segment 1 atom name
c                  2       : link segment 2 atom name
c
      if(.not.ma_push_get(mt_int,3*mlnk,'llnk',l_llnk,i_llnk))
     + call errquit('Memory allocation failed for llnk',9999)
      if(.not.ma_push_get(mt_byte,8*mlnk,'clnk',l_clnk,i_clnk))
     + call errquit('Memory allocation failed for clnk',9999)
c
c     read the sequence file
c     ----------------------
c
      if(.not.pre_rdseq(lfnout,lfnseq,filseq,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call errquit('pre_rdseq failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1002)
 1002 format('sequence is read')
      endif
c
c     determine termini
c     -----------------
c
      if(.not.pre_term(lfnout,lfnsgm,dir_s,dir_x,dir_u,dir_t,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call errquit('pre_term failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1003)
 1003 format('termini are determined')
      endif
c
c     determine array dimensions
c     --------------------------
c
      if(.not.pre_dimens(lfnout,lfnsgm,dir_s,dir_x,dir_u,dir_t,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + natm,nbnd,nang,ndih,nimp))
     + call errquit('pre_dimens failed',9999)
c
c     adjust dimension to account for links
c     -------------------------------------
c
      matm=natm+1
      mbnd=nbnd+3*nseq+1
      mang=nang+27*nseq+1
      mdih=ndih+45*nseq+1
      mimp=nimp+6*nseq+1
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1004) mseq,matm,mbnd,mang,mdih,mimp
 1004 format('list dimensions are determined:',/,
     + ' number of segments  ',i10,/,
     + ' number of atoms     ',i10,/,
     + ' number of bonds     ',i10,/,
     + ' number of angles    ',i10,/,
     + ' number of dihedrals ',i10,/,
     + ' number of impropers ',i10)
      endif
c
c     allocate memory for atomic lists
c     --------------------------------
c
c     integer latm(1,matm) : charge group
c                  2       : polarization group
c                  3       : link number
c                  4       : center type
c                  5       : segment number
c                  6       : molecule number
c                  7       : atomic number
c                  8       : fraction
c
c
c     char*6  catm(1,matm) : atom name 
c                  2       : atom type set 1
c                  3       : atom type set 2
c                  4       : atom type set 3
c
      if(.not.ma_push_get(mt_int,9*matm,'latm',l_latm,i_latm))
     + call errquit('Memory allocation failed for latm',9999)
      if(.not.ma_push_get(mt_byte,24*matm,'catm',l_catm,i_catm))
     + call errquit('Memory allocation failed for catm',9999)
c
c     real*8 qatm(1,matm) : partial atomic charges set 1
c                 2       : polarizability set 1
c                 3       : partial atomic charges set 2
c                 4       : polarizability set 2
c                 5       : partial atomic charges set 3
c                 6       : polarizability set 3
c
      mqu=3*matm
      if(.not.ma_push_get(mt_dbl,6*mqu,'qatm',l_qatm,i_qatm))
     + call errquit('Memory allocation failed for qatm',9999)
c
c     allocate memory for bonded lists
c     --------------------------------
c
c     integer lbnd(1:2,mbnd) : bond indices
c                  3         : parameter source
c                  4         : constraint type
c     real*8  rbnd(1,mbnd)   : bond length set 1
c                  2         : force constant set 1
c                  3         : bond length set 2
c                  4         : force constant set 2
c                  5         : bond length set 3
c                  6         : force constant set 3
c
      if(.not.ma_push_get(mt_int,4*mbnd,'lbnd',l_lbnd,i_lbnd))
     + call errquit('Memory allocation failed for lbnd',9999)
      if(.not.ma_push_get(mt_dbl,6*mbnd,'rbnd',l_rbnd,i_rbnd))
     + call errquit('Memory allocation failed for rbnd',9999)
c
c     integer lang(1:3,mang) : angle indices
c                  4         : parameter source
c                  5         : constraint type
c     real*8  rang(1,mang)   : angle set 1
c                  2         : force constant set 1
c                  3         : angle set 2
c                  4         : force constant set 2
c                  5         : angle set 3
c                  6         : force constant set 3
c
      if(.not.ma_push_get(mt_int,5*mang,'lang',l_lang,i_lang))
     + call errquit('Memory allocation failed for lang',9999)
      if(.not.ma_push_get(mt_dbl,6*mang,'rang',l_rang,i_rang))
     + call errquit('Memory allocation failed for rang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c                  5         : parameter source
c                  6         : constraint type
c                  7         : number of functions set 1
c                  8         : number of functions set 2
c                  9         : number of functions set 3
c
c     integer kdih(6,1,mimp) : multiplicity set 1
c                    2       : multiplicity set 2
c                    3       : multiplicity set 3
c
c     real*8  rdih(6,1,mdih) : phase angle set 1
c                    2       : force constant set 1
c                    3       : phase angle set 2
c                    4       : force constant set 2
c                    5       : pahse angle set 3
c                    6       : force constant set 3 
c
      if(.not.ma_push_get(mt_int,9*mdih,'ldih',l_ldih,i_ldih))
     + call errquit('Memory allocation failed for ldih',9999)
      if(.not.ma_push_get(mt_int,18*mdih,'kdih',l_kdih,i_kdih))
     + call errquit('Memory allocation failed for ldih',9999)
      if(.not.ma_push_get(mt_dbl,36*mdih,'rdih',l_rdih,i_rdih))
     + call errquit('Memory allocation failed for rdih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c                  5         : parameter source
c                  6         : constraint type
c                  7         : multiplicity set 1
c                  8         : multiplicity set 2
c                  9         : multiplicity set 3
c
c     real*8  rimp(1,mimp)   : phase angle set 1
c                  2         : force constant set 1
c                  3         : phase angle set 2
c                  4         : force constant set 2
c                  5         : pahse angle set 3
c                  6         : force constant set 3
c
      if(.not.ma_push_get(mt_int,9*mimp,'limp',l_limp,i_limp))
     + call errquit('Memory allocation failed for limp',9999)
      if(.not.ma_push_get(mt_dbl,6*mimp,'rimp',l_rimp,i_rimp))
     + call errquit('Memory allocation failed for rimp',9999)
c
c     generate the atomic and bonded lists
c     ------------------------------------
c
      natm=0
      nbnd=0
      nang=0
      ndih=0
      nimp=0
c
      nval=max(mbnd,mang,mdih,mimp,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call errquit('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call errquit('Memory allocation failed for ndx',9999)
      if(.not.ma_push_get(mt_dbl,nval,'rval',l_rval,i_rval))
     + call errquit('Memory allocation failed for rval',9999)
      if(.not.pre_mklist(lfnout,lfnsgm,lfnmat,
     + ffield,dir_s,dir_x,dir_u,dir_t,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,
     + int_mb(i_ldih),int_mb(i_kdih),dbl_mb(i_rdih),mdih,ndih,
     + int_mb(i_limp),dbl_mb(i_rimp),mimp,nimp,
     + natmt,nbndt,nangt,ndiht,nimpt,wcorr,
     + nval,int_mb(i_ival),dbl_mb(i_rval),int_mb(i_ndx)))
     + call errquit('pre_mklist failed',9999)
      if(.not.ma_pop_stack(l_rval))
     + call errquit('Memory deallocation failed for rval',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call errquit('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call errquit('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1005)
 1005 format('lists are generated')
      endif
c
c     apply topology modifications to atom types
c     ------------------------------------------
c
      if(.not.pre_modify(1,lfnout,lfnmod,filmod,ffield,int_mb(i_latm),
     + byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,nangt,
     + int_mb(i_ldih),int_mb(i_kdih),dbl_mb(i_rdih),mdih,ndih,ndiht,
     + int_mb(i_limp),dbl_mb(i_rimp),mimp,nimp,nimpt))
     + call errquit('pre_modify failed',9999)
c
c     get number of atom types
c     ------------------------
c
      if(.not.pre_natyps(lfnout,byte_mb(i_catm),matm,natmt,natt))
     + call errquit('pre_natyps failed',9999)
      matt=natt
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1006) natt
 1006 format('number of atom types is ',i10)
      endif
c
c     allocate memory for atom type lists
c     -----------------------------------
c
      if(.not.ma_push_get(mt_int,3*matt,'latt',l_latt,i_latt))
     + call errquit('Memory allocation failed for latt',9999)
      if(.not.ma_push_get(mt_int,3*matt,'lats',l_lats,i_lats))
     + call errquit('Memory allocation failed for lats',9999)
      if(.not.ma_push_get(mt_byte,36*matt,'catt',l_catt,i_catt))
     + call errquit('Memory allocation failed for catt',9999)
      if(.not.ma_push_get(mt_dbl,12*matt*matt,'patt',l_patt,i_patt))
     + call errquit('Memory allocation failed for patt',9999)
      if(.not.ma_push_get(mt_dbl,3*matt,'ratt',l_ratt,i_ratt))
     + call errquit('Memory allocation failed for ratt',9999)
c
c     substitute force field parameters
c     ---------------------------------
c
      if(.not.pre_params(lfnpar,filpar,lfnout,
     + dir_s,dir_x,dir_u,dir_t,ffield,releps,q14fac,ignore,
     + int_mb(i_latm),byte_mb(i_catm),matm,natmt,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nangt,
     + int_mb(i_ldih),int_mb(i_kdih),dbl_mb(i_rdih),mdih,ndiht,
     + int_mb(i_limp),dbl_mb(i_rimp),mimp,nimpt,
     + int_mb(i_latt),int_mb(i_lats),byte_mb(i_catt),dbl_mb(i_patt),
     + dbl_mb(i_ratt),matt,natt,nats))
     + call errquit('pre_params failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1007)
 1007 format('force field parameters are substituted')
      endif
c
c     allocate memory for third neighbor list
c     ---------------------------------------
c
      m3rd=6*mang
      n3rd=0
      if(.not.ma_push_get(mt_int,2*m3rd,'l3rd',l_l3rd,i_l3rd))
     + call errquit('Memory allocation failed for l3rd',9999)
      nval=max(m3rd,matm,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call errquit('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call errquit('Memory allocation failed for ndx',9999)
c
c     construct third neighbor list
c     -----------------------------
c
      if(.not.pre_third(int_mb(i_lbnd),mbnd,nbnd,int_mb(i_lang),
     + mang,nang,int_mb(i_l3rd),m3rd,n3rd,1,1,1,
     + nval,int_mb(i_ival),int_mb(i_ndx),matm))
     + call errquit('pre_third failed',9999)
c
      n3rdt=n3rd
      if(nbndt.gt.nbnd.and.nangt.gt.nang) then
      if(.not.pre_third(int_mb(i_lbnd),mbnd,nbndt,int_mb(i_lang),
     + mang,nangt,int_mb(i_l3rd),m3rd,n3rdt,nbnd+1,nang+1,n3rd+1,
     + nval,int_mb(i_ival),int_mb(i_ndx),matm))
     + call errquit('pre_third failed',9999)
      endif
      if(.not.ma_pop_stack(l_ndx))
     + call errquit('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call errquit('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1008) n3rd
 1008 format(' Solute third neighbor list length is',i10)
      endif
c
c     allocate memory for excluded pair list
c     --------------------------------------
c
      mexc=n3rd+nbnd+nang
      nexc=0
      if(.not.ma_push_get(mt_int,2*mexc,'lexc',l_lexc,i_lexc))
     + call errquit('Memory allocation failed for lexc',9999)
c
c     construct excluded pair list
c     ----------------------------
c
      nval=max(mexc,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call errquit('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call errquit('Memory allocation failed for ndx',9999)
      if(.not.pre_excl(int_mb(i_l3rd),m3rd,n3rd,int_mb(i_lbnd),
     + mbnd,nbnd,int_mb(i_lang),mang,nang,int_mb(i_lexc),mexc,nexc,
     + nval,int_mb(i_ival),int_mb(i_ndx)))
     + call errquit('pre_excl failed',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call errquit('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call errquit('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1009) nexc
 1009 format(' Solute excluded pair list length is',i10)
      endif
c
c     allocate memory for solvent non-bonded list
c     -------------------------------------------
c
      mnon=(natmt-natm)*(natmt-natm)
      nnon=0
      if(.not.ma_push_get(mt_int,mnon,'lnon',l_lnon,i_lnon))
     + call errquit('Memory allocation failed for lnon',9999)
c
c     construct solvent non-bonded list
c     ---------------------------------
c
      if(.not.pre_nonbon(natm+1,natmt,int_mb(i_lbnd),mbnd,nbnd+1,nbndt,
     + int_mb(i_lang),mang,nang+1,nangt,
     + int_mb(i_l3rd),m3rd,n3rd+1,n3rdt,
     + int_mb(i_lnon),mnon,nnon))
     + call errquit('pre_nonbon failed',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1010) nnon
 1010 format(' Solvent non-bonded list length is',i10)
      endif
c
c     apply topology modifications to bonded parameters
c     -------------------------------------------------
c
      if(.not.pre_modify(2,lfnout,lfnmod,filmod,ffield,int_mb(i_latm),
     + byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,nangt,
     + int_mb(i_ldih),int_mb(i_kdih),dbl_mb(i_rdih),mdih,ndih,ndiht,
     + int_mb(i_limp),dbl_mb(i_rimp),mimp,nimp,nimpt))
     + call errquit('pre_modify failed',9999)
c
c     write the topology file
c     -----------------------
c
      if(.not.pre_wrttop(lfnout,ffield,lfntop,filtop,releps,q14fac,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,mqu,
     + int_mb(i_latt),int_mb(i_lats),byte_mb(i_catt),dbl_mb(i_patt),
     + dbl_mb(i_ratt),matt,natt,nats,int_mb(i_latm),byte_mb(i_catm),
     + dbl_mb(i_qatm),matm,natm,int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,
     + nbnd,int_mb(i_lang),
     + dbl_mb(i_rang),mang,nang,int_mb(i_ldih),int_mb(i_kdih),
     + dbl_mb(i_rdih),mdih,ndih,
     + int_mb(i_limp),dbl_mb(i_rimp),mimp,nimp,int_mb(i_l3rd),
     + m3rd,n3rd,int_mb(i_lexc),mexc,nexc,int_mb(i_lnon),mnon,nnon,
     + natmt,nbndt,nangt,ndiht,nimpt,n3rdt,wcorr,slvnam))
     + call errquit('pre_wrttop failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1011)
 1011 format(' Topology file is written')
      endif
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_lnon))
     + call errquit('Memory deallocation failed for lnon',9999)
      if(.not.ma_pop_stack(l_lexc))
     + call errquit('Memory deallocation failed for lexc',9999)
      if(.not.ma_pop_stack(l_l3rd))
     + call errquit('Memory deallocation failed for l3rd',9999)
      if(.not.ma_pop_stack(l_ratt))
     + call errquit('Memory deallocation failed for ratt',9999)
      if(.not.ma_pop_stack(l_patt))
     + call errquit('Memory deallocation failed for patt',9999)
      if(.not.ma_pop_stack(l_catt))
     + call errquit('Memory deallocation failed for catt',9999)
      if(.not.ma_pop_stack(l_lats))
     + call errquit('Memory deallocation failed for lats',9999)
      if(.not.ma_pop_stack(l_latt))
     + call errquit('Memory deallocation failed for latt',9999)
      if(.not.ma_pop_stack(l_rimp))
     + call errquit('Memory deallocation failed for rimp',9999)
      if(.not.ma_pop_stack(l_limp))
     + call errquit('Memory deallocation failed for limp',9999)
      if(.not.ma_pop_stack(l_rdih))
     + call errquit('Memory deallocation failed for rdih',9999)
      if(.not.ma_pop_stack(l_kdih))
     + call errquit('Memory deallocation failed for kdih',9999)
      if(.not.ma_pop_stack(l_ldih))
     + call errquit('Memory deallocation failed for ldih',9999)
      if(.not.ma_pop_stack(l_rang))
     + call errquit('Memory deallocation failed for rang',9999)
      if(.not.ma_pop_stack(l_lang))
     + call errquit('Memory deallocation failed for lang',9999)
      if(.not.ma_pop_stack(l_rbnd))
     + call errquit('Memory deallocation failed for rbnd',9999)
      if(.not.ma_pop_stack(l_lbnd))
     + call errquit('Memory deallocation failed for lbnd',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call errquit('Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call errquit('Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_latm))
     + call errquit('Memory deallocation failed for latm',9999)
      if(.not.ma_pop_stack(l_clnk))
     + call errquit('Memory deallocation failed for clnk',9999)
      if(.not.ma_pop_stack(l_llnk))
     + call errquit('Memory deallocation failed for llnk',9999)
      if(.not.ma_pop_stack(l_cseq))
     + call errquit('Memory deallocation failed for cseq',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call errquit('Memory deallocation failed for lseq',9999)
c
      pre_mktop=.true.
      return
      end
