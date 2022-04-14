# $Id$


  OBJ = atim.o cio.o  lucia.o

  OBJ_OPTIMIZE = aixopn.o aixrd.o baserou.o cdafile.o crop_simple.o \
	daname.o daname_main.o ddafile.o densi2_ab.o determinant.o \
	diis.o diis_simple.o dim_pa.o dummies.o \
	eispack.o ekt.o filehand.o gasdir2.o gasrou.o gucc_fock.o i2addr.o \
	idafile_.o idafile.o linpack1.o lucia1.o lucia3.o lucia4.o \
	lucia5.o lucia6.o lucia7.o lucia_akbk.o lucia_anders.o lucia_caspt.o \
	lucia_cc.o lucia_ccrsp.o lucia_csf.o lucia_gencc2.o lucia_gencc3.o \
	lucia_gencc4.o lucia_gencc5.o lucia_gencc.o lucia_gtbce.o lucia_hf.o \
	lucia_ic2.o lucia_icci0.o lucia_icei.o lucia_ic.o lucia_mcscf.o \
	lucia_nortci.o lucia_occ.o lucia_oscc.o lucia_prodexp.o lucia_supsym.o \
	luciaw.o lucifc.o matml7.o micdv6o.o numint.o optcont.o optcont_lambda.o \
	optcont_orbopt.o optc_routines.o prtim.o qenter.o \
	strnln.o tra_kappa.o tra_simple.o triprt.o \
	ucc_grad.o vg.o wr_motra_info.o  \
        rsbb2ben.o rsbb2b_m.o rsbb2bn_ms.o
#	ucc_grad.o vg.o wr_motra_info.o  lucia.o \
#       rsbb2ben.o rsbb2b_m.o rsbb2bn_ms.o

      LIBRARY =	liblucia.a

    USES_BLAS = baserou.F gucc_fock.F lucia5.F lucia.F linpack1.F lucia_icei.F matml7.F

HEADERS = akbk.inc cands.inc cc_exc.inc ceblk.inc cecore.inc cei.inc cfinal_e.inc \
	cgas.inc chktyp.inc cicisp.inc cintfo.inc clunit.inc comjep.inc constants.inc \
	corbex.inc cprnt.inc crun.inc csfbas.inc cshift.inc csm.inc \
	csmprd.inc cstate.inc ctcc.inc ctl.inc cxpt.inc dbgsym.inc extbas.inc fio.inc \
	fragmol.inc freeunits.inc frorbs.inc gasstr.inc gicci.inc glbbas.inc glstate.inc \
	gtbce.inc implicit.inc intform.inc ioparm.inc irat.inc kprdwvf.inc lorr.inc \
	lucinp.inc machine.inc multd2h.inc mv7task.inc MxFile.inc mxpdim.inc newccp.inc \
	oper.inc opti.inc orbinp.inc prdwvf.inc priunit.inc proto_sp_mat.inc rou_stat.inc \
	spinfo.inc stinf.inc strbas.inc strinp.inc switch.inc symrdc.inc \
	vb.inc wrkspc.inc wrkspc-static.inc

include ../config/makefile.h
include ../config/makelib.h
