#!/bin/sh
#
cd /tmp
runtest.unix procs 2 \
 $NWCHEM_TEST/na_k/nak \
 $NWCHEM_TEST/na_k/nak_md \
 $NWCHEM_TEST/crown/crown_md \
 $NWCHEM_TEST/ethanol/ethanol_md \
 $NWCHEM_TEST/ethanol/ethanol_ti \
 $NWCHEM_TEST/had/had_em \
 $NWCHEM_TEST/had/had_md \
 $NWCHEM_TEST/prep/a3n \
 $NWCHEM_TEST/prep/aal \
 $NWCHEM_TEST/prep/fsc \
 $NWCHEM_TEST/water/water_md
# $NWCHEM_TEST/h2o_qmd/h2o_scf 
