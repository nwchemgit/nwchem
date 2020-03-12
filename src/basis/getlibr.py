#!/usr/bin/python3
# This script downloads the basis set library data from www.basissetexchange.org
# into the directory $NWCHEM_TOP/src/basis/libraries.bse
# to use, set the env. variable NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries.bse/
# Requires the installation of the python env. from
# https://github.com/MolSSI-BSE/basis_set_exchange
# See https://molssi-bse.github.io/basis_set_exchange/
#
# names changed
# def2-universal-jfit was weigend_coulomb_fitting
# dgauss-a1-dftjfit   was dgauss_a1_dft_coulomb_fitting
# dgauss-a2-dftjfit   was dgauss_a2_dft_coulomb_fitting
#
import basis_set_exchange as bse
from datetime import datetime
today = datetime.now().isoformat(timespec='minutes')
print(today)
all_bs = bse.get_all_basis_names()
md = bse.get_metadata()
for bas_name in all_bs:
    #get version and list of elements
    version_bs = md[bas_name]['latest_version']
    elements_list = md[bas_name]['versions'][version_bs]['elements']
    #open file
    # get rid of asterisks
    file_name = bas_name.replace("*","s")
    #get rid of parenthesis
    file_name = file_name.replace("(","")
    file_name = file_name.replace(")","")
    #replace commas with underscore
    file_name = file_name.replace(",","_")
    #replace whitespace with underscore
    file_name = file_name.replace(" ","_")
    #replace forward slash with underscore
    file_name = file_name.replace("/","_")
    print(' file name is '+file_name+"\n")
    output_file = open(file_name,'w')
    output_file.write('# BSE Version '+bse.version()+'\n')
    output_file.write('# Data downloaded at '+today+'\n')
    output_file.write('# '+bas_name+' version number '+version_bs+'\n')
    output_file.write('# Description: '+md[bas_name]['description']+'\n')
    output_file.write('# Role: '+md[bas_name]['role']+'\n')
    output_file.write('# '+bse.get_references(bas_name,fmt='txt').replace('\n','\n# '))
    output_file.write('# \n')
    for element in elements_list:
        #element='h'
        try:
            bs_str=bse.get_basis(bas_name, header=False, elements=element, fmt='nwchem', optimize_general=True, uncontract_general=True)
        except:
#            print("failed for"+element)
            pass
        else:
            bs_str=bs_str.replace("BASIS","basis")
            bs_str=bs_str.replace("END","end")
            bs_str=bs_str.replace("PRINT","")
            element_str=bse.misc.compact_elements([element])
            bs_str=bs_str.replace("ao basis",element_str+"_"+bas_name)
            #ECP
            bs_str=bs_str.replace("ECP","ecp \""+element_str+"_"+bas_name+"\"")
            output_file.write(bs_str)
            #
            print(bas_name+" "+element_str)
print("end")

