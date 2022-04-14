#!/usr/bin/python3
# This script downloads the basis set library data from www.basissetexchange.org
# into the directory $NWCHEM_TOP/src/basis/libraries.bse
# To run, cd $NWCHEM_TOP/src/basis/libraries.bse/ && ../getlibr.py
# this will update the content of $NWCHEM_TOP/src/basis/libraries.bse
# To use the updates library, set the env. variable NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries.bse/
# Requires the installation of the python env. from
# https://github.com/MolSSI-BSE/basis_set_exchange
# e.g. python3 -m pip install --user basis_set_exchange
# See https://molssi-bse.github.io/basis_set_exchange/
#
# names changed
# def2-universal-jfit was weigend_coulomb_fitting
# dgauss-a1-dftjfit   was dgauss_a1_dft_coulomb_fitting
# dgauss-a2-dftjfit   was dgauss_a2_dft_coulomb_fitting

import basis_set_exchange as bse
from datetime import datetime
today = datetime.now().isoformat(timespec='minutes')
print(today)
all_bs = bse.get_all_basis_names()
md = bse.get_metadata()
summary_file = open('summary.txt','w')

def writebs(md, bas_name, summary_file, get_aux=0):
    md_bas_name = bas_name.lower()
    md_bas_name = md_bas_name.replace("*","_st_")
    md_bas_name = md_bas_name.replace("/","_sl_")
    print(' md_bas_name '+md_bas_name+"\n")
    print(' bas_name '+bas_name+"\n")
    version_bs = md[md_bas_name]['latest_version']
    elements_list = md[md_bas_name]['versions'][version_bs]['elements']
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
    #lowercase
    file_name = file_name.lower()
    if get_aux==1:
        file_name = file_name + "-autoaux"
    print(' file name is '+file_name+"\n")
    output_file = open(file_name,'w')
    output_file.write('# BSE Version '+bse.version()+'\n')
    output_file.write('# Data downloaded on '+today+'\n')
    
    if get_aux==0:
        output_file.write('# '+bas_name+' version number '+version_bs+'\n')
        output_file.write('# Description: '+md[md_bas_name]['description']+'\n')
        output_file.write('# Role: '+md[md_bas_name]['role']+'\n')
        output_file.write('# '+bse.get_references(bas_name,fmt='txt').replace('\n','\n# '))
        output_file.write('# \n')
    elif get_aux==1:
        output_file.write('# '+bas_name+' version number '+version_bs+' AutoAux \n')
        output_file.write('# Role: JK Fitting \n')
        output_file.write('# Stoychev GL, Auer AA, Neese F. \n# Automatic Generation of Auxiliary Basis Sets.\n# J Chem Theory Comput. 2017 Feb 14;13(2):554-562.\n# doi: 10.1021/acs.jctc.6b01041.\n')
        output_file.write('# \n')
        
    n_elements=0
    for element in elements_list:
        n_elements = n_elements + 1
    if get_aux==1:
        summary_file.write('Basis set \"'+bas_name+'-autoaux\" (number of atoms '+str(n_elements)+')\n')
    else:
        summary_file.write('Basis set \"'+bas_name+'\" (number of atoms '+str(n_elements)+')\n')
    for element in elements_list:
        #element='h'
        try:
            bs_str=bse.get_basis(bas_name, header=False, elements=element, fmt='nwchem', optimize_general=True, uncontract_general=True, get_aux=get_aux)
        except:
#            print("failed for"+element)
            pass
        else:
            bs_str=bs_str.replace("BASIS","basis")
            bs_str=bs_str.replace("END","end")
            bs_str=bs_str.replace("PRINT","")
            element_str=bse.misc.compact_elements([element])
            if get_aux==1:
                bs_str=bs_str.replace("ao basis",element_str+"_"+bas_name+"-autoaux")
            else:
                bs_str=bs_str.replace("ao basis",element_str+"_"+bas_name)
            #ECP
            bs_str=bs_str.replace("ECP","ecp \""+element_str+"_"+bas_name+"\"")
            output_file.write(bs_str)
            #
            print(bas_name+" "+element_str)
    return

for bas_name in all_bs:
    md_bas_name = bas_name.lower()
    md_bas_name = md_bas_name.replace("*","_st_")
    md_bas_name = md_bas_name.replace("/","_sl_")
    writebs(md, bas_name, summary_file)
    if md[md_bas_name]['role'] == 'orbital':
        writebs(md, bas_name, summary_file, get_aux=1)
print("end")
