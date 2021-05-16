#!/usr/bin/env python
"""
nwchem_tddft_spectrum.py - Extract the tddft data from an NWChem output
                           file and plot it in the form of a absorption
                           spectrum.

Usage:
nwchem_tddft_spectrum.py (options) (files)

Options:
 -a  value  Use value as the Gaussian exponent to generate the plot
            (default=10000)
 -b         Print both spectra
 -h         Print this help
 -s         Print the singlet spectrum  (default)
 -t         Print the triplet spectrum
 -o         Print the open-shell spectrum
 -u  value  Use value as the units. Must be in [h,eV,kcal/mol]; default eV

Send any errors or additions to Rick Muller, rmuller@sandia.gov.
"""

defaults = {
    'alpha' : 10000,
    'plot_singlet': True,
    'plot_triplet': False,
    'plot_unrestricted': False,
    'plot_zeros': False, # Show zero oscillator strengths as slight negative blips
    'units': 'eV',
    }

from pylab import *

def parse_nwchem_output(fname):
    import re
    from pprint import pprint
    singlet_pat = re.compile('\s+Root\s+\d+\s+singlet')
    triplet_pat = re.compile('\s+Root\s+\d+\s+triplet')
    unrestricted_pat = re.compile('\s+Root\s+\d')
    singlets = []
    triplets = []
    unrestricted = []
    f = open(fname)
    while 1:
        line = f.readline()
        if not line: break
        if singlet_pat.search(line):
            energy,osc,trans_mom = parse_element(line,f)
            singlets.append((energy,osc))
        if triplet_pat.search(line):
            energy,osc,trans_mom = parse_element(line,f)
            triplets.append((energy,osc))
        if unrestricted_pat.search(line):
            energy,osc,trans_mom = parse_element(line,f)
            unrestricted.append((energy,osc))
    #pprint(singlets)
    #pprint(triplets)
    return singlets,triplets,unrestricted

def parse_element(line,f):
    import re
    oscillator_pat = re.compile('\s+Oscillator')
    words = line.split()
    try:
        energy = float(words[4])
    except:
        energy = float(words[3])
    line = f.readline() # skip ----
    line = f.readline()
    words = line.split()
    try:
        txyz = float(words[3]),float(words[5]),float(words[7])
    except:
        txyz = (0,0,0)
    line = f.readline()
    words = line.split()
    while None == oscillator_pat.search(line):
      line = f.readline()
      words = line.split()
    try:
        osc = float(words[3])
    except:
        osc = 0
    return energy,osc,txyz

def to_plot(data,**kwargs):
    alpha = kwargs.get('alpha',defaults['alpha'])
    plot_zeros = kwargs.get('plot_zeros',defaults['plot_zeros'])
    abmin_floor = kwargs.get('abmin',0.001)
    
    Emin = min(E for E,ab in data)
    Emax = max(E for E,ab in data)
    abmin = min([0] + [ab for E,ab in data if ab > 0])
    abmin = max(abmin,abmin_floor)
    Edel = Emax-Emin
    Es = linspace(Emin-Edel/10.,Emax+Edel/10.,2000)
    Spectrum = zeros(Es.shape,'d')
    for E,ab in data:
        if ab > 0:
            Spectrum += ab*exp(-alpha*(Es-E)**2)
        elif plot_zeros:
            Spectrum -= abmin*exp(-alpha*(Es-E)**2)
    return Es,Spectrum

def nwchem_tddft_spectrum(fname,**kwargs):
    from os.path import basename
    thetitle = kwargs.get('title',basename(fname))
    units = kwargs.get('units',defaults['units'])
    assert units in ['h','eV','ev','kcal/mol','kcm']
    if units in ['eV','ev']:
        conv = 27.211
    elif units in ['kcal/mol','kcm']:
        conv = 627.51
    else:
        conv = 1.0 # h
    singlet,triplet,unrestricted = parse_nwchem_output(fname)
    if thetitle: title(thetitle)
    if kwargs['plot_singlet']:
        Es0,Spectrum0 = to_plot(singlet,**kwargs)
        plot(Es0*conv,Spectrum0,label='sing')
    if kwargs['plot_triplet']:
        Es1,Spectrum1 = to_plot(triplet,**kwargs)
        plot(Es1*conv,Spectrum1,label='trip')
    if kwargs['plot_unrestricted']:
        Es1,Spectrum1 = to_plot(unrestricted,**kwargs)
        plot(Es1*conv,Spectrum1,label='open')
    legend()
    xlabel('Energy (%s)' % units)
    ylabel('Oscillator Strength')
    if fname.endswith('.out'):
        oname = fname.replace('.out','.png')
    else:
        oname = fname + ".png"
    savefig(oname)
    show()
    return

if __name__ == '__main__':
    from getopt import getopt
    opts,args = getopt(sys.argv[1:],'hstoa:b')
    kwargs = defaults
    for key,val in opts:
        if key == '-a':
            kwargs['alpha'] = float(val)
        if key == '-b':
            kwargs['plot_singlet'] = True
            kwargs['plot_triplet'] = True
            kwargs['plot_unrestricted'] = False
        if key == '-s':
            kwargs['plot_singlet'] = True
            kwargs['plot_triplet'] = False
            kwargs['plot_unrestricted'] = False
        if key == '-t':
            kwargs['plot_singlet'] = False
            kwargs['plot_triplet'] = True
            kwargs['plot_unrestricted'] = False
        if key == '-o':
            kwargs['plot_singlet'] = False
            kwargs['plot_triplet'] = False
            kwargs['plot_unrestricted'] = True
        if key == '-h':
            print(__doc__)
            sys.exit()
    if len(args) == 0:
        print(__doc__)
        sys.exit()
    for fname in args:
        nwchem_tddft_spectrum(fname,**kwargs)
    

        
