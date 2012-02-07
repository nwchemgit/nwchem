'''
testing pydoc
'''
import re
import string
import sys

class AtomDictionary():
    '''
    Stores static information about atoms
    such as their full name (e.g. Hydrogen)
    Van der Waals radius, covalent radius, and weights.
    Private dictionaries are used.
    Also includes helper function to get element name from atom name. e.g. to infer that
    3HW or 3hW is actually is H.
    '''
    @staticmethod 
    def genProperty(a,d):
        '''
        Helper function to assist accessing atom
        dictionaries by digesting "nonstandard" atom name
        such O2 or 3HW through succession of steps
        1. converting into all caps (h -> H)
        2. removing all digits (1HW -> HW)
        3. stripping to just the first found letter (1HW -> H)
        '''
        
        a0 = string.upper(a)
        try:
            return d[a0.strip()],a0
        except KeyError:
            pass
                
        s1=re.match(r" *[0-9]*([A-Z]+) *", a0)
        try:
            return d[s1.group(1)],s1.group(1)
        except KeyError:
            pass
        
        s1=re.match(r" *[0-9]*([A-Z] *)", a0)
        try:
            return d[s1.group(1)],s1.group(1)
        except KeyError:
            print "Cannot match name", a
            sys.exit(1)       
    
    @staticmethod 
    def fullName(s):
        return AtomDictionary.genProperty(s, AtomDictionary.__fulename)[0]
    
    @staticmethod 
    def elementName(s):
        return AtomDictionary.genProperty(s, AtomDictionary.__fulename)[1]
    
    @staticmethod 
    def vdwRadius(s):
        return AtomDictionary.genProperty(s, AtomDictionary.__vdw)[0]
    
    @staticmethod 
    def covRadius(s):
        return AtomDictionary.genProperty(s, AtomDictionary.__rc)[0]
    
    @staticmethod 
    def weight(s):
        return AtomDictionary.genProperty(s, AtomDictionary.__weight)[0]


    
    __weight = {}
    __weight['H']=1.008000
    __weight['He']=4.003000
    __weight['Li']=6.941000
    __weight['Be']=9.012000
    __weight['B']=10.811000
    __weight['C']=12.011000
    __weight['N']=14.007000
    __weight['O']=15.999000
    __weight['F']=18.998000
    __weight['Ne']=20.180000
    __weight['Na']=22.991000
    __weight['Mg']=24.305000
    __weight['Al']=26.982000
    __weight['Si']=28.086000
    __weight['P']=30.974000
    __weight['S']=32.066000
    __weight['Cl']=35.453000
    __weight['Ar']=39.948000
    __weight['K']=39.098000
    __weight['Ca']=40.078000
    __weight['Sc']=44.956000
    __weight['Ti']=47.867000
    __weight['V']=50.942000
    __weight['Cr']=51.996000
    __weight['Mn']=54.938000
    __weight['Fe']=55.845000
    __weight['Co']=58.933000
    __weight['Ni']=58.693000
    __weight['Cu']=63.546000
    __weight['Zn']=65.390000
    __weight['Ga']=69.723000
    __weight['Ge']=72.610000
    __weight['As']=74.922000
    __weight['Se']=78.960000
    __weight['Br']=79.904000
    __weight['Kr']=83.800000
    __weight['Rb']=85.468000
    __weight['Sr']=87.620000
    __weight['Y']=88.906000
    __weight['Zr']=91.224000
    __weight['Nb']=92.906000
    __weight['Mo']=95.940000
    __weight['Tc']=0.000000
    __weight['Ru']=101.070000
    __weight['Rh']=102.906000
    __weight['Pd']=106.420000
    __weight['Ag']=107.868000
    __weight['Cd']=112.411000
    __weight['In']=114.818000
    __weight['Sn']=118.710000
    __weight['Sb']=121.760000
    __weight['Te']=127.600000
    __weight['I']=126.904000
    __weight['Xe']=131.290000
    __weight['Cs']=132.905000
    __weight['Ba']=137.327000
    __weight['Lu']=174.967000
    __weight['Hf']=178.490000
    __weight['Ta']=180.948000
    __weight['W']=183.840000
    __weight['Re']=186.207000
    __weight['Os']=190.230000
    __weight['Ir']=192.217000
    __weight['Pt']=195.078000
    __weight['Au']=196.967000
    __weight['Hg']=200.590000
    __weight['Tl']=204.383000
    __weight['Pb']=207.200000
    __weight['Bi']=208.980000
    __weight['Po']=0.000000
    __weight['At']=0.000000
    __weight['Rn']=0.000000
    __weight['Ce']=140.116000
    __weight['Dy']=162.500000
    __weight['Er']=167.260000
    __weight['Eu']=151.964000
    __weight['Gd']=157.250000
    __weight['Ho']=164.930000
    __weight['La']=138.906000
    __weight['Nd']=144.240000
    __weight['Pm']=0.000000
    __weight['Pr']=140.908000
    __weight['Sm']=150.360000
    __weight['Tb']=158.925000
    __weight['Tm']=168.934000
    __weight['Yb']=173.040000
    __weight['Fr']=0.000000
    __weight['Ra']=0.000000
    __weight['Lr']=0.000000
    __weight['Rf']=0.000000
    __weight['Db']=0.000000
    __weight['Sg']=0.000000
    __weight['Bh']=0.000000
    __weight['Hs']=0.000000
    __weight['Mt']=0.000000
    __weight['Ds']=0.000000
    __weight['Ac']=0.000000
    __weight['Am']=0.000000
    __weight['Bk']=0.000000
    __weight['Cf']=0.000000
    __weight['Cm']=0.000000
    __weight['Es']=0.000000
    __weight['Fm']=0.000000
    __weight['Md']=0.000000
    __weight['No']=0.000000
    __weight['Np']=0.000000
    __weight['Pa']=231.036000
    __weight['Pu']=0.000000
    __weight['Th']=232.038000
    __weight['U']=238.029000

    __rc = {}
    __rc['H']=0.230000
    __rc['He']=1.500000
    __rc['Li']=1.280000
    __rc['Be']=0.960000
    __rc['B']=0.830000
    __rc['C']=0.680000
    __rc['N']=0.680000
    __rc['O']=0.680000
    __rc['F']=0.640000
    __rc['Ne']=1.500000
    __rc['Na']=1.660000
    __rc['Mg']=1.410000
    __rc['Al']=1.210000
    __rc['Si']=1.200000
    __rc['P']=1.050000
    __rc['S']=1.020000
    __rc['Cl']=0.990000
    __rc['Ar']=1.510000
    __rc['K']=2.030000
    __rc['Ca']=1.760000
    __rc['Sc']=1.700000
    __rc['Ti']=1.600000
    __rc['V']=1.530000
    __rc['Cr']=1.390000
    __rc['Mn']=1.610000
    __rc['Fe']=1.520000
    __rc['Co']=1.260000
    __rc['Ni']=1.240000
    __rc['Cu']=1.320000
    __rc['Zn']=1.220000
    __rc['Ga']=1.220000
    __rc['Ge']=1.170000
    __rc['As']=1.210000
    __rc['Se']=1.220000
    __rc['Br']=1.210000
    __rc['Kr']=1.500000
    __rc['Rb']=2.200000
    __rc['Sr']=1.950000
    __rc['Y']=1.900000
    __rc['Zr']=1.750000
    __rc['Nb']=1.640000
    __rc['Mo']=1.540000
    __rc['Tc']=1.470000
    __rc['Ru']=1.460000
    __rc['Rh']=1.450000
    __rc['Pd']=1.390000
    __rc['Ag']=1.450000
    __rc['Cd']=1.440000
    __rc['In']=1.420000
    __rc['Sn']=1.390000
    __rc['Sb']=1.390000
    __rc['Te']=1.470000
    __rc['I']=1.400000
    __rc['Xe']=1.500000
    __rc['Cs']=2.440000
    __rc['Ba']=2.150000
    __rc['Lu']=1.870000
    __rc['Hf']=1.750000
    __rc['Ta']=1.700000
    __rc['W']=1.620000
    __rc['Re']=1.510000
    __rc['Os']=1.440000
    __rc['Ir']=1.410000
    __rc['Pt']=1.360000
    __rc['Au']=1.500000
    __rc['Hg']=1.320000
    __rc['Tl']=1.450000
    __rc['Pb']=1.460000
    __rc['Bi']=1.480000
    __rc['Po']=1.400000
    __rc['At']=1.210000
    __rc['Rn']=1.500000
    __rc['Ce']=2.040000
    __rc['Dy']=1.920000
    __rc['Er']=1.890000
    __rc['Eu']=1.980000
    __rc['Gd']=1.960000
    __rc['Ho']=1.920000
    __rc['La']=2.070000
    __rc['Nd']=2.010000
    __rc['Pm']=1.990000
    __rc['Pr']=2.030000
    __rc['Sm']=1.980000
    __rc['Tb']=1.940000
    __rc['Tm']=1.900000
    __rc['Yb']=1.870000
    __rc['Fr']=2.600000
    __rc['Ra']=2.210000
    __rc['Lr']=1.500000
    __rc['Rf']=1.500000
    __rc['Db']=1.500000
    __rc['Sg']=1.500000
    __rc['Bh']=1.500000
    __rc['Hs']=1.500000
    __rc['Mt']=1.500000
    __rc['Ds']=1.500000
    __rc['Ac']=2.150000
    __rc['Am']=1.800000
    __rc['Bk']=1.540000
    __rc['Cf']=1.830000
    __rc['Cm']=1.690000
    __rc['Es']=1.500000
    __rc['Fm']=1.500000
    __rc['Md']=1.500000
    __rc['No']=1.500000
    __rc['Np']=1.900000
    __rc['Pa']=2.000000
    __rc['Pu']=1.870000
    __rc['Th']=2.060000
    __rc['U']=1.960000
    
    __fulename = {}
    __fulename['H']='Hydrogen'
    __fulename['He']='Helium'
    __fulename['Li']='Lithium'
    __fulename['Be']='Beryllium'
    __fulename['B']='Boron'
    __fulename['C']='Carbon'
    __fulename['N']='Nitrogen'
    __fulename['O']='Oxygen'
    __fulename['F']='Fluorine'
    __fulename['Ne']='Neon'
    __fulename['Na']='Sodium'
    __fulename['Mg']='Magnesium'
    __fulename['Al']='Aluminium'
    __fulename['Si']='Silicon'
    __fulename['P']='Phosphorus'
    __fulename['S']='Sulphur'
    __fulename['Cl']='Chlorine'
    __fulename['Ar']='Argon'
    __fulename['K']='Potassium'
    __fulename['Ca']='Calcium'
    __fulename['Sc']='Scandium'
    __fulename['Ti']='Titanium'
    __fulename['V']='Vanadium'
    __fulename['Cr']='Chromium'
    __fulename['Mn']='Manganese'
    __fulename['Fe']='Iron'
    __fulename['Co']='Cobalt'
    __fulename['Ni']='Nickel'
    __fulename['Cu']='Copper'
    __fulename['Zn']='Zinc'
    __fulename['Ga']='Gallium'
    __fulename['Ge']='Germanium'
    __fulename['As']='Arsenic'
    __fulename['Se']='Selenium'
    __fulename['Br']='Bromine'
    __fulename['Kr']='Krypton'
    __fulename['Rb']='Rubidium'
    __fulename['Sr']='Strontium'
    __fulename['Y']='Yttrium'
    __fulename['Zr']='Zirconium'
    __fulename['Nb']='Niobium'
    __fulename['Mo']='Molybdenum'
    __fulename['Tc']='Technetium'
    __fulename['Ru']='Ruthenium'
    __fulename['Rh']='Rhodium'
    __fulename['Pd']='Palladium'
    __fulename['Ag']='Silver'
    __fulename['Cd']='Cadmium'
    __fulename['In']='Indium'
    __fulename['Sn']='Tin'
    __fulename['Sb']='Antimony'
    __fulename['Te']='Tellurium'
    __fulename['I']='Iodine'
    __fulename['Xe']='Xenon'
    __fulename['Cs']='Caesium'
    __fulename['Ba']='Barium'
    __fulename['Lu']='Lutetium'
    __fulename['Hf']='Hafnium'
    __fulename['Ta']='Tantalum'
    __fulename['W']='Tungsten'
    __fulename['Re']='Rhenium'
    __fulename['Os']='Osmium'
    __fulename['Ir']='Iridium'
    __fulename['Pt']='Platinum'
    __fulename['Au']='Gold'
    __fulename['Hg']='Mercury'
    __fulename['Tl']='Thallium'
    __fulename['Pb']='Lead'
    __fulename['Bi']='Bismuth'
    __fulename['Po']='Polonium'
    __fulename['At']='Astatine'
    __fulename['Rn']='Radon'
    __fulename['Ce']='Cerium'
    __fulename['Dy']='Dysprosium'
    __fulename['Er']='Erbium'
    __fulename['Eu']='Europium'
    __fulename['Gd']='Gadolinium'
    __fulename['Ho']='Holmium'
    __fulename['La']='Lanthanum'
    __fulename['Nd']='Neodymium'
    __fulename['Pm']='Promethium'
    __fulename['Pr']='Praseodymium'
    __fulename['Sm']='Samarium'
    __fulename['Tb']='Terbium'
    __fulename['Tm']='Thulium'
    __fulename['Yb']='Ytterbium'
    __fulename['Fr']='Francium'
    __fulename['Ra']='Radium'
    __fulename['Lr']='Lawrencium'
    __fulename['Rf']='Rutherfordium'
    __fulename['Db']='Dubnium'
    __fulename['Sg']='Seaborgium'
    __fulename['Bh']='Bohrium'
    __fulename['Hs']='Hassium'
    __fulename['Mt']='Meitnerium'
    __fulename['Ds']='Darmstadtium'
    __fulename['Ac']='Actinium'
    __fulename['Am']='Americium'
    __fulename['Bk']='Berkelium'
    __fulename['Cf']='Californium'
    __fulename['Cm']='Curium'
    __fulename['Es']='Einsteinium'
    __fulename['Fm']='Fermium'
    __fulename['Md']='Mendelevium'
    __fulename['No']='Nobelium'
    __fulename['Np']='Neptunium'
    __fulename['Pa']='Protactinium'
    __fulename['Pu']='Plutonium'
    __fulename['Th']='Thorium'
    __fulename['U']='Uranium'
    
    __vdw={}
    __vdw['H']=1.090000
    __vdw['He']=1.400000
    __vdw['Li']=1.820000
    __vdw['Be']=2.000000
    __vdw['B']=2.000000
    __vdw['C']=1.700000
    __vdw['N']=1.550000
    __vdw['O']=1.520000
    __vdw['F']=1.470000
    __vdw['Ne']=1.540000
    __vdw['Na']=2.270000
    __vdw['Mg']=1.730000
    __vdw['Al']=2.000000
    __vdw['Si']=2.100000
    __vdw['P']=1.800000
    __vdw['S']=1.800000
    __vdw['Cl']=1.750000
    __vdw['Ar']=1.880000
    __vdw['K']=2.750000
    __vdw['Ca']=2.000000
    __vdw['Sc']=2.000000
    __vdw['Ti']=2.000000
    __vdw['V']=2.000000
    __vdw['Cr']=2.000000
    __vdw['Mn']=2.000000
    __vdw['Fe']=2.000000
    __vdw['Co']=2.000000
    __vdw['Ni']=1.630000
    __vdw['Cu']=1.400000
    __vdw['Zn']=1.390000
    __vdw['Ga']=1.870000
    __vdw['Ge']=2.000000
    __vdw['As']=1.850000
    __vdw['Se']=1.900000
    __vdw['Br']=1.850000
    __vdw['Kr']=2.020000
    __vdw['Rb']=2.000000
    __vdw['Sr']=2.000000
    __vdw['Y']=2.000000
    __vdw['Zr']=2.000000
    __vdw['Nb']=2.000000
    __vdw['Mo']=2.000000
    __vdw['Tc']=2.000000
    __vdw['Ru']=2.000000
    __vdw['Rh']=2.000000
    __vdw['Pd']=1.630000
    __vdw['Ag']=1.720000
    __vdw['Cd']=1.580000
    __vdw['In']=1.930000
    __vdw['Sn']=2.170000
    __vdw['Sb']=2.000000
    __vdw['Te']=2.060000
    __vdw['I']=1.980000
    __vdw['Xe']=2.160000
    __vdw['Cs']=2.000000
    __vdw['Ba']=2.000000
    __vdw['Lu']=2.000000
    __vdw['Hf']=2.000000
    __vdw['Ta']=2.000000
    __vdw['W']=2.000000
    __vdw['Re']=2.000000
    __vdw['Os']=2.000000
    __vdw['Ir']=2.000000
    __vdw['Pt']=1.720000
    __vdw['Au']=1.660000
    __vdw['Hg']=1.550000
    __vdw['Tl']=1.960000
    __vdw['Pb']=2.020000
    __vdw['Bi']=2.000000
    __vdw['Po']=2.000000
    __vdw['At']=2.000000
    __vdw['Rn']=2.000000
    __vdw['Ce']=2.000000
    __vdw['Dy']=2.000000
    __vdw['Er']=2.000000
    __vdw['Eu']=2.000000
    __vdw['Gd']=2.000000
    __vdw['Ho']=2.000000
    __vdw['La']=2.000000
    __vdw['Nd']=2.000000
    __vdw['Pm']=2.000000
    __vdw['Pr']=2.000000
    __vdw['Sm']=2.000000
    __vdw['Tb']=2.000000
    __vdw['Tm']=2.000000
    __vdw['Yb']=2.000000
    __vdw['Fr']=2.000000
    __vdw['Ra']=2.000000
    __vdw['Lr']=2.000000
    __vdw['Rf']=2.000000
    __vdw['Db']=2.000000
    __vdw['Sg']=2.000000
    __vdw['Bh']=2.000000
    __vdw['Hs']=2.000000
    __vdw['Mt']=2.000000
    __vdw['Ds']=2.000000
    __vdw['Ac']=2.000000
    __vdw['Am']=2.000000
    __vdw['Bk']=2.000000
    __vdw['Cf']=2.000000
    __vdw['Cm']=2.000000
    __vdw['Es']=2.000000
    __vdw['Fm']=2.000000
    __vdw['Md']=2.000000
    __vdw['No']=2.000000
    __vdw['Np']=2.000000
    __vdw['Pa']=2.000000
    __vdw['Pu']=2.000000
    __vdw['Th']=2.000000
    __vdw['U']=1.860000


    
if __name__ == '__main__':
    print AtomDictionary.fullName(' H1  ')
    print AtomDictionary.vdwRadius('h')
    

