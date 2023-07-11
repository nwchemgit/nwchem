#!/bin/bash

for atom in \
          'H' 'He' 'Li' 'Be' 'B' 'C' 'N' 'O' 'F' 'Ne'       \
          'Na' 'Mg' 'Al' 'Si' 'P' 'S' 'Cl' 'Ar' 'K' 'Ca'       \
          'Sc' 'Ti' 'V' 'Cr' 'Mn' 'Fe' 'Co' 'Ni' 'Cu' 'Zn'       \
          'Ga' 'Ge' 'As' 'Se' 'Br' 'Kr' 'Rb' 'Sr' 'Y' 'Zr'       \
          'Nb' 'Mo' 'Tc' 'Ru' 'Rh' 'Pd' 'Ag' 'Cd' 'In' 'Sn'       \
          'Sb' 'Te' 'I' 'Xe' 'Cs' 'Ba' 'La' 'Ce' 'Pr' 'Nd'       \
          'Pm' 'Sm' 'Eu' 'Gd' 'Tb' 'Dy' 'Ho' 'Er' 'Tm' 'Yb'       \
          'Lu' 'Hf' 'Ta' 'W' 'Re' 'Os' 'Ir' 'Pt' 'Au' 'Hg'       \
          'Tl' 'Pb' 'Bi' 'Po' 'At' 'Rn' 'Fr' 'Ra' 'Ac' 'Th'       \
          'Pa' 'U' 'Np' 'Pu' 'Am' 'Cm' 'Bk' 'Cf' 'Es' 'Fm'       \
          'Md' 'No' 'Lr' 'Rf' 'Db' 'Sg' 'Bh' 'Hs' 'Mt' 'Ds'       \
          'Rg' 'Cn' 'Nh' 'Fl' 'Mc' 'Lv' 'Ts' 'Og' 'Uue' 'Ubn'       \
          'Hydrogen' 'Helium' 'Lithium' 'Beryllium' 'Boron'            \
          'Carbon' 'Nitrogen' 'Oxygen' 'Fluorine' 'Neon' 'Sodium'     \
          'Magnesium' 'Aluminium' 'Silicon' 'Phosphorous'               \
          'Sulphur' 'Chlorine' 'Argon' 'Potassium' 'Calcium'           \
          'Scandium' 'Titanium' 'Vanadium' 'Chromium' 'Manganese'      \
          'Iron' 'Cobalt' 'Nickel' 'Copper' 'Zinc' 'Gallium'          \
          'Germanium' 'Arsenic' 'Selenium' 'Bromine' 'Krypton'         \
          'Rubidium' 'Strontium' 'Yttrium' 'Zirconium' 'Niobium'       \
          'Molybdenum' 'Technetium' 'Ruthenium' 'Rhodium'               \
          'Palladium' 'Silver' 'Cadmium' 'Indium' 'Tin'                \
          'Antinomy' 'Tellurium' 'Iodine' 'Xenon' 'Caesium'            \
          'Barium' 'Lanthanum' 'Cerium' 'Praseodymium' 'Neodymium'     \
          'Promethium' 'Samarium' 'Europium' 'Gadolinium'               \
          'Terbium' 'Dysprosium' 'Holmium' 'Erbium' 'Thulium'          \
          'Ytterbium' 'Lutetium' 'Hafnium' 'Tantalum' 'Tungsten'       \
          'Rhenium' 'Osmium' 'Iridium' 'Platinum' 'Gold'               \
          'Mercury' 'Thallium' 'Lead' 'Bismuth' 'Polonium'             \
          'Astatine' 'Radon' 'Francium' 'Radium' 'Actinium'            \
          'Thorium' 'Protoactinium' 'Uranium' 'Neptunium'               \
          'Plutonium' 'Americium' 'Curium' 'Berkelium'                  \
          'Californium' 'Einsteinium' 'Fermium' 'Mendelevium'           \
          'Nobelium' 'Lawrencium' 'Rutherfordium' 'Dubnium'             \
          'Seaborgium' 'Bohrium' 'Hassium' 'Meitnerium'                 \
          'Darmstadtium' 'Roentgenium' 'Copernicium' 'Nihonium'         \
          'Flerovium' 'Moscovium' 'Livermorium' 'Tennessine'            \
          'Oganesson' 'Ununennium' 'Unbinilium'
do

echo "
start scratch_${atom}

basis spherical bse
 * library AHGBS-5
end

geometry units bohr
  ${atom} 0.0 0.0 0.0
end

relativistic
  douglas-kroll dkh
end

dft
  odft
  xc slater
end

#print rtdb
#print rtdbvalues

task dft ignore
" > ${atom}.nw || echo "${atom}" failed

done
