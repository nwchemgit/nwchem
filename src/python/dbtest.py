from nwchem import *

# nwchem_init has to be called to set up tcgmsg
db = nwchem_init("h2o")
db.ls()
