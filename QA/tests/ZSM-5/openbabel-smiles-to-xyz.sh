#!/bin/bash
obabel -h --gen3d -ismi SMILES/h2o.smi        -oxyz -Oxyz/h2o.xyz
obabel -h --gen3d -ismi SMILES/1-propanol.smi -oxyz -Oxyz/1-propanol.xyz
obabel -h --gen3d -ismi SMILES/2-propanol.smi -oxyz -Oxyz/2-propanol.xyz
obabel -h --gen3d -ismi SMILES/propene.smi    -oxyz -Oxyz/propene.xyz
