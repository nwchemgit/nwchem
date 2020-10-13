# README PDB_LARGE

NWChem has a special file format for large PDB files. When PDB files are
the residue numbers can no longer be listed with just 4 digits. In the
large format 6 digits can be used. However, for applications it is 
useful to be able to convert to and from the large format.

The large format is specified by providing the PDB keyword "LRGPDB"
before the first "ATOM" or "HETATM" keywork.

This converter detects what format the input PDB file is in and 
generates a PDB output file in the other format.
