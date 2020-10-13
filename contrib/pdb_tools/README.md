Tools to manipulate PDB files
=============================

General PDB files need tools to apply transformations to them. As this need is
quite common there a number of tools out there already and where applicable
you should use them. Examples of such tools are:

* convpdb.pl
* pdb4amber
* charmmlipid2amber.py
* pdbfixer

Unfortunately some of these tools try to be too helpful. For example `convpdb.pl`
will not only scale the coordinates in a PDB file if you ask for it, but it will
also drop all TER records, rename HOH to TIP3, drop the CRYST1 record, and drop
the chemical symbol at the end of the line. This might leave you with a useless
result. 

The tools in this directory do exactly what it says on the tin, nothing less
and nothing more. The tools provided here are:

* pdb_amber2nwchem 
  * Convert an Amber PDB file into an NWChem PDB file by renaming atoms and
    residues.
* pdb_large
  * Interconvert a PDB file between the regular PDB file format and the
    NWChem PDB file format for large molecular systems. In essence the 
    large PDB format allows 6 digits for the residue number instead of 4.
* pdb_scale.F
  * Scale the atomic coordinates and the box size with a specified factor.
* pdb_supercell.py
  * Build a supercell of a given PDB file by doubling the initial molecular
    system in the requested direction.
