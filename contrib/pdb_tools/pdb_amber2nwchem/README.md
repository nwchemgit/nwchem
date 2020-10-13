# README PDB_AMBER2NWCHEM

`pdb_amber2nwchem` is a tool to convert PDB files from the Amber
dialect into the NWChem dialect. PDB files have standardized atom
names for many residues of biological relevance. However, most PDB
files do not contain Hydrogen atoms and different conventions have
been adopted by different code communities. This tool converts
a PDB file in the Amber convention to one following the NWChem
conventions, so that NWChem can read it successfully.

Apart from atom names the tool also converts some residue names:
- "WAT" --> "HOH"
- "Na+" --> "Na"
- "K+"  --> "K"
- "Cl-" --> "Cl"

to name a few examples.

# References

- [1] "IUPAC-IUB Commission on Biochemical Nomenclature. 
      Abbreviations and Symbols for the Description of the
      Conformation of Polypeptide Chains. Tentative Rules (1969)".
      Biochemistry (1970) 9, 3471-3479.
      DOI: [10.1021/bi00820a001](https://doi.org/10.1021/bi00820a001)

- [2] Charles Hoogstraten, "Correlation of hydrogen atom naming
      systems, including diastereotopic protons." Web-page:
      http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl 
      [accessed Jan 6, 2018]

- [3] IUPAC-IUBMB Joint Commission on Biochemical Nomenclature
      and Nomenclature Commission of IUBMB, "Biochemical
      Nomenclature and Related Documents (White Book)." 2nd Edition,
      Portland Press, 1992, pp. 39-69. ISBN: 1-85578-005-4.
      Part: "Amino acids, peptides and proteins", Section 3AA-1
      "Names of Common alpha-Amino Acids". 
      DOI: [10.1351/pac198456050595](https://doi.org/10.1351/pac198456050595)
