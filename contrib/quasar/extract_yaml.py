
emitter_yaml = False
emitter_ruamel = False

import sys
try:
    try:
        import ruamel.yaml as ruamel
    except ImportError:
        import ruamel_yaml as ruamel
    emitter_ruamel = True
except ImportError:
    import yaml
    emitter_yaml = True
    

preamble="""
"$schema": https://microsoft.com/qchem-0.1.schema.json
"""

def extract_fields():
    data = {} #yaml.load(initial_input)
    data['format'] = {'version' : '0.1'}
    data['bibliography'] = [{'url' : 'https://www.nwchem-sw.org'}]
    data['generator'] = {'source' : 'nwchem',
                         'version' : '6.8'}
    coulomb_repulsion = None
    scf_energy = None
    scf_energy_offset = None
    energy_offset = None
    one_electron_integrals = None
    two_electron_integrals = None
    n_electrons_alpha = None
    n_electrons_beta = None
    n_orbitals = None
    reader_mode = ""
    for line in sys.stdin.readlines():
        ln = line.strip()
        ln_segments = ln.split()
        if len(ln) == 0 or ln[0]=="#": #blank or comment line
            continue
        if reader_mode=="":
            if ln == "============================== echo of input deck ==============================":
                reader_mode = "input_deck"
            elif ln_segments[:2] == ["enrep_tce", "="]:
                coulomb_repulsion = {
                    'units' : 'hartree',
                    'value' : float(ln_segments[2])
                }
            elif ln_segments[:2] == ["EHF(total)", "="]:
                scf_energy = {
                    'units' : 'hartree',
                    'value' : float(ln_segments[2])
                }
            elif ln_segments[:3] == ["Shift", "(HFtot-HFA)", "="]:
                scf_energy_offset = {
                    'units' : 'hartree',
                    'value' : float(ln_segments[3])
                }
                energy_offset = {
                    'units' : 'hartree',
                    'value' : float(ln_segments[3])
                }
            elif ln == "begin_one_electron_integrals":
                reader_mode = "one_electron_integrals"
                one_electron_integrals = {
                    'units' : 'hartree',
                    'format' : 'sparse',
                    'values' : []
                }
            elif ln == "begin_two_electron_integrals":
                reader_mode = "two_electron_integrals"
                two_electron_integrals = {
                    'units' : 'hartree',
                    'format' : 'sparse',
                    'index_convention' : 'mulliken',
                    'values' : []
                }
            elif ln_segments[0] == "q_ele_a":
                assert len(ln_segments) == 2
                n_electrons_alpha = int(ln_segments[1])
            elif ln_segments[0] == "q_ele_b":
                assert len(ln_segments) == 2
                n_electrons_beta = int(ln_segments[1])
            elif ln_segments[0] == "q_orb":
                assert len(ln_segments) == 2
                n_orbitals = int(ln_segments[1])
        elif reader_mode == "input_deck":
            if ln == "================================================================================":
                reader_mode = ""
            elif ln_segments[0:2]== ["geometry", "units"]:
                reader_mode = "input_geometry"
                assert len(ln_segments) == 3
                geometry = {'coordinate_system': 'cartesian'}
                geometry['units'] = ln_segments[2]
            elif ln == "basis":
                reader_mode = "input_basis"
        elif reader_mode=="input_geometry":
            if ln=="end":
                reader_mode = "input_deck"
            elif ln_segments[0] == "symmetry":
                assert len(ln_segments) == 2
                geometry['symmetry'] = ln_segments[1]
            else: #atom description line
                assert len(ln_segments) == 4
                if not 'atoms' in geometry:
                    geometry['atoms'] = []
                geometry['atoms'] += [{"name":ln_segments[0],
                                       "coords":
                                       [float(ln_segments[1]), float(ln_segments[2]), float(ln_segments[3])]}]
        elif reader_mode == "input_basis":
            if ln == "end":
                reader_mode = "input_deck"
            else:
                assert len(ln_segments) == 3
                assert ln_segments[1] == "library"
                basis_set = {}
                basis_set['name'] = ln_segments[2]
                basis_set['type'] = 'gaussian'
        elif reader_mode == "one_electron_integrals":
            if ln == "end_one_electron_integrals":
                reader_mode = ""
            else:
                assert len(ln_segments) == 3
                one_electron_integrals['values'] += [[
                    int(ln_segments[0]),
                    int(ln_segments[1]),
                    float(ln_segments[2])
                ]]
        elif reader_mode == "two_electron_integrals":
            if ln == "end_two_electron_integrals":
                reader_mode = ""
            else:
                assert len(ln_segments) == 5
                two_electron_integrals['values'] += [[
                     int(ln_segments[0]),
                      int(ln_segments[1]),
                      int(ln_segments[2]),
                      int(ln_segments[3]),
                      float(ln_segments[4])
                    ]]
    fci_energy = {"units" : "hartree", "value" : 0.0, "upper":0.0, "lower":0.0}

    assert one_electron_integrals is not None, "one_electron_integrals is missing from NWChem output. Required to extract YAML"
    assert two_electron_integrals is not None, "two_electron_integrals is missing from NWChem output. Required to extract YAML"
    assert geometry is not None, "geometry information is missing from NWChem output. Required to extract YAML"
    assert basis_set  is not None, "basis_set is missing from NWChem output. Required to extract YAML"
    assert coulomb_repulsion  is not None, "coulomb_repulsion is missing from NWChem output. Required to extract YAML"
    assert scf_energy  is not None, "scf_energy is missing from NWChem output. Required to extract YAML"
    assert scf_energy_offset  is not None, "scf_energy_offset is missing from NWChem output. Required to extract YAML"
    assert energy_offset  is not None, "energy_offset is missing from NWChem output. Required to extract YAML"
    assert fci_energy  is not None, "fci_energy is missing from NWChem output. Required to extract YAML"
    assert n_orbitals  is not None, "n_orbitals is missing from NWChem output. Required to extract YAML"
    assert n_electrons_alpha  is not None, "n_electrons_alpha is missing from NWChem output. Required to extract YAML"
    assert n_electrons_beta  is not None, "n_electrons_beta is missing from NWChem output. Required to extract YAML"
    hamiltonian = {'one_electron_integrals' : one_electron_integrals,
                   'two_electron_integrals' : two_electron_integrals}
    integral_sets =  [{"metadata": { 'molecule_name' : 'unknown'},
                       "geometry":geometry,
                       "basis_set":basis_set,
                       "coulomb_repulsion" : coulomb_repulsion,
                       "scf_energy" : scf_energy,
                       "scf_energy_offset" : scf_energy_offset,
                       "energy_offset" : energy_offset,
                       "fci_energy" : fci_energy,
                       "hamiltonian" : hamiltonian,
                       "n_orbitals" : n_orbitals,
                       "n_electrons" : n_electrons_alpha + n_electrons_beta }]
    data['integral_sets'] = integral_sets
    return data
    
def emitter_ruamel_func():
    yaml = ruamel.YAML(typ="safe")
    #yaml.default_flow_style = False
    print(preamble)
    data = extract_fields()
    yaml.dump(data, sys.stdout)

def emitter_yaml_func():
    print(preamble)
    data = extract_fields()
    print(yaml.dump(data, default_flow_style=False))

    
# def emitter_yaml_func():
#     print(preamble)
#     yformat = {"version" : "0.1"}
#     generator = {"source" : "nwchem",
#                  "version" : "6.8"}
#     bibliography = [ {"url" : "https://www.nwchem-sw.org"}] # "doi": "10.1016/j.cpc.2010.04.018"}
#     geometry = { "coordinate_system": "cartesian"}
#     basis_set = {"type" : "gaussian"}
#     metadata = {"molecule" : "unknown"}
#     coulomb_repulsion = {"units" : "hartree"}
#     scf_energy = {"units" : "hartree"}
#     scf_energy_offset = {"units" : "hartree"}
#     energy_offset = {"units" : "hartree"}
#     fci_energy = {"units" : "hartree", "value" : 0.0, "upper":0.0, "lower":0.0}
#     hamiltonian = {}
#     one_electron_integrals = {"units" : "hartree",
#                               "format" : "sparse",
#                               "values" : []}
#     two_electron_integrals = {"units" : "hartree",
#                               "format" : "sparse",
#                               "index_convention": "mulliken",
#                               "values" : []}
#     n_orbitals = 0
#     n_electrons = None

#     reader_mode = ""
#     for line in sys.stdin.readlines():
#         ln = line.strip()
#         ln_segments = ln.split()
#         if len(ln) == 0 or ln[0]=="#": #blank or comment line
#             continue
#         if reader_mode=="":
#             if ln == "============================== echo of input deck ==============================":
#                 reader_mode = "input_deck"
#             elif ln_segments[:2] == ["enrep_tce", "="]:
#                 coulomb_repulsion["value"] = float(ln_segments[2])
#             elif ln_segments[:2] == ["EHF(total)", "="]:
#                 scf_energy["value"] = float(ln_segments[2])
#             elif ln_segments[:3] == ["Shift", "(HFtot-HFA)", "="]:
#                 scf_energy_offset["value"] = float(ln_segments[3])
#                 energy_offset["value"] = float(ln_segments[3])
#             elif ln == "begin_one_electron_integrals":
#                 reader_mode = "one_electron_integrals"
#             elif ln == "begin_two_electron_integrals":
#                 reader_mode = "two_electron_integrals"
#             elif ln_segments[0] == "q_ele_a":
#                 assert len(ln_segments) == 2
#                 if n_electrons is None:
#                     n_electrons = int(ln_segments[1])
#                 else:
#                     n_electrons += int(ln_segments[1])
#             elif ln_segments[0] == "q_ele_b":
#                 assert len(ln_segments) == 2
#                 if n_electrons is None:
#                     n_electrons = int(ln_segments[1])
#                 else:
#                     n_electrons += int(ln_segments[1])
#             elif ln_segments[0] == "q_orb":
#                 assert len(ln_segments) == 2
#                 n_orbitals = int(ln_segments[1])
#         elif reader_mode == "input_deck":
#             if ln == "================================================================================":
#                 reader_mode = ""
#             elif ln_segments[0:2]== ["geometry", "units"]:
#                 reader_mode = "input_geometry"
#                 assert len(ln_segments) == 3
#                 geometry["units"] = ln_segments[2]
#                 #geometry["symmetry"] = ln.split()[1]
#             elif ln == "basis":
#                 reader_mode = "input_basis"
#         elif reader_mode=="input_geometry":
#             if ln=="end":
#                 reader_mode = "input_deck"
#             elif ln_segments[0] == "symmetry":
#                 assert len(ln_segments) == 2
#                 geometry["symmetry"] = ln_segments[1]
#             else: #atom description line
#                 assert len(ln_segments) == 4
#                 if not "atoms" in geometry:
#                     geometry["atoms"] = []
#                 geometry["atoms"] += [{"name":ln_segments[0], "coords":
#                                        [float(ln_segments[1]), float(ln_segments[2]), float(ln_segments[3])]}]
#         elif reader_mode == "input_basis":
#             if ln == "end":
#                 reader_mode = "input_deck"
#             else:
#                 assert len(ln_segments) == 3
#                 assert ln_segments[1] == "library"
#                 basis_set["name"] = ln_segments[2]
#         elif reader_mode == "one_electron_integrals":
#             if ln == "end_one_electron_integrals":
#                 reader_mode = ""
#             else:
#                 assert len(ln_segments) == 3
#                 one_electron_integrals["values"] += [[
#                     int(ln_segments[0]),
#                     int(ln_segments[1]),
#                     float(ln_segments[2])
#                 ]]
#         elif reader_mode == "two_electron_integrals":
#             if ln == "end_two_electron_integrals":
#                 reader_mode = ""
#             else:
#                 assert len(ln_segments) == 5
#                 two_electron_integrals["values"] += [[
#                      int(ln_segments[0]),
#                       int(ln_segments[1]),
#                       int(ln_segments[2]),
#                       int(ln_segments[3]),
#                       float(ln_segments[4])
#                     ]]
#     hamiltonian["one_electron_integrals"] = one_electron_integrals
#     hamiltonian["two_electron_integrals"] = two_electron_integrals 
#     output = {"format" : yformat,
#               "generator" : generator,
#               "bibliography" : bibliography,
#               "integral_sets":
#               [{"metadata":metadata,
#                 "geometry":geometry,
#                 "basis_set":basis_set,
#                 "coulomb_repulsion" : coulomb_repulsion,
#                 "scf_energy" : scf_energy,
#                 "scf_energy_offset" : scf_energy_offset,
#                 "energy_offset" : energy_offset,
#                 "fci_energy" : fci_energy,
#                 "hamiltonian" : hamiltonian,
#                 "n_orbitals" : n_orbitals,
#                 "n_electrons" : n_electrons }]}
#     print(yaml.dump(output, default_flow_style=False))

def main():
    assert emitter_yaml or emitter_ruamel, "Extraction failed: could not import YAML or RUAMEL packages."
    if emitter_yaml:
        emitter_yaml_func()
    elif emitter_ruamel:
        emitter_ruamel_func()
    else:
        assert False #unreachable code

if __name__ == "__main__":
    main()
