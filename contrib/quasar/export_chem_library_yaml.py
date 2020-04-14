
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
"$schema": https://raw.githubusercontent.com/Microsoft/Quantum/master/Chemistry/Schema/broombridge-0.1.schema.json
"""

def is_integer(str):
    try:
        int(str)
        return True
    except ValueError:
        return False

def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

def extract_fields():
    data = {} #yaml.load(initial_input)
    data['format'] = {'version' : '0.1'}
    data['bibliography'] = [{'url' : 'https://www.nwchem-sw.org'}]
    data['generator'] = {'source' : 'nwchem',
                         'version' : '6.8'}
    skip_input_geometry = False
    geometry = None
    coulomb_repulsion = None
    scf_energy = None
    scf_energy_offset = None
    energy_offset = None
    one_electron_integrals = None
    two_electron_integrals = None
    n_electrons_alpha = None
    basis_set = None
    n_electrons_beta = None
    n_orbitals = None
    initial_state = None
    ccsd_energy = None
    reader_mode = ""
    excited_state_count = 1
    excitation_energy = 0.0
    fci_energy = None
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
            elif ln.find("Number of active alpha electrons") != -1:
                n_electrons_alpha = int(ln_segments[-1])
            elif ln.find("Number of active beta electrons") != -1:
                n_electrons_beta = int(ln_segments[-1])
            elif ln.find("Number of active orbitals") != -1:
                n_orbitals = int(ln_segments[-1])
            elif ln.find("Total MCSCF energy") != -1:
                fci_val = float(ln.split("=")[1].strip())
                fci_energy = {"units" : "hartree", "value" : fci_val, "upper": fci_val+0.1, "lower":fci_val-0.1}
            #elif ln.split("=")[0].strip() == "CCSD total energy / hartree":
            #    if fci_energy is None:
            #        fci_val = float(ln.split("=")[1].strip())
            #        fci_energy = {"units" : "hartree", "value" : fci_val, "upper": fci_val*0.99, "lower":fci_val*1.01}
            elif ln == "Ground state specification:":
                reader_mode = "initial_state"
                assert ccsd_energy is not None, "CCSD energy should be available before the groud state specification."
                if initial_state is None:
                    initial_state = []
                initial_state += [ {'state':{
                        'label' : '|G>',
                        'energy' : { 'units' : 'hartree', 
                                     'value' : ccsd_energy} ,
                        'superposition' : []
                        }}]
            elif ln == "Excited state specification:":
                reader_mode = "initial_state"
                assert ccsd_energy is not None, "CCSD energy should be available before the excited state specification."
                if initial_state is None:
                    initial_state = []
                initial_state += [ {'state':{
                        'label' : '|E%d>'%(excited_state_count),
                        'energy' : { 'units' : 'hartree', 
                                     'value' : ccsd_energy+excitation_energy} ,
                        'superposition' : []
                        }}]
                excited_state_count += 1
            elif ln.split("=")[0].strip() == "Excitation energy / hartree":
                excitation_energy = float(ln.split("=")[1].strip())
            elif ln.split("=")[0].strip() == "CCSD total energy / hartree":
                ccsd_energy = float(ln.split("=")[1].strip())
            elif ln == '''Geometry "geometry" -> ""''':
                reader_mode = "cartesian_geometry"
                if geometry is None:
                    geometry = {'coordinate_system': 'cartesian'}
                geometry['atoms'] = []
        elif reader_mode == "cartesian_geometry":
            if ln_segments[0] == "Output":
                if 'angstroms' in ln_segments:
                    geometry['units'] = 'angstrom'
                elif 'a.u.' in ln_segments:
                    geometry['units'] = 'bohr'
            elif ln_segments[0] == 'Atomic':
                reader_mode = ""
            elif is_integer(ln_segments[0]):
                assert 'atoms' in geometry
                #if not 'atoms' in geometry:
                #    geometry['atoms'] = []
                geometry['atoms'] += [{"name":ln_segments[1],
                                       "coords":
                                       [float(ln_segments[3]), float(ln_segments[4]), float(ln_segments[5])]}]
        elif reader_mode == "initial_state":
            if ln == "-------------------------------------":
                reader_mode = ""
            elif 'norm' in ln or 'string' in ln or 'exp' in ln:
                continue
            else:
                #vals = ln.split(":")
                if len(initial_state[-1]['state']['superposition']) > 0 and initial_state[-1]['state']['superposition'][-1][-1] != '|vacuum>':
                    initial_state[-1]['state']['superposition'][-1] += ln.replace('|0>', '|vacuum>').split()
                else:
                    vals = ln.replace('|0>', '|vacuum>').split(":")
                    initial_state[-1]['state']['superposition'] += [
                        [float(vals[0].strip())] +
                        vals[1].split()[:]
                    ]
                    '''initial_state[-1]['state']['superposition'] += [
                        [float(vals[0].strip())] +
                        vals[1].split()[:-1] +
                        ['|vacuum>']
                    ]'''
        elif reader_mode == "input_deck":
            if ln == "================================================================================":
                reader_mode = ""
            elif ln_segments[0:2]== ["geometry", "units"]:
                reader_mode = "input_geometry"
                assert len(ln_segments) >= 3
                geometry = {'coordinate_system': 'cartesian'}
                geometry['units'] = ln_segments[2]
            elif ln_segments[0].lower() == "basis":
                reader_mode = "input_basis"
                basis_set = {}
                basis_set['name'] = 'unknown'
                basis_set['type'] = 'gaussian'
        elif reader_mode=="input_geometry":
            if ln.lower()=="end":
                reader_mode = "input_deck"
            elif ln_segments[0] == "symmetry":
                assert len(ln_segments) == 2
                geometry['symmetry'] = ln_segments[1]
            elif skip_input_geometry == False: #atom description line
                if len(ln_segments) != 4 or not is_float(ln_segments[1]) or not is_float(ln_segments[2]) or not is_float(ln_segments[3]):
                    skip_input_geometry = True
                else:
                    if not 'atoms' in geometry:
                        geometry['atoms'] = []
                        geometry['atoms'] += [{"name":ln_segments[0],
                                               "coords":
                                                   [float(ln_segments[1]), float(ln_segments[2]), float(ln_segments[3])]}]
        elif reader_mode == "input_basis":
            if ln.lower() == "end":
                reader_mode = "input_deck"
            else:
                if ln.find('library') != -1:
                    assert len(ln_segments) == 3
                    assert ln_segments[1] == "library"
                    basis_set['name'] = ln_segments[2]
                    basis_set['type'] = 'gaussian'
        elif reader_mode == "one_electron_integrals":
            if ln == "end_one_electron_integrals":
                reader_mode = ""
            else:
                assert len(ln_segments) == 3
                if int(ln_segments[0]) >= int(ln_segments[1]):
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
    if fci_energy is None:
        assert ccsd_energy is not None
        fci_energy = {"units" : "hartree", "value" : ccsd_energy, "upper": ccsd_energy+0.1, "lower": ccsd_energy-0.1}

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
    if initial_state is not None:
        integral_sets[-1]["initial_state_suggestions"] = initial_state
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

    

def main():
    assert emitter_yaml or emitter_ruamel, "Extraction failed: could not import YAML or RUAMEL packages."
    if emitter_yaml:
        emitter_yaml_func()
    elif emitter_ruamel:
        emitter_ruamel_func()
    else:
        assert False, "Unreachable code"

if __name__ == "__main__":
    main()
