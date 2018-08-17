 #Copyright (c) 2016 Kyle Barlow
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys
import os
import pybel # Not using this import, but available in case you want to use the alternate pybel interface into open babel
import openbabel as ob

"""
A Python script that calls the Confab rotamer generator (within OpenBabel).
The script allows for rotamers to be generated only for certain parts of the molecule
if the user specifies atoms that should be "unfrozen" (allowed to move) during conformer generation.
Script by: Kyle Barlow
Confab software:
O'Boyle NM, Vandermeersch T, Flynn CJ, Maguire AR, Hutchison GR. Confab - Systematic generation of diverse low-energy conformers. Journal of Cheminformatics. 2011;3:8. doi:10.1186/1758-2946-3-8.
Open Babel: http://openbabel.org/
"""

def read_molecules(filepath, single = False):
    in_format = filepath.strip().split( '.' )[-1]
    obconversion = ob.OBConversion()
    obconversion.SetInFormat( in_format )
    obmol = ob.OBMol()

    molecules = []

    notatend = obconversion.ReadFile( obmol, filepath )
    while notatend:
        molecules.append( obmol )
        obmol = ob.OBMol()
        notatend = obconversion.Read( obmol )

    if single:
        assert( len(molecules) == 1 )
        return molecules[0]
    else:
        return molecules

def count_aromatic_bonds( mol ):
    aromatic_bond_count = 0
    total_bond_count = 0
    for bond in ob.OBMolBondIter( mol ):
        total_bond_count += 1
        if bond.IsAromatic():
            aromatic_bond_count += 1
    print( '%d bonds are aromatic (out of %d total)' % ( aromatic_bond_count, total_bond_count ))


# Change confab settings if you want to
def generate_limited_rotamer_set(
        input_path, output_path, atoms_to_unfreeze,
        output_format = 'mol2',
        include_input_conformation_in_output = True,
        # Note: the Confab settings below (rmsd_cutoff, conf_cutoff, energy_cutoff, verbose)
        # are currently set to the default values from Confab.
        # You you might need to change these for your project.
        rmsd_cutoff = 0.5,
        conf_cutoff = 10000,
        energy_cutoff = 25.0,
        confab_verbose = True
):
    assert( output_path.endswith(output_format) )

    assert( os.path.isfile(input_path) )
    mol = read_molecules( input_path, single = True )
    # force field for open babel applied - mmff94
    pff = ob.OBForceField_FindType( "mmff94" )
    assert( pff.Setup(mol) ) # Make sure setup works OK

    count_aromatic_bonds( mol )

    print( "Molecule " + mol.GetTitle())
    print( "Number of rotatable bonds = " + str( mol.NumRotors() ))

    if len(atoms_to_unfreeze) > 0:
        print( '%d atoms will be allowed to move; freezing others' % len( atoms_to_unfreeze ))
        constraints = ob.OBFFConstraints()
        for atom in ob.OBMolAtomIter(mol):
            atom_id = atom.GetIndex() + 1
            if atom_id not in atoms_to_unfreeze:
                constraints.AddAtomConstraint(atom_id)
        pff.SetConstraints( constraints )

    # Run Confab conformer generation
    pff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, confab_verbose)

    pff.GetConformers(mol);
    if include_input_conformation_in_output:
        confs_to_write = mol.NumConformers()
    else:
        confs_to_write = mol.NumConformers() - 1

    print( "Generated %d conformers total (%d will be written)" % (mol.NumConformers(), confs_to_write))

    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)

    output_strings = []

    for conf_num in range(confs_to_write):
      mol.SetConformer(conf_num);
      output_strings.append( obconversion.WriteString(mol) )

    print( 'Writing %d conformations to %s' % (len(output_strings), output_path))
    with open(output_path, 'w') as f:
        for output_string in output_strings:
            f.write(output_string)

if __name__ == '__main__':
    # You can display atom numbers (IDs) in pymol with:
    # label all, " %s" % ID
    # This allows you to manually choose which atoms to "unfreeze"
    # Rotamers will only be generated if they involve only unfrozen atoms

    assert( len(sys.argv) >= 3 )
    input_path = sys.argv[1]
    assert( os.path.isfile(input_path) )

    output_path = sys.argv[2]
    output_format = output_path.strip().split('.')[-1]

    if len( sys.argv ) > 3:
        atoms_to_unfreeze = [int(x) for x in sys.argv[3:]]
    else:
        atoms_to_unfreeze = []

    generate_limited_rotamer_set( input_path, output_path, atoms_to_unfreeze, output_format = output_format )
