import gemmi
import argparse
import sys

def convert_pdb_to_cif(pdb_file, cif_file, entry_id):
    # Read the PDB file
    structure = gemmi.read_structure(pdb_file)
    
    # Create a new CIF document
    doc = gemmi.cif.Document()
    block = doc.add_new_block('converted_structure')
    
    # Add _entry.id
    block.set_pair('_entry.id', entry_id)
    
    # Add atom_site loop with the specified fields
    atom_site_loop = block.init_loop('_atom_site.', [
        'group_PDB',
        'id',
        'type_symbol',
        'label_atom_id',
        'label_alt_id',
        'label_comp_id',
        'label_asym_id',
        'label_entity_id',
        'label_seq_id',
        'pdbx_PDB_ins_code',
        'Cartn_x',
        'Cartn_y',
        'Cartn_z',
        'occupancy',
        'B_iso_or_equiv',
        'pdbx_formal_charge',
        'auth_seq_id',
        'auth_comp_id',
        'auth_asym_id',
        'auth_atom_id',
        'pdbx_PDB_model_num'
    ])
    
    # Populate the atom_site loop with atom data
    for model_index, model in enumerate(structure, start=1):
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_site_loop.add_row([
                        'ATOM',
                        str(atom.serial),
                        atom.element.name,
                        atom.name,
                        '.',
                        residue.name,
                        chain.name if chain.name else 'A',  # Set to 'A' if missing
                        '1',  # entity_id (assuming single entity)
                        str(residue.seqid.num),
                        '.', 
                        f'{atom.pos.x:.3f}',
                        f'{atom.pos.y:.3f}',
                        f'{atom.pos.z:.3f}',
                        f'{atom.occ:.2f}',
                        f'{atom.b_iso:.2f}',
                        '?' if atom.charge == 0 else str(atom.charge),  # formal_charge
                        str(residue.seqid.num),
                        residue.name,
                        chain.name if chain.name else 'A',  # Set to 'A' if missing
                        atom.name,
                        str(model_index)
                    ])
    
    # Add entity_poly loop
    entity_poly_loop = block.init_loop('_entity_poly.', [
        'entity_id',
        'type',
        'nstd_linkage',
        'nstd_monomer',
        'pdbx_strand_id'
    ])
    
    # Add data to entity_poly loop
    # Note: This is example data. You may need to adjust this based on your specific PDB file
    entity_poly_loop.add_row([
        '1',
        'polyribonucleotide',
        'no',
        'yes',
        'A'
    ])
    
    # Write the CIF file
    doc.write_file(cif_file)

def main():
    parser = argparse.ArgumentParser(description='Convert PDB file to CIF format')
    parser.add_argument('input', help='Input PDB file path')
    parser.add_argument('output', help='Output CIF file path')
    parser.add_argument('-id', '--entry_id', default='NULL', help='Entry ID for the structure (default: NULL)')
    args = parser.parse_args()

    try:
        convert_pdb_to_cif(args.input, args.output, args.entry_id)
        print(f"Successfully converted {args.input} to {args.output} with entry ID: {args.entry_id}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()