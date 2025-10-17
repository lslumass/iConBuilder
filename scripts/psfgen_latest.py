from psfgen import PsfGen
from iConBuilder import utils
import argparse, os
import MDAnalysis as mda


def decompose_complex(pdb, idx, i, gen):
    u = mda.Universe(pdb)
    segids = u.residues.segments.segids
    segnum = len(segids)
    for j, segid in enumerate(segids):
        sel = u.select_atoms(f"segid {segid}")
        tmp_pdb = f'tmp_{segid}.pdb'
        sel.atoms.write(tmp_pdb)
        if segid.startswith("P"):
            segid = f"PC{chr(65+idx)}{j+i*segnum}"
            gen.add_segment(segid=segid, pdbfile=tmp_pdb, auto_angles=False)
        elif segid.startswith("R"):
            segid = f"RC{chr(65+idx)}{j+i*segnum}"
            gen.add_segment(segid=segid, pdbfile=tmp_pdb, auto_angles=False, auto_dihedrals=False)
        else:
            print("Error: Only protein-protein or protein-RNA complex is supported.")
            exit(1)

def set_terminus(gen, segid, charge_status):
    if segid.startswith("P"):
        nter, cter = gen.get_resids(segid)[0], gen.get_resids(segid)[-1]
        if charge_status == 'charged':
            gen.set_charge(segid, nter, "CA", 1.00)
            gen.set_charge(segid, cter, "CO", -1.00)
        elif charge_status == 'NT':
            gen.set_charge(segid, nter, "CA", 1.00)
        elif charge_status == 'CT':
            gen.set_charge(segid, cter, "CO", -1.00)
        elif charge_status == 'positive':
            gen.set_charge(segid, nter, "CA", -1.00)
            gen.set_charge(segid, cter, "CO", -1.00)
        else:
            print("Error: Only 'neutral', 'charged', 'NT', and 'CT' charge status are supported.")
            exit(1)


def main():
    parser = argparse.ArgumentParser(description="generate PSF for Hyres systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_pdb_files",  help="Hyres PDB file(s), it should be the pdb of monomer", required=True, nargs="+")
    parser.add_argument("-o", "--output_psf_file", help="output name/path for Hyres PSF", required=True)
    parser.add_argument("-n", "--num_of_chains", type=int, default=[1,], nargs="+", 
                        help="Number of copies for each pdb; it should have the same length as the given pdb list specified in the '-i' argument")
    parser.add_argument("-m", "--molecule_type", type=str, default=['P'], nargs="+",
                        help="select from 'P', 'R', 'D', 'C', 'Mg' for protein, RNA, DNA, complex, Magnesium, respectively")
    parser.add_argument("-t", "--ter", choices=['neutral', 'charged', 'NT', 'CT', 'positive'], 
                        help="Terminal charged status (choose from ['neutral', 'charged', 'NT', 'CT', 'positive'])", default='neutral')
    args = parser.parse_args()
   
    pdb_list = args.input_pdb_files
    outpsf = args.output_psf_file
    num_list = args.num_of_chains
    
    assert len(pdb_list) == len(num_list), "chain number list must have the same length as the pdb list (specified by the '-n' argument)"
    type_list = args.molecule_type
    assert len(pdb_list) == len(type_list), "molecule type list must have the same length as the pdb list (specified by the '-n' argument)"
    ter = args.ter

    # load topology files
    RNA_topology, _ = utils.load_ff('RNA')
    protein_topology, _ = utils.load_ff('Protein')

    # generate psf
    gen = PsfGen()
    gen.read_topology(RNA_topology)
    gen.read_topology(protein_topology)

    for idx, pdb in enumerate(pdb_list): 
        # loop through each pdb and make copies
        num = num_list[idx]
        mol_type = type_list[idx]
        if mol_type == 'P':
            for i in range(num):
                segid = f"P{chr(65+idx)}{i}" 
                gen.add_segment(segid=segid, pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
        elif mol_type == 'R':
            for i in range(num):
                segid = f"R{chr(65+idx)}{i}" 
                gen.add_segment(segid=segid, pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
        elif mol_type == 'C':
            for i in range(num):
                decompose_complex(pdb, idx, i, gen)
        elif mol_type == 'Mg':
            for i in range(num):
                segid = f"MG{i}" 
                gen.add_segment(segid=segid, pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
        else:
            print("Error: Only type of 'P', 'R', 'D', 'C', and 'Mg' are supported.")
            exit(1)

    # re-set the charge status of terminus
    for segid in gen.get_segids():
        if ter != "neutral":
            set_terminus(gen, segid, ter)       
            
    gen.write_psf(filename=outpsf)
    os.system("rm -rf tmp_*.pdb")

if __name__ == '__main__':
    main()
