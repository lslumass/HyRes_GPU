import sys
from psfgen import PsfGen
import argparse
import os

# Global varibale
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
hyres_topology_input = os.path.join(project_dir, "hyres_forcefield", "top_hyres_GPU.inp")

def main():
    parser = argparse.ArgumentParser(description="generate PSF for Hyres systems",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_pdb_files", help="Hyres PDB file(s)", required=True, nargs="+")
    parser.add_argument("-o", "--output_psf_file", help="output name/path for Hyres PSF", required=True, default="output.psf")
    parser.add_argument("-n", "--num_of_chains", help="Number of copies for each pdb; it should have the same length as the given pdb list specified in the '-i' argument", default=[1,], nargs="+")
    parser.add_argument("-t", "--ter", choices=['neutral', 'charged'], help="Terminal charged status (choose from ['neutral', 'charged'])", default='neutral')
    args = parser.parse_args()
   
    pdb_list = args.input_pdb_files
    outpsf = args.output_psf_file
    num_list = [1,]*len(pdb_list) if len(pdb_list) > 1 and args.num_of_chains == [1,] else [int(i) for i in args.num_of_chains] 
    assert len(pdb_list) == len(num_list), "pdb file list must have the same length as the chain number list (specified by the '-n' argument)"
    ter = args.ter

    gen = PsfGen()
    gen.read_topology(hyres_topology_input)


    if len(pdb_list) == 1:  # copies of singe chain 
        pdb = pdb_list[0]
        num = num_list[0]
        for i in range(num):
            segid = f"A{i}"
            gen.add_segment(segid=segid, pdbfile=pdb)
            if ter == 'charged':
                res_start, res_end = gen.get_resids(segid)[0], gen.get_resids(segid)[-1]
                gen.set_charge(segid, res_start, "N", 1.00)
                gen.set_charge(segid, res_end, "O", -1.00)
        gen.write_psf(filename=outpsf)

    else:   # combine different pdb files with different copies
        for idx, (pdb, num) in enumerate(zip(pdb_list, num_list), 1):
            # loop through each pdb and make copies
            for i in range(num):
                segid = f"{chr(64+idx)}{i}" 
                gen.add_segment(segid=segid, pdbfile=pdb)
                if ter == 'charged':
                    res_start, res_end = gen.get_resids(segid)[0], gen.get_resids(segid)[-1]
                    gen.set_charge(segid, res_start, "N", 1.00)
                    gen.set_charge(segid, res_end, "O", -1.00)
        gen.write_psf(filename=outpsf)

if __name__ == '__main__':
    main()
