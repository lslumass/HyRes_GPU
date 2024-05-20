from openmm import *
import numpy as np
import mdtraj as md
from itertools import combinations

# TODO 
# define restraint functions for folded proteins
# mdtraj need to be the lastest version (>=1.9.9)

def apply_internal_restraints(system, pdb_file, domain_list=[], K_internal_cons = 10*unit.kilojoule_per_mole/unit.nanometers**2, cutoff=1.2):
    """
    system: omm system
    pdb_file: path to the pdb
    domain_list: [('A', (1, 132)), ('A', (280, 376)), ('B', (1, 132)), ...]
    K_internal_cons: K constant
    cutoff: applying internal restraints within this cutoff of a residue, 1.2 nm
    """
    pdb = md.load_pdb(pdb_file)
    internal_force = HarmonicBondForce()

    chainid_dict = {}
    for chain in pdb.topology.chains:
        chainid_dict[chain.chain_id] = chain.index

    if len(domain_list) == 0:
        domain_list = []
    else:
        for domain in domain_list:
            chainid, (starting_resid, ending_resid) = domain
            # change chainID to mdtraj int-based chainid
            chainID = chainid_dict[chainid]
        
            # slice: get type: mdtraj.core.trajectory.Trajectory
            selection = f"chainid {chainID} and residue {starting_resid} to {ending_resid}" 
            selected_domain = pdb.atom_slice(pdb.topology.select(selection))
            # calculate secondary structure
            dssp = md.compute_dssp(selected_domain, simplified=True)
            # use "1" for structured amino acid, and "0" for loops
            folded = np.where(dssp!='C', 1, 0)[0]
            
            # get the C-alpha atom index
            resid = [selected_domain.topology.residue(idx).resSeq for idx, i in enumerate(folded) if i==1 ]
            folded_CA_idx = [pdb.topology.select(f"chainid {chainID} and residue {i} and name CA ")[0] for i in resid]

            pairs = list(combinations(folded_CA_idx, 2))
            internal_pairs = np.array(pairs)
            pairs_num = 0
            for index in internal_pairs:
                r1=np.squeeze(pdb.xyz[:,int(index[0]),:])
                r2=np.squeeze(pdb
                _md.xyz[:,int(index[1]),:])
                dist0=np.sqrt((r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2)
                if dist0 < cutoff:
                    pairs_num += 1
                    internal_force.addBond(int(index[0]),int(index[1]),dist0*unit.nanometers, K_internal_cons)
            print(f"Number of internal pairs of resid {str(domain)}: {pairs_num}")
        system.addForce(internal_force)
