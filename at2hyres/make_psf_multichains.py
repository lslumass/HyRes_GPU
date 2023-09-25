#!/usr/bin/env python
import sys, math
import numpy as np


## this script is used for creating psf file of multichains from psf of single chain
## author: Shanlong Li
## create date: Jan 29, 2023
## modify date: Sep 21, 2023

# define variables
psf_file = sys.argv[1]
num_chain = int(sys.argv[2])
psf_out = str(num_chain)+'.psf'
name = sys.argv[3]

# read the psf file of one protein
with open(psf_file, "r") as psf:
    lines = psf.readlines()

title = lines[0]
count = 0
for line in lines:
    l = line.split()
    if "!NATOM" in l:
        num_atom = int(l[0])
        atom_start = count + 1
        atom_end = atom_start + num_atom
    elif "!NBOND:" in l:
        num_bond = int(l[0])
        bond_start = count + 1
        bond_end = bond_start + math.ceil(num_bond/4)
    elif "!NTHETA:" in l:
        num_angle = int(l[0])
        angle_start = count + 1
        angle_end = angle_start + math.ceil(num_angle/3)
    elif '!NPHI:' in l:
        num_dihed = int(l[0])
        dihed_start = count + 1
        dihed_end = dihed_start + math.ceil(num_dihed/2)
    elif '!NIMPHI:' in l:
        num_impr = int(l[0])
        impr_start = count + 1
        impr_end = impr_start + math.ceil(num_impr/2)
    elif '!NDON:' in l:
        num_donor = int(l[0])
        donor_start = count + 1
        donor_end = donor_start + math.ceil(num_donor/4)
    elif '!NACC:' in l:
        num_accp = int(l[0])
        accp_start = count + 1
        accp_end = accp_start + math.ceil(num_accp/4)
    elif '!NNB' in l:
        num_NNB = num_atom
    elif '!NCRTERM:' in l:
        num_cross = int(l[0])
        cross_start = count + 1
        cross_end = cross_start + num_cross
    elif '!NGRP' in l:
        num_grp = int(l[0])
        grp_start = count + 1
        grp_end = grp_start + math.ceil(num_grp/3)
    count += 1

with open(psf_out, 'w') as f:
    print(title, file=f)
    print('       2 !NTITLE', file=f)
    print('* PSF OF SIMULATION BOX', file=f)
    print('* CREATE THRUGH MAKE_PSF.PY Shanlong Li\n', file=f)
    print(str(num_atom*num_chain).rjust(8), '!NATOM', file=f)
    # atom section
    atoms = []
    for line in lines[atom_start:atom_end]:
        l = line.split()
        one = ' ' + l[2].ljust(5) + l[3].ljust(5) + l[4].ljust(5) + l[5].rjust(4) + l[6].rjust(11) + l[7].rjust(14) + l[8].rjust(12)
        atoms.append(one)
    for i in range(num_chain):
        segment = name + str(i)
        for j in range(num_atom):
            idx = i*num_atom + j + 1
            print(str(idx).rjust(8), segment.ljust(4) + atoms[j], file=f)
    # bond section
    print('',file=f)
    print(str(num_bond*num_chain).rjust(8), '!NBOND: bonds', file=f)
    bond, bonds = [], []
    for line in lines[bond_start:bond_end]:
        l = iter(line.split())
        for idx in l:
            bond.append([int(idx), int(next(l))])
    for i in range(num_chain):
        for b in bond:
            idx1, idx2 = i*num_atom + b[0], i*num_atom + b[1]
            bonds.append(str(idx1).rjust(8) + str(idx2).rjust(8))
    for i in range(0, len(bonds), 4):
        bond = bonds[i:i+4]
        print(''.join([b for b in bond]), file=f)
    # angle section
    print('',file=f)
    print(str(num_angle*num_chain).rjust(8), '!NTHETA: angles', file=f)
    ang, angs = [], []
    for line in lines[angle_start:angle_end]:
        for i in range(0, len(line.split()), 3):
            l = line.split()[i:i+3]
            ang.append([int(l[0]), int(l[1]), int(l[2])])
    for i in range(num_chain):
        for a in ang:
            idx1, idx2, idx3 = i*num_atom + a[0], i*num_atom + a[1], i*num_atom + a[2]
            angs.append(str(idx1).rjust(8) + str(idx2).rjust(8) + str(idx3).rjust(8))
    for i in range(0, len(angs), 3):
        ang = angs[i:i+3]
        print(''.join([a for a in ang]), file=f)

    # dihedral section
    print('',file=f)
    print(str(num_dihed*num_chain).rjust(8), '!NPHI: dihedrals', file=f)
    dihed, diheds = [], []
    for line in lines[dihed_start:dihed_end]:
        for i in range(0, len(line.split()), 4):
            l = line.split()[i:i+4]
            dihed.append([int(l[0]), int(l[1]), int(l[2]), int(l[3])])
    for i in range(num_chain):
        for d in dihed:
            idx1, idx2, idx3, idx4 = i*num_atom + d[0], i*num_atom + d[1], i*num_atom + d[2], i*num_atom + d[3]
            diheds.append(str(idx1).rjust(8) + str(idx2).rjust(8) + str(idx3).rjust(8) + str(idx4).rjust(8))
    for i in range(0, len(diheds), 2):
        dihed = diheds[i:i+2]
        print(''.join([d for d in dihed]), file=f)

    # improper section
    print('',file=f)
    print(str(num_impr*num_chain).rjust(8), '!NIMPHI: impropers', file=f)
    impr, imprs = [], []
    for line in lines[impr_start:impr_end]:
        for i in range(0, len(line.split()), 4):
            l = line.split()[i:i+4]
            impr.append([int(l[0]), int(l[1]), int(l[2]), int(l[3])])
    for i in range(num_chain):
        for im in impr:
            idx1, idx2, idx3, idx4 = i*num_atom + im[0], i*num_atom + im[1], i*num_atom + im[2], i*num_atom + im[3]
            imprs.append(str(idx1).rjust(8) + str(idx2).rjust(8) + str(idx3).rjust(8) + str(idx4).rjust(8))
    for i in range(0, len(imprs), 2):
        impr = imprs[i:i+2]
        print(''.join([im for im in impr]), file=f)
    
    # donor section
    print('',file=f)
    print(str(num_donor*num_chain).rjust(8), '!NDON: donors', file=f)
    donor, donors = [], []
    for line in lines[donor_start:donor_end]:
        l = iter(line.split())
        for idx in l:
            donor.append([int(idx), int(next(l))])
    for i in range(num_chain):
        for d in donor:
            idx1, idx2 = i*num_atom + d[0], i*num_atom + d[1]
            donors.append(str(idx1).rjust(8) + str(idx2).rjust(8))
    for i in range(0, len(donors), 4):
        donor = donors[i:i+4]
        print(''.join([d for d in donor]), file=f)

    # acceptor section
    print('',file=f)
    print(str(num_accp*num_chain).rjust(8), '!NACC: acceptors', file=f)
    accp, accps = [], []
    for line in lines[accp_start:accp_end]:
        l = iter(line.split())
        for idx in l:
            accp.append([int(idx), int(next(l))])
    for i in range(num_chain):
        for a in accp:
            idx1, idx2 = i*num_atom + a[0], i*num_atom + a[1]
            accps.append(str(idx1).rjust(8) + str(idx2).rjust(8))
    for i in range(0, len(accps), 4):
        accp = accps[i:i+4]
        print(''.join([a for a in accp]), file=f)
    
    # NNB section
    print('',file=f)
    print('       0 !NNB\n', file=f)
    mod, div = num_atom*num_chain//8, num_atom*num_chain%8
    for i in range(mod):
        print('       0       0       0       0       0       0       0       0', file=f)
    nnb = ''
    for j in range(div):
        nnb += '       0'
    print(nnb, file=f)

    # GRP setction
    print('',file=f)
    print(str(num_grp*num_chain).rjust(8), '      0 !NGRP NST2', file=f)
    grp, grps = [], []
    for line in lines[grp_start:grp_end]:
        for i in range(0, len(line.split()), 3):
            l = line.split()[i:i+3]
            grp.append([int(l[0]), l[1].rjust(8) + l[2].rjust(8)])
    for i in range(num_chain):
        for g in grp:
            idx, cont = i*num_atom + g[0], g[1]
            grps.append(str(idx).rjust(8) + cont)
    for i in range(0, len(grps), 3):
        grp = grps[i:i+3]
        print(''.join([g for g in grp]), file=f)
    
    # MOLNT
    print('',file=f)
    print(str(num_chain).rjust(8), '!MOLNT', file=f)
    mlts = []
    for i in range(num_chain):
        for j in range(num_atom):
            mlts.append(str(i+1).rjust(8))
    for i in range(0, len(mlts), 8):
        mlt = mlts[i:i+8]
        print(''.join([m for m in mlt]), file=f)
    
    print('',file=f)
    print('       0       0 !NUMLP NUMLPH', file=f)

    # cros-term section
    print('',file=f)
    print(str(num_cross*num_chain).rjust(8), '!NCRTERM: cross-terms', file=f)
    ct, cts = [], []
    for line in lines[cross_start:cross_end]:
        l = line.split()
        ct.append([int(x) for x in l])
    for i in range(num_chain):
        for c in ct:
            one = ''
            for x in c:
                one += str(x+i*num_atom).rjust(8)
            print(one, file=f)
    
print('Done!')
