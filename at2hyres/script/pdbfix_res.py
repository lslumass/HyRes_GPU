import sys, os

## this script is used to fix the unstandard PDB file including:
## 1) change the chain name to 'X', regardless of original ones
## 2) re-order the resid starting from 1 for each segment
## 3) rename the ACE and CT3 terminus as AMN and CBX
## 4) set Occupancy (55-60) as 1.00 and set Temperature factor (61-66) as 0.00

## Author: Shanlong Li @ UMass
## Latest: Jan 11, 2024

pdb_file = sys.argv[1]  
out_file = sys.argv[2]
capping = sys.argv[3]

## read all atoms
with open(pdb_file, "r") as pdb:
    lines = pdb.readlines()

## rename the chain name
atoms = []
for line in lines:
    if line.startswith("ATOM"):
        new = line[:21] + 'X' + line[22:]
        atoms.append(new)

## get the length of each chain and number of residue of each chain
segname = atoms[0].split()[-1]
chains = []
ress = []
count = 0
resnum = 1
resid = int(atoms[0][22:26])
for atom in atoms:
    if atom.split()[-1] == segname:
        count += 1
        if int(atom[22:26]) != resid:
            resnum += 1
            resid = int(atom[22:26])
    else:
        chains.append(count)
        ress.append(resnum)
        segname = atom.split()[-1]
        count = 1
        resnum = 1
        resid = int(atom[22:26])
chains.append(count)
ress.append(resnum)
#print(chains, ress)

if capping == '1':
    ## re-number the resid of N-T/C-T
    count = 0
    for i in range(len(chains)):
        nt_start, ct_start = count, count + chains[i] - 6
        nt_end, ct_end = nt_start + 6, ct_start + 6
        for j in range(nt_start, ct_end):
            if j < nt_end:
                atoms[j] = atoms[j][:17] + 'AMN' + atoms[j][20:22] + '   1' + atoms[j][26:]
            elif j >= ct_start:
                atoms[j] = atoms[j][:17] + 'CBX' + atoms[j][20:22] + '   2' + atoms[j][26:]
            else:
                atoms[j] = atoms[j][:22] + str(int(atoms[j][22:26]) + 2).rjust(4) + atoms[j][26:] 
        count += chains[i]
else:
    print('no cap')
    ## re-name the Oxygen of C=O to "O"
    count = 0
    for i in range(len(chains)):
        first_atom = count
        last_atom = count + chains[i] - 1
        for j in range(last_atom, first_atom, -1):
            if str(atoms[j][12:16]) == ' C  ':
                atoms[j+1] = atoms[j+1][:12] + ' O  ' + atoms[j+1][16:]
                break
        count += chains[i]

n_res = 1
resid = int(atoms[0][22:26])
segid = atoms[0].split()[-1]
new = []
for atom in atoms:
    if atom.split()[-1] == segid:
        if int(atom[22:26]) != resid:
            resid = int(atom[22:26])
            n_res += 1
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
        else:
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
    else:
        segid = atom.split()[-1]
        resid = int(atom[22:26])
        n_res = 1
        if int(atom[22:26]) != resid:
            resid = int(atom[22:26])
            n_res += 1
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
        else:
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
    new.append(l)

with open(out_file, "w") as out:
    for line in new:
        print(line, file=out)
    print('END', file=out)
