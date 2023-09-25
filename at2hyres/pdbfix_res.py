import sys

# fix the resid for different chain in PDB file

pdb_file = sys.argv[1]  
out_file = sys.argv[2]
capping = sys.argv[3]

## read all atoms
with open(pdb_file, "r") as pdb:
    lines = pdb.readlines()

atoms = []
for line in lines:
    if line.startswith("ATOM"):
        atoms.append(line)

if capping == '1':
    ## get the length of each chain
    segname = atoms[0].split()[-1]
    chains = []
    count = 0
    for atom in atoms:
        if atom.split()[-1] == segname:
            count += 1
        else:
            chains.append(count)
            segname = atom.split()[-1]
            count = 1
    chains.append(count)
    print(chains)
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
                atoms[j] = atoms[j][:22] + str(int(atoms[j].split()[5]) + 2).rjust(4) + atoms[j][26:] 
        count += chains[i]
else:
    print('no cap')

n_res = 1
resid = atoms[0].split()[5]
segid = atoms[0].split()[-1]
new = []
for atom in atoms:
    if atom.split()[-1] == segid:
        if atom.split()[5] != resid:
            resid = atom.split()[5]
            n_res += 1
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
        else:
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
    else:
        segid = atom.split()[-1]
        resid = atom.split()[5]
        n_res = 1
        if atom.split()[5] != resid:
            resid = atom.split()[5]
            n_res += 1
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
        else:
            l = atom[:22] + str(n_res).rjust(4) + atom[26:54] + str(1.00).rjust(6) + str(0.00).rjust(6) + atom[66:-1]
    new.append(l)

with open(out_file, "w") as out:
    for line in new:
        print(line, file=out)
    print('END', file=out)