# Prepare Hyres PSF

we created a `psfgen_hyres.py` in the **generate_psf** folder to generate Hyres PSF for all kinds of scenarios.

```
# see help
python psfgen_hyres.py -h

generate PSF for Hyres systems

optional arguments:
  -h, --help            Show this help message and exit
  -i INPUT_PDB_FILES    Hyres PDB file(s) of mononer (default: None)
  -o OUTPUT_PSF_FILE    Output name/path for Hyres PSF (default: output.psf)
  -n NUM_OF_CHAINS      Number of copies for each pdb (default: [1,])
  -t {neutral,charged}  Terminal charged status (default: neutral)

```
1. generate PSF for a single-chain protein
```
python psfgen_hyres.py -i protein.pdb -o protein.psf
python psfgen_hyres.py -i protein.pdb -o protein.psf -ter charged
```
2. generate PSF for a multi-chain protein
```
python psfgen_hyres.py -i chainA.pdb chainB.pdb -o complex.psf
python psfgen_hyres.py -i chainA.pdb chainB.pdb -o complex.psf -ter charged
```
3. generate PSF for LLPS simulation of one kind of protein
```
python psfgen_hyres.py -i chainA.pdb -n 100 -o llps.psf
python psfgen_hyres.py -i chainA.pdb -n 100 -o llps.psf -ter charged
```
4. generate PSF for LLPS simulation of multiple kinds of protein
```
# 20 copies of chain A + 30 copies of chain B
python psfgen_hyres.py -i chainA.pdb chainB.pdb -n 20 30 -o llps.psf
python psfgen_hyres.py -i chainA.pdb chainB.pdb -n 20 30 -o llps.psf -ter charged
```

