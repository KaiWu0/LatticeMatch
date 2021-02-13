LatticeMatch is a program to search and create heterostructures. It match two 2D unit cells to construct a new supercell.
Its first fortran version is written by Wei Hu, I rewrite it by python and add some new functions.
Its first version is provided in the paper "Highly-efficient heterojunction solar cells based on two-dimensional tellurene and transition metal dichalcogenides".

Usage
1. Prepare two POSCARs (vasp structure file), we call POSCAR1 and POSCAR2.
2. Input "python match.py POSCAR1 POSCAR2". You can see something like:
1 -------------
mismatch: len: 1.850%, angle: 0.133%, total atoms: 105
 0  2 -7  3, lattice constant: 11.425 34.206, angle: 59.934, atoms: 42
 1  3 -6  3, lattice constant: 11.217 33.652, angle: 60.000, atoms: 63
2 -------------
mismatch: len: 1.850%, angle: 0.133%, total atoms: 105
 0  2 -7  3, lattice constant: 11.425 34.206, angle: 59.934, atoms: 42
 1 -2  9  3, lattice constant: 11.217 33.652, angle: 60.000, atoms: 63
3 -------------
mismatch: len: 1.850%, angle: 0.133%, total atoms: 105
 0  2 -7  3, lattice constant: 11.425 34.206, angle: 59.934, atoms: 42
 2  3 -3  6, lattice constant: 11.217 33.652, angle: 60.000, atoms: 63

If there is no output in this step, you can increase the max_atoms argument in match.py.
If there are too many outputs in this step, you can decrease the max_atoms argument in match.py.

3. Select a structure you like (lower mismatch, less total atoms), such the second. Input "python match.py POSCAR1 POSCAR2 2", you will get 4 POSCARs, their two layers are slide with respect to each other, you can choose any of one.
4. Verify the generated structure by VESTA or Materials Studio.
