#!/usr/bin/env python
import sys
import numpy as np
from my_cp2k_packages.cp2k_to_pymatgen import Cp2kBand
from pymatgen.electronic_structure.core import Spin


# the structure of EIGENVAL file:
# 1) # of ions , # of ions, the number of loops after which the averaged pair correlation functions and the DOS are written, ISPIN
# 2) the volume of the cell (in Ã…^3) and the lattice parameters of the box (in m)
# 3) T
# 4) the string 'CAR'
# 5) the header
# 6) # of electrons, # of k-points, #of bands
# 7) the k-point and its weight
try:
    inp_file = sys.argv[1]
    bs_file = sys.argv[2]
except IndexError:
    print("\nusage: bs_to_eigenval.py [inp_filename] [bs_filename]\n")
    sys.exit()

system = Cp2kBand()
bandstructure = system.get_band_structure(inp_file, bs_file)
n_electrons = system.nelectron
volume = bandstructure.structure.lattice.volume
a, b, c = bandstructure.structure.lattice.abc
n_ion = len(bandstructure.structure.sites)
ISPIN = 1 # the spin polarization is not included at this stage
T = 1.000000000000000E-004
header = "generated from cp2k.bs"
n_bands, n_kpoints = np.shape(bandstructure.bands[Spin(1)])
kpoints = bandstructure.kpoints
eigenvals = np.transpose(bandstructure.bands[Spin(1)])  # shape: (n_kpoints, n_bands)

kpoint_style = "    {:>17.14f}    {:>17.14f}    {:>17.14f}        {:>12.6e}\n"
eigenval_style = "{:>5d}      {:>12.8f}\n"

with open("EIGENVAL", 'w') as fout:
    fout.write(f"   {n_ion}   {n_ion}     1     {ISPIN}\n")
    fout.write("  {:>13.7e}  {:>13.7e}  {:>13.7e}  {:>13.7e}  {:>13s}\n".format(
        volume, a, b, c, str(0.5000000E-15)
        ))
    fout.write("  1.000000000000000E-004\n")
    fout.write("  CAR\n")
    fout.write(f" {header}\n")
    fout.write(f"  {n_electrons}  {n_kpoints}  {n_bands}\n")
    fout.write("\n")

    for k in range(n_kpoints):
        fout.write(kpoint_style.format(
            kpoints[k].a, kpoints[k].b, kpoints[k].c, 1/n_kpoints
        ))
        for band in range(n_bands):
            fout.write(eigenval_style.format(
                band+1, eigenvals[k][band]
            ))
        fout.write("\n")



