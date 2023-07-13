#!/usr/bin/env python
import os
import numpy as np
from my_cp2k_packages.cp2k_to_pymatgen import Cp2kBand
from pymatgen.io.cp2k.outputs import Cp2kOutput
import matplotlib.pyplot as plt

pwd = os.getcwd()
strains = ["-0.010","-0.005","0.000", "+0.005","+0.010"]

energies = np.zeros(len(strains))

def check_convergence_and_get_energy(filename):
    out = Cp2kOutput(filename)
    out.convergence()
    out.parse_energies()
    return out.final_energy

def get_volume(filename):
    out = Cp2kBand()
    out.parse_structure(filename)
    return out.structure.lattice.volume 

def save_fig(strains, energies, Cij):
    fit = np.polyfit(strains, energies,2)
    x = np.linspace(np.min(strains),np.max(strains))
    y = fit[0]*x**2 + fit[1]*x + fit[2]
    plt.figure(figsize = (7,6))
    plt.plot(x,y,'--',c = 'black')
    plt.scatter(strains,energies,c = 'b',s = 50)
    plt.xlabel("Strain")
    plt.xticks(strains)
    plt.ylabel("Energy (eV/cell)")
    plt.title("Elastic constant: {:.2f} GPa".format(Cij))
    name = str(os.getcwd()).split("/")[-1] + '.pdf'
    plt.savefig(name, format = 'pdf')

style = "    {:>6s}\t{:>10.6f}"
print("\n    {:>6s}\t{:>10s}".format("Strain","Energy (eV)"))
for i,strain in enumerate(strains):
    foldname = pwd + '/' + "strain_" + strain
    outfile = foldname + '/' + 'cp2k.out'
    inpfile = foldname + '/' + 'cp2k.inp'
    energies[i] = check_convergence_and_get_energy(outfile)
    if float(strain) == 0.0:
        base_energy = energies[i]
        volume = get_volume(inpfile)

    print(style.format(strain, energies[i]))

strains = np.array(strains, dtype = float)
energies = energies - base_energy

Cij = np.polyfit(strains, energies,2)[0] * 1.602176634E-19 / volume / 1E-30 / 1E9  # Gpa

print("\n    Elastic constant (Gpa): {:.2f}\n".format(Cij))

save_fig(strains, energies, Cij)
