from my_cp2k_packages.cp2k_to_pymatgen import Cp2kBand
from amset.plot.electronic_structure import ElectronicStructurePlotter
from amset.tools.plot import get_kpath, _log_band_stats
import warnings

warnings.simplefilter("ignore")

inpfile = "cp2k.inp"
bsfile  = "cp2k.bs"
job = Cp2kBand()
bandstructure = job.get_band_structure(inpfile, bsfile)
structure = job.structure
nelectron = job.nelectron

#kpt_list =  [ [[0.,0.5, 0.], [0., 0., 0]],[[0., 0., 0], [0.5, 0., 0.]] ]
#labels = [ ['Y', r'$\Gamma$'], [r'$\Gamma$','X'] ]

kpt_list =  [ [[0.,0., 0.5], [0., 0., 0]],[[0., 0., 0], [0., 0.5, 0.]] ]
labels = [ ['Z', r'$\Gamma$'], [r'$\Gamma$','Y'] ]

kpath = get_kpath(structure, mode = "seekpath", kpt_list = kpt_list, labels = labels)

plotter = ElectronicStructurePlotter(bandstructure, nelectron, interpolation_factor=5)
plt, bs_plotter = plotter.get_plot(
    plot_band_structure=True, 
    plot_dos=True, 
    vbm_cbm_marker=False, 
    zero_to_efermi=True,
    kpath = kpath, 
    emin=-2., 
    emax=1.,
    dos_aspect=3,
    dos_estep=0.01,
    dos_kpoints=50,
    line_density=100,
    height=6,
    width=None,
    elabel="Energy (eV)",
    style='~/xcl/bin/band_style.mplstyle',
    #style=None,
    no_base_style=True,
    fonts=None,
    return_plotter=True, 
)

plt.savefig("band.png", bbox_inches="tight", dpi=400)

_log_band_stats(bs_plotter._bs[0])
