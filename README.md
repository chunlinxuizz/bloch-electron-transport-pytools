# TE_Calc_Toolkit
- Small codes to treat [CP2K](https://www.cp2k.org/)
outputs for transport calculations using [AMSET](https://github.com/hackingmaterials/Amset)
- They are still under development, any cooperation or suggestions are welcome😊
- `my_cp2k_packages.cp2k_to_pymatgen` can parse CP2K input and output files and initialize
   useful objects of [pymatgen](https://pymatgen.org/), such as
  `<pymatgen.electronic_structure.Bandstructure>`
- `cp2k_transport.cp2k_elastic_constant` can fit and plot the elastic constant in a specific direction 
- `cp2k_transport.cp2k_deformation_potential` can fit and plot the deformation potential of VBM, CBM and
  Fermi level for metallic systems using core levels or vacuum level as reference
- A pymatgen `Bandstructure` object is needed to initialize the amset `Runner`, and to obtain accurate bandstructure,
  accurate eigenvalues based on an ultra dense k-mesh is necessary before interpolation. However, it is memory intensive to
  calculate them in a single calculation. The purpose of `cp2k_transport.kmesh_to_cp2kinp` and `cp2k_transport.merge_cp2k_band` are
  spread the dense Monkhost-Pack k-mesh into several input files and merge the .bs files in to one, respectively.
- The elastic constants, deformation potentials, as well as bandstructures are the input of AMSET to get thermoelectric transport
  coefficients(https://hackingmaterials.lbl.gov/amset/)
