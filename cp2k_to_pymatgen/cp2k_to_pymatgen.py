import linecache, re, pathlib
import numpy as np
import typing

from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.core.periodic_table import Element


class Cp2kBand:
    """
    Get Structure and BandStructure from CP2K .inp and .bs file

    """
    
    def __init__(self):
        
    def parse_structure(self, Structure_File: str = "cp2k.inp"):
        self._Structure_File = Structure_File
        self.charge = 0

        lattice = []
        species = []
        coords = []
        pbcs = {
            "XYZ" : (1,1,1),
            "XY"  : (1,1,0),
            "YZ"  : (0,1,1),
            "XZ"  : (1,0,1),
            "X"   : (1,0,0),
            "Y"   : (0,1,0),
            "Z"   : (0,0,1),
        }

        with open(self._Structure_File,'r') as f:
            lines = f.readlines()
            for linendx, line in enumerate(lines):
                if re.search("&CELL", line):
                    for i in range(1,4):
                        xyz = lines[linendx + i].strip().split()[1:4]
                        lattice.append(list(map(float,xyz)))
                    pbc = pbcs[lines[linendx + 4].strip().split()[1]]
                elif re.search("&COORD", line):
                    atom = 1
                    while True:
                        if re.search("&END COORD", lines[linendx + atom]):
                            break
                        atom_coord = lines[linendx + atom].strip().split()
                        species.append(atom_coord[0])
                        coords.append(list(map(float,atom_coord[1:4])))
                        atom += 1
                    coords = np.array(coords)
                elif re.search("CHARGE", line):
                    self.charge = int(line.strip().split()[1])
                
                elif re.search("MULTIPLICITY", line):
                    self.multiplicity = int(line.strip().split()[1])
                
        self.lattice = Lattice(lattice, pbc)
        self.lattice_rec = self.lattice.reciprocal_lattice
        self.structure = Structure(
            self.lattice, 
            species, 
            coords, 
            charge = self.charge,
            validate_proximity = False,
            to_unit_cell = False,
            coords_are_cartesian = True,
            site_properties = None
        )
        
    def parse_bandstructure(self, Band_File: str = "cp2k.bs") -> BandStructure:
        self._Band_File = Band_File
        
        with open(self._Band_File,'r') as f:
            firstline = f.readline().strip().split()
            nspecialkpoints = int(firstline[firstline.index("special") - 1])
            nkpoints =  int(firstline[firstline.index("k-points,") - 1])
            nbands = int(firstline[-2])
        
        kpoints = np.zeros((nkpoints, 3))
        eigenvals = np.zeros((nbands, nkpoints))
        
        nelectron = 0
        for band in range(nbands):
            line_num = 2 + nspecialkpoints + 2 + band
            occ = float(linecache.getline(self._Band_File, line_num).split()[2])
            if occ:
                nelectron += int(occ)
            else:
                cb_ndx = band
                vb_ndx = band - 1
                break
        self.nelectron = nelectron
        
        for kpt in range(nkpoints):
            line_num = nspecialkpoints + 2 * (kpt + 1) + nbands * kpt
            kpos = (linecache.getline(self._Band_File, line_num)).split()[5:8]
            kpoints[kpt][0], kpoints[kpt][1], kpoints[kpt][2] = \
                float(kpos[0]), float(kpos[1]), float(kpos[2])

            for band in range(nbands):
                line_num = 2 + nspecialkpoints + 2 * (kpt + 1) + band + nbands * kpt
                ene = float(linecache.getline(self._Band_File, line_num).split()[1])
                eigenvals[band][kpt] = ene
        
        self.kpoints = kpoints
        self.eigenvals = eigenvals
        self.vbm = eigenvals[vb_ndx,:].max()
        self.cbm = eigenvals[cb_ndx,:].min()
        self.efermi = (self.cbm + self.vbm)/2  # set fermi level in the gap middle
        self.bandstructure = BandStructure(
            self.kpoints,
            {Spin(1):self.eigenvals},
            self.lattice_rec,
            self.efermi,
            structure = self.structure
        )
        
    def set_fermi(self, efermi):
        '''
        efermi: accurate fermi level
        
        '''
        self.efermi = efermi
        self.bandstructure.efermi = efermi

    def get_band_structure(self, Structure_File = 'cp2k.inp', Band_File = 'cp2k.bs'):
        self.parse_structure(Structure_File)
        self.parse_bandstructure(Band_File)
        return self.bandstructure
        
    def get_core_level(self):
        #return self.eigenvals.min()
        return np.average(self.eigenvals[:9,:])
