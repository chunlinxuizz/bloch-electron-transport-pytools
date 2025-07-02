import re
import h5py
import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

str_to_spin = {"up": Spin.up, "down": Spin.down}

def load_deformation_potentials(filename):
    deformation_potentials = {}
    with h5py.File(filename, "r") as f:
        deform_keys = [k for k in list(f.keys()) if "deformation_potentials" in k]
        for key in deform_keys:
            spin = str_to_spin[key.split("_")[-1]]
            deformation_potentials[spin] = np.array(f[key])

        structure_str = np.string_(np.array(f["structure"])).decode()
        structure = Structure.from_str(structure_str, fmt="json")
        kpoints = np.array(f["kpoints"])

    return deformation_potentials, kpoints, structure

def write_deformation_potentials(
    deformation_potentials, kpoints, structure, filename="deformation.h5"
):
    with h5py.File(filename, "w") as f:
        for spin, spin_deform in deformation_potentials.items():
            name = "deformation_potentials_{}".format(spin.name)
            f.create_dataset(name, data=spin_deform, compression="gzip")
        f["structure"] = np.string_(structure.to_json())
        f["kpoints"] = kpoints
    return filename

def monkhorst_pack(size, time_inv = True):
    """Construct a uniform sampling of k-space of given size."""
    if np.less_equal(size, 0).any():
        raise ValueError('Illegal size: %s' % list(size))
    kpts = np.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    kpts = (kpts + 0.5) / size - 0.5
    if time_inv:
        kpts = kpts[len(kpts)//2:]
    return kpts

def generate_deformation_potentials(nband, deformation_potential_tensor, kpoints):

    if len(deformation_potential_tensor) != nband:
        raise ValueError(
            "The number of deformation potential tensor should equal to `nband`."
        )
    
    deformation_potential_tensor = np.array(deformation_potential_tensor)
    deformation_potentials = np.zeros((nband, len(kpoints),3,3))
    
    for band in range(nband):
        if deformation_potential_tensor[band].shape[0] == 3:
            if len(deformation_potential_tensor[band].shape) == 1:
                deformation_potentials[band,:] = np.diag(deformation_potential_tensor[band])
            else:
                deformation_potentials[band,:] = np.asarray(deformation_potential_tensor[band])
                if deformation_potential_tensor[band].shape != (3, 3):
                    raise ValueError("Unsupported deformation potential tensor shape.")
        else:
            raise ValueError("Unsupported deformation potential tensor shape.")
    return deformation_potentials

def parse_structure_cp2k(
    structure_file: str = "cp2k.inp",
):
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
    with open(structure_file,'r') as f:
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
                charge = int(line.strip().split()[1])
            elif re.search("MULTIPLICITY", line):
                multiplicity = int(line.strip().split()[1])
    lattice = Lattice(lattice, pbc)
    
    return Structure(
        lattice, 
        species, 
        coords, 
        charge = charge,
        validate_proximity = False,
        to_unit_cell = False,
        coords_are_cartesian = True,
        site_properties = None
    )

def parse_structure_vasp(
    structure_file: str = "POSCAR",
):
    from pymatgen.io.vasp.inputs import Poscar
    poscar = Poscar.from_file(structure_file)
    
    return poscar.structure
    
    
def main(
    structure_file: str,
    structure_file_style: str,
    kmesh_size,
    nband: int,
    deformation_potential_tensor: list,
    filename: str = 'deformation.h5',
    spin: str = 'up',
):
    structure_parser = 'parse_structure_' + structure_file_style
    structure = eval(structure_parser + '(structure_file)')
    kpoints = monkhorst_pack(kmesh_size)
    spin = str_to_spin[spin]
    deformation_potentials = {
        spin: generate_deformation_potentials(
                            nband, 
                            deformation_potential_tensor, 
                            kpoints,
    )}
    
    write_deformation_potentials(
        deformation_potentials, 
        kpoints, 
        structure, 
        filename=filename,
    )

#----
kwargs={
"kmesh_size": [3,3,3],
"nband": 2,
"deformation_potential_tensor":[
# band 1 [Dxx,Dyy,Dzz]:
    [8.1, 0.5,2],
# band 2 [Dxx,Dyy,Dzz]:
    [81, 0.5, 2],
# band N... 
],
"spin": "up",
"structure_file": 'cp2k.inp',
"structure_file_style": 'cp2k',
"filename": 'deformation.h5',
}

if __name__ == '__main__':
    main(**kwargs)
