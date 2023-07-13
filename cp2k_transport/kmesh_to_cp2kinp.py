#!/usr/bin/env python
import sys, linecache, os, re
from pathlib import Path
import numpy as np

def monkhorst_pack(size, time_inv = True):
    """Construct a uniform sampling of k-space of given size."""
    if np.less_equal(size, 0).any():
        raise ValueError('Illegal size: %s' % list(size))
    kpts = np.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    kpts = (kpts + 0.5) / size - 0.5
    if time_inv:
        kpts = kpts[len(kpts)//2:]
    return kpts

def write_kmesh_to_cp2kinp(
    kpts: np.ndarray,
    Out_File: str,
    start_kpoint: int,
    end_kpoint: int,
    CP2k_inp_File: str = "template.inp",
) -> None:

    """
    generate CP2K .inp file from a template .inp file
    and IBZKPT file with dense k-mesh from VASP
    """

    tailstr = []

    with open(CP2k_inp_File, 'r') as finp:
        lines = finp.readlines()
        for linendx, line in enumerate(lines):
            if re.search("&KPOINT_SET", line):
                kpoint_line = linendx
                kpoint_col = line.index("&") + 2
            if re.search("&END KPOINT_SET", line):
                for l in range(linendx, len(lines)):
                    tailstr.append(lines[l])
                break

    with open(Out_File, 'w') as fout:
        for line in range(0, kpoint_line):
            fout.write(lines[line])

        fout.write(" " * (kpoint_col-2) + "&KPOINT_SET\n")
        fout.write(" " * kpoint_col + "UNITS B_VECTOR\n")

        for line in range(start_kpoint, end_kpoint):
            kpoint = kpts[line]
            space = ' ' * (kpoint_col)
            line_to_write = "{}SPECIAL_POINT  {:>16.10f}  {:>16.10f}  {:>16.10f}\n".format(
            space, kpoint[0], kpoint[1], kpoint[2]
            )
            fout.write(line_to_write)

        fout.write(" "*kpoint_col + "NPOINTS 1\n")

        for item in tailstr:
            fout.write(item)

def write_kmesh_to_many_cp2kinps(
    kpts: np.ndarray,
    Out_File: str,
    n_kpoint_per_file: int = 1000,
    CP2k_inp_File: str = "template.inp",
):
    n_tot_kpoints = len(kpts)
    if n_tot_kpoints <= n_kpoint_per_file:
        out_file = Out_File + ".inp"
        write_kmesh_to_cp2kinp(
            kpts, 
            out_file, 
            0, 
            n_tot_kpoints, 
            CP2k_inp_File,
        )
    else:
        nkpoints = [0,]
        for i in range(1, n_tot_kpoints // n_kpoint_per_file + 1):
            nkpoints.append(nkpoints[0] + n_kpoint_per_file*i)
        nkpoints.append(n_tot_kpoints)

        for i in range(len(nkpoints)-1):
            folder_name = Out_File + "_" + str(i+1)
            os.mkdir(folder_name)
            out_file = folder_name + "/" + "cp2k.inp"
            write_kmesh_to_cp2kinp(
                kpts,
                out_file, 
                nkpoints[i], 
                nkpoints[i+1], 
                CP2k_inp_File,
            )

def _main():
    print('='*60)

    k_mesh = (input("Please input the k-mesh (i.e., kx ky kz):\n")).split()
    if len(k_mesh) == 3:
        k_mesh = [int(n) for n in k_mesh]
    else:
        print("K-mesh should be three integers!")
        sys.exit(0)
    kpts = monkhorst_pack(k_mesh)
    n_kpoint_per_file = int(input(f"The number of kpoints in each calculation? There are {len(kpts)} kpoints in total\n"))
    
    file_name = input("Please input the filename of the calculations? Enter to the default one 'band'\n")
    if not file_name:
        file_name = "band"
        
    print("Start to generate inputs (y/n)? Enter to start")
    start = input()
    if not start == 'n' or start == 'N':
        write_kmesh_to_many_cp2kinps(
        kpts,
        file_name,
        n_kpoint_per_file,
        )
        np.savetxt("KMESH", kpts, fmt = '%.10f')
    else:
        sys.exit(0)
            
    print("All done!")
#------------------------------------------

if __name__ == '__main__':
    _main()
