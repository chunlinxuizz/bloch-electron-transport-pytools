#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:15:32 2023

@author: chunlinxu
"""

import json
import numpy as np
import sys

def read_transport(
    transport_file: str,
) -> dict:

    with open(transport_file, "r") as ftrans:
        transport = json.load(ftrans)
    
    return transport

def do_transpose(vector, tensor):
    vector = np.array(vector)
    vector = vector / np.linalg.norm(vector)
    tensor = np.array(tensor)  
  
    return np.dot(np.dot(vector.T,tensor),vector)

def calc_power_factor(cond, seebeck):
    cond = np.array(cond)
    seebeck = np.array(seebeck)
    shape = np.shape(cond)
    ndop, nT = shape[0], shape[1]
    PF = np.zeros(shape)
    for i in range(ndop):
        for j in range(nT):
            sig_nT = cond[i][j]
            S_nT = seebeck[i][j]
            PF[i][j] = (S_nT**2) * sig_nT * 1E-6 # uW/(m*K^2)
    return list(PF)


def write_transport(
    transport: dict,
    lattice_vec,
    transpose: bool = False,
):
    """
    transport_data: a list with shape (ndop, nT, 3, 3) 
    """
    
    transports = {"conductivity", "seebeck", "electronic_thermal_conductivity", "power_factor"}
    
    dop_list = transport["doping"]
    T_list = transport["temperatures"]
    
    for key in transport.keys():
        if key in transports:
            
            transport_data = transport[key]

            if key == "conductivity":
                transport_data = list(np.array(transport_data)/100) # S/m to S/cm

            style = "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\n"

            for Tndx, T in enumerate(T_list):
                outfile_name = key + "_" + str(T) + "K.dat"
                with open(outfile_name, "w") as fout:
                    for dndx, dop in enumerate(dop_list):
                        if transpose:
                            fout.write("{:.5e}\t".format(dop))
                            for a in range(np.shape(lattice_vec)[0]):
                                data = do_transpose(lattice_vec[a], transport_data[dndx][Tndx])
                                fout.write("{:.5e}\t".format(data))
                            fout.write("\n")
                        else:
                            fout.write(
                                style.format(
                                dop,
                                transport_data[dndx][Tndx][0][0], transport_data[dndx][Tndx][0][1], transport_data[dndx][Tndx][0][2],
                                transport_data[dndx][Tndx][1][0], transport_data[dndx][Tndx][1][1], transport_data[dndx][Tndx][1][2],
                                transport_data[dndx][Tndx][2][0], transport_data[dndx][Tndx][2][1], transport_data[dndx][Tndx][2][2],
                                            )
                                )
        elif key == "mobility" and transport[key]:
            mobility = transport[key]
                        
            for scattering_type in mobility.keys():
                transport_data = mobility[scattering_type]
                
                for Tndx, T in enumerate(T_list):
                    outfile_name = key + "_" + scattering_type + "_" + str(T) + "K.dat"
                    
                    with open(outfile_name, "w") as fout:
                        for dndx, dop in enumerate(dop_list):
                            if transpose:
                                fout.write("{:.5e}\t".format(dop))
                                for a in range(np.shape(lattice_vec)[0]):
                                    data = do_transpose(lattice_vec[a], transport_data[dndx][Tndx])
                                    fout.write("{:.5e}\t".format(data))
                                fout.write("\n")
                            else:
                                fout.write(
                                    style.format(
                                        dop,
                                        transport_data[dndx][Tndx][0][0], transport_data[dndx][Tndx][0][1], transport_data[dndx][Tndx][0][2],
                                        transport_data[dndx][Tndx][1][0], transport_data[dndx][Tndx][1][1], transport_data[dndx][Tndx][1][2],
                                        transport_data[dndx][Tndx][2][0], transport_data[dndx][Tndx][2][1], transport_data[dndx][Tndx][2][2],
                                                )
                                )

#---------------------

if __name__ == "__main__":
    try:
        transport_file = sys.argv[1]
    except IndexError:
        print("\nusage: write_transport_json.py [json_filename] [structure_filename](optional)\n")
        sys.exit()
        
    try:
        structure_file = sys.argv[2]
        print("\ntransports along cell lattices(abc) will be written...\n")
        if "inp" in structure_file:
            from my_cp2k_packages.cp2k_to_pymatgen import Cp2kBand
            system = Cp2kBand()
            system.parse_structure(structure_file)
            lattice = system.structure.lattice.matrix
            transpose = True
        elif "POSCAR" in structure_file:
            from pymatgen.io.vasp.inputs import Poscar
            poscar = Poscar.from_file(structure_file)
            lattice = poscar.structure.lattice.matrix
            transpose = True
        else:
            print("\nunknown structure file...\n")

    except IndexError:
        print("\nfull transport tensor will be written...\n")
        transpose = False
        lattice = []
        
    transport = read_transport(transport_file)
    transport['power_factor'] = calc_power_factor(
            transport['conductivity'],
            transport['seebeck'],
        )
    write_transport(transport, lattice, transpose)
